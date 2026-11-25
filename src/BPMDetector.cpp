// SPDX-License-Identifier: MIT
#include "BPMDetector.h"

#include "MagnitudeFFT.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

namespace
{
    constexpr double pi = 3.14159265358979323846;

    constexpr int envOrder = 10;             // 1024-point FFT
    constexpr int envSize = 1 << envOrder;
    constexpr int envHop = 256;              // 86.13 envelope frames per second at 22.05 kHz
    constexpr double maxSeconds = 120.0;     // analyse at most the middle two minutes
    constexpr double minBpm = 60.0, maxBpm = 200.0;
    constexpr double priorCentre = 120.0;    // BPM
    constexpr double priorWidth = 0.75;      // octaves (standard deviation)

    std::vector<double> onsetEnvelope (const float* samples, size_t count, const AnalysisCancellation& cancelled)
    {
        std::vector<double> flux;
        if (count <= size_t (envSize))
            return flux;

        std::vector<float> window (static_cast<size_t> (envSize));
        for (int n = 0; n < envSize; ++n)
            window[size_t (n)] = float (0.5 - 0.5 * std::cos (2.0 * pi * n / envSize));

        attune_detail::MagnitudeFFT fft (envOrder);
        std::vector<float> buffer (size_t (envSize) * 2);
        std::vector<double> previous (size_t (envSize / 2), 0.0), current (previous.size());
        flux.reserve (count / envHop);
        bool first = true;
        for (size_t start = 0; start + envSize <= count; start += envHop)
        {
            if ((flux.size() % 256) == 0 && cancelled)
                return {};
            std::fill (buffer.begin(), buffer.end(), 0.0f);
            for (int n = 0; n < envSize; ++n)
                buffer[size_t (n)] = samples[start + size_t (n)] * window[size_t (n)];
            fft.transform (buffer.data());
            double sum = 0.0;
            for (size_t bin = 0; bin < current.size(); ++bin)
            {
                current[bin] = std::log1p (1000.0 * double (buffer[bin]) / double (envSize / 2));
                sum += std::max (0.0, current[bin] - previous[bin]);
            }
            flux.push_back (first ? 0.0 : sum); // the first frame has nothing to compare with
            first = false;
            std::swap (previous, current);
        }

        // Subtract a +-0.25 s moving average so sustained loudness does not count as rhythm.
        const int n = int (flux.size());
        const int radius = 21;
        std::vector<double> prefix (size_t (n) + 1, 0.0), envelope (static_cast<size_t> (n), 0.0);
        for (int i = 0; i < n; ++i)
            prefix[size_t (i) + 1] = prefix[size_t (i)] + flux[size_t (i)];
        for (int i = 0; i < n; ++i)
        {
            const int lo = std::max (0, i - radius), hi = std::min (n - 1, i + radius);
            const double mean = (prefix[size_t (hi) + 1] - prefix[size_t (lo)]) / double (hi - lo + 1);
            envelope[size_t (i)] = std::max (0.0, flux[size_t (i)] - mean);
        }
        return envelope;
    }

    double interpolate (const std::vector<double>& values, double x)
    {
        if (x < 0.0 || values.size() < 2 || x > double (values.size() - 1))
            return 0.0;
        const auto i = std::min (size_t (x), values.size() - 2);
        const double f = x - double (i);
        return values[i] * (1.0 - f) + values[i + 1] * f;
    }

    // Magnitude spectrum of the (Hann-tapered, zero-padded) envelope, read at any tempo.
    class PulseSpectrum
    {
    public:
        PulseSpectrum (const std::vector<double>& envelope, double framesPerSecond)
        {
            int order = 12;
            while ((size_t (1) << order) < envelope.size() * 8 && order < 20)
                ++order;
            const size_t size = size_t (1) << order;
            std::vector<float> buffer (size * 2, 0.0f);
            const double n = double (envelope.size());
            for (size_t i = 0; i < envelope.size() && i < size; ++i)
                buffer[i] = float (envelope[i] * (0.5 - 0.5 * std::cos (2.0 * pi * (double (i) + 0.5) / n)));
            attune_detail::MagnitudeFFT fft (order);
            fft.transform (buffer.data());
            magnitudes.assign (buffer.begin(), buffer.begin() + std::ptrdiff_t (size / 2));
            binsPerHz = double (size) / framesPerSecond;
        }

        double at (double hz) const { return interpolate (magnitudes, hz * binsPerHz); }

        // Pulse strength at a tempo: the strongest of the beat, eighth and sixteenth rates. Hats
        // often outweigh kicks in the envelope, so the beat rate alone is not dependable, and
        // taking the strongest level keeps this term neutral about the octave.
        double pulse (double bpm) const
        {
            const double hz = bpm / 60.0;
            return std::max ({ at (hz), at (2.0 * hz), at (4.0 * hz) });
        }

    private:
        std::vector<double> magnitudes;
        double binsPerHz = 1.0;
    };

    double prior (double bpm)
    {
        const double octaves = std::log2 (bpm / priorCentre) / priorWidth;
        return std::exp (-0.5 * octaves * octaves);
    }

    // True when two tempi are the same pulse counted at another level (x2, x3, x1.5 ...).
    bool metricallyRelated (double a, double b)
    {
        const double ratio = std::max (a, b) / std::min (a, b);
        for (const double r : { 1.0, 1.5, 2.0, 3.0, 4.0 })
            if (std::abs (ratio - r) / r < 0.04)
                return true;
        return false;
    }
}

TempoEstimate BPMDetector::analyze (const std::vector<float>& audio, double sampleRate,
                                    const AnalysisCancellation& cancelled)
{
    TempoEstimate result;
    if (audio.empty() || cancelled || ! std::isfinite (sampleRate) || sampleRate < 4000.0 || sampleRate > 192000.0)
        return result;

    // Intros and outros often have no drums: use the middle of long files.
    const auto maxSamples = size_t (maxSeconds * sampleRate);
    const size_t count = std::min (audio.size(), maxSamples);
    const size_t offset = (audio.size() - count) / 2;
    const auto envelope = onsetEnvelope (audio.data() + offset, count, cancelled);
    const int n = int (envelope.size());
    const double fps = sampleRate / envHop;
    const double seconds = double (n) / fps;
    // At least a couple of beats at the slowest tempo are needed for a pulse to exist.
    if (cancelled || seconds < 1.5 || *std::max_element (envelope.begin(), envelope.end()) <= 0.0)
        return result;
    {
        auto sorted = envelope;
        const auto top = sorted.size() / 10;
        std::nth_element (sorted.begin(), sorted.begin() + std::ptrdiff_t (top), sorted.end(), std::greater<double>());
        const double total = std::accumulate (sorted.begin(), sorted.end(), 0.0);
        const double peaks = std::accumulate (sorted.begin(), sorted.begin() + std::ptrdiff_t (top), 0.0);
        result.onsetSpikiness = total > 0.0 ? peaks / total : 0.0;
    }

    // Autocorrelation up to four beats at the slowest tempo (or half the excerpt).
    const int maxLag = std::min (n / 2, int (std::ceil (4.0 * 60.0 / minBpm * fps)) + 2);
    std::vector<double> ac (size_t (maxLag) + 1, 0.0);
    for (int lag = 0; lag <= maxLag; ++lag)
    {
        if (cancelled)
            return {};
        double dot = 0.0;
        for (int i = 0; i < n - lag; ++i)
            dot += envelope[size_t (i)] * envelope[size_t (i + lag)];
        ac[size_t (lag)] = dot / double (n - lag);
    }
    if (! (ac[0] > 0.0))
        return result;
    for (int lag = maxLag; lag >= 0; --lag)
        ac[size_t (lag)] /= ac[0];

    const PulseSpectrum spectrum (envelope, fps);
    if (cancelled)
        return {};

    // Repetition one to four beats apart; lags beyond the excerpt are simply left out.
    const auto repetition = [&ac, fps, maxLag] (double bpm)
    {
        const double beat = 60.0 * fps / bpm;
        double sum = 0.0;
        int used = 0;
        for (int k = 1; k <= 4; ++k)
        {
            const double lag = k * beat;
            if (lag > double (maxLag))
                break;
            // A peak may sit between two integer lags; take the best of a +-1 frame window.
            double best = 0.0;
            for (const double shift : { -1.0, -0.5, 0.0, 0.5, 1.0 })
                best = std::max (best, interpolate (ac, lag + shift));
            sum += best;
            ++used;
        }
        return used > 0 ? sum / used : 0.0;
    };
    const auto evidence = [&] (double bpm) { return spectrum.pulse (bpm) * repetition (bpm); };

    const int steps = int (std::round ((maxBpm - minBpm) * 10.0));
    std::vector<double> scores (size_t (steps) + 1, 0.0);
    int best = -1;
    for (int i = 0; i <= steps; ++i)
    {
        const double bpm = minBpm + i * 0.1;
        scores[size_t (i)] = evidence (bpm) * prior (bpm);
        if (best < 0 || scores[size_t (i)] > scores[size_t (best)])
            best = i;
    }
    if (best < 0 || ! (scores[size_t (best)] > 0.0))
        return result;
    double bpm = minBpm + best * 0.1;

    // Refine: the envelope spectrum resolves tempo to a fraction of 60/seconds BPM. Stay within
    // the width of that peak so the chosen octave and metre are kept.
    {
        const double span = std::clamp (45.0 / seconds, 0.15, 0.03 * bpm);
        const double step = span / 200.0;
        double fine = bpm, fineScore = -1.0;
        for (double candidate = bpm - span; candidate <= bpm + span; candidate += step)
        {
            const double score = spectrum.pulse (candidate);
            if (score > fineScore)
            {
                fineScore = score;
                fine = candidate;
            }
        }
        bpm = fine;
    }

    // Strongest reading that is not the same tempo; used for confidence and as the alternative.
    double rival = 0.0, rivalScore = 0.0, unrelatedScore = 0.0;
    for (int i = 0; i <= steps; ++i)
    {
        const double other = minBpm + i * 0.1;
        if (std::abs (other - bpm) / bpm < 0.06)
            continue;
        if (scores[size_t (i)] > rivalScore)
        {
            rivalScore = scores[size_t (i)];
            rival = other;
        }
        if (! metricallyRelated (other, bpm))
            unrelatedScore = std::max (unrelatedScore, scores[size_t (i)]);
    }
    const double top = scores[size_t (best)];
    // Mostly "is there a clear pulse at all" (unrelated tempi are weak), partly "is the octave clear".
    const double clarity = 1.0 - unrelatedScore / top;
    const double octave = 1.0 - rivalScore / top;
    result.bpm = bpm;
    result.confidence = std::clamp (0.75 * clarity + 0.25 * octave, 0.0, 1.0);
    result.alternativeBpm = (rival > 0.0 && rivalScore > 0.5 * top) ? rival
                          : (bpm >= 100.0 ? bpm * 0.5 : bpm * 2.0);
    return result;
}
