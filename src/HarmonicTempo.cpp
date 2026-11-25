// SPDX-License-Identifier: MIT
#include "HarmonicTempo.h"

#include "MagnitudeFFT.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>

namespace
{
    constexpr double pi = 3.14159265358979323846;
    constexpr int fftOrder = 13;               // 8192 points: 2.7 Hz bins at 22.05 kHz
    constexpr int fftSize = 1 << fftOrder;
    constexpr int hop = 1024;                  // 21.5 chroma frames per second
    constexpr double minHz = 70.0, maxHz = 2000.0;
    constexpr double contextSeconds = 0.45;    // harmony compared over this much before / after
    constexpr double minPeriod = 0.9, maxPeriod = 8.0; // seconds per chord

    using Chroma = std::array<double, 12>;

    std::vector<Chroma> chromagram (const std::vector<float>& audio, double rate, const AnalysisCancellation& cancelled, double tuningCents)
    {
        std::vector<Chroma> frames;
        if (audio.size() < size_t (fftSize))
            return frames;

        std::vector<float> window (static_cast<size_t> (fftSize));
        for (int n = 0; n < fftSize; ++n)
            window[size_t (n)] = float (0.5 - 0.5 * std::cos (2.0 * pi * n / fftSize));

        // Bin -> pitch class around the music's own tuning, so detuned material (A4 = 450 Hz)
        // does not smear every note across two classes.
        const int firstBin = std::max (1, int (std::ceil (minHz * fftSize / rate)));
        const int lastBin = std::min (fftSize / 2 - 1, int (std::floor (maxHz * fftSize / rate)));
        std::vector<int> pitchClass (static_cast<size_t> (lastBin + 1), 0);
        for (int bin = firstBin; bin <= lastBin; ++bin)
        {
            const double midi = 69.0 + 12.0 * std::log2 (bin * rate / fftSize / 440.0) - tuningCents / 100.0;
            pitchClass[size_t (bin)] = ((int (std::lround (midi)) % 12) + 12) % 12;
        }

        attune_detail::MagnitudeFFT fft (fftOrder);
        std::vector<float> buffer (size_t (fftSize) * 2);
        for (size_t start = 0; start + fftSize <= audio.size(); start += hop)
        {
            if ((frames.size() % 64) == 0 && cancelled)
                return {};
            std::fill (buffer.begin(), buffer.end(), 0.0f);
            for (int n = 0; n < fftSize; ++n)
                buffer[size_t (n)] = audio[start + size_t (n)] * window[size_t (n)];
            fft.transform (buffer.data());

            Chroma chroma {};
            for (int bin = firstBin; bin <= lastBin; ++bin)
                chroma[size_t (pitchClass[size_t (bin)])] += std::log1p (100.0 * double (buffer[size_t (bin)]) / (fftSize / 2));
            double norm = 0.0;
            for (const double v : chroma)
                norm += v * v;
            norm = std::sqrt (norm);
            if (norm > 1.0e-9)
                for (auto& v : chroma)
                    v /= norm;
            frames.push_back (chroma);
        }
        return frames;
    }

    // How different the harmony after each frame is from the harmony before it (0 = same).
    std::vector<double> harmonicNovelty (const std::vector<Chroma>& frames, double framesPerSecond)
    {
        const int n = int (frames.size());
        const int w = std::max (2, int (std::lround (contextSeconds * framesPerSecond)));
        std::vector<double> novelty (static_cast<size_t> (n), 0.0);
        if (n < 2 * w + 2)
            return novelty;

        std::vector<Chroma> prefix (static_cast<size_t> (n) + 1, Chroma {});
        for (int i = 0; i < n; ++i)
            for (size_t c = 0; c < 12; ++c)
                prefix[size_t (i) + 1][c] = prefix[size_t (i)][c] + frames[size_t (i)][c];

        for (int i = w; i + w <= n; ++i)
        {
            double dot = 0.0, a2 = 0.0, b2 = 0.0;
            for (size_t c = 0; c < 12; ++c)
            {
                const double a = prefix[size_t (i)][c] - prefix[size_t (i - w)][c];
                const double b = prefix[size_t (i + w)][c] - prefix[size_t (i)][c];
                dot += a * b;
                a2 += a * a;
                b2 += b * b;
            }
            novelty[size_t (i)] = (a2 > 0.0 && b2 > 0.0) ? 1.0 - dot / std::sqrt (a2 * b2) : 0.0;
        }

        // Keep only the peaks above the local level: those are the chord changes. A peak must
        // also be a real change of notes: a held chord only wobbles (beating, vibrato, filter
        // movement) by a tiny fraction of that, and its wobble must not read as changes.
        constexpr double minChange = 0.02;
        const int radius = std::max (3, int (std::lround (0.6 * framesPerSecond)));
        std::vector<double> peaks (novelty.size(), 0.0);
        for (int i = 0; i < n; ++i)
        {
            const int lo = std::max (0, i - radius), hi = std::min (n - 1, i + radius);
            double localMax = 0.0, sum = 0.0;
            for (int j = lo; j <= hi; ++j)
            {
                localMax = std::max (localMax, novelty[size_t (j)]);
                sum += novelty[size_t (j)];
            }
            const double mean = sum / double (hi - lo + 1);
            if (novelty[size_t (i)] >= localMax && novelty[size_t (i)] > mean && novelty[size_t (i)] >= minChange)
                peaks[size_t (i)] = novelty[size_t (i)] - mean;
        }
        // Blur each peak over a few frames so slightly irregular changes still line up.
        std::vector<double> blurred (peaks.size(), 0.0);
        for (int i = 0; i < n; ++i)
            for (int d = -2; d <= 2; ++d)
                if (i + d >= 0 && i + d < n)
                    blurred[size_t (i + d)] += peaks[size_t (i)] * (1.0 - std::abs (d) / 3.0);
        return blurred;
    }
}

double HarmonicTempo::fold (double bpm)
{
    if (! (bpm > 0.0))
        return 0.0;
    while (bpm < foldLow)
        bpm *= 2.0;
    while (bpm >= 2.0 * foldLow)
        bpm *= 0.5;
    return bpm;
}

HarmonicTempoEstimate HarmonicTempo::analyze (const std::vector<float>& audio, double rate,
                                              const AnalysisCancellation& cancelled, double tuningCents)
{
    HarmonicTempoEstimate result;
    if (! (rate > 0.0))
        return result;
    const auto frames = chromagram (audio, rate, cancelled, tuningCents);
    if (cancelled)
        return {};
    const double fps = rate / hop;
    const auto novelty = harmonicNovelty (frames, fps);
    const int n = int (novelty.size());
    const double seconds = n / fps;
    // At least three changes are needed before a regular interval means anything.
    if (seconds < 2.5 * minPeriod)
        return result;

    // Autocorrelation of the change curve over chord periods that fit the excerpt twice.
    const int minLag = int (std::floor (minPeriod * fps));
    const int maxLag = std::min (int (std::ceil (maxPeriod * fps)), n / 2);
    if (maxLag <= minLag + 2)
        return result;
    double energy = 0.0;
    for (const double v : novelty)
        energy += v * v;
    if (energy <= 1.0e-12)
        return result;
    std::vector<double> ac (static_cast<size_t> (maxLag) + 2, 0.0);
    for (int lag = 1; lag <= maxLag + 1 && lag < n; ++lag)
    {
        double dot = 0.0;
        for (int i = 0; i + lag < n; ++i)
            dot += novelty[size_t (i)] * novelty[size_t (i + lag)];
        ac[size_t (lag)] = dot / energy * double (n) / double (n - lag); // unbiased, 1 = perfect repeat
    }

    // A regular change every P seconds also repeats at 2P, 3P...: score a period by its
    // multiples and keep the shortest period that explains the curve nearly as well.
    const auto combScore = [&] (int lag)
    {
        double sum = 0.0;
        int used = 0;
        for (int k = 1; k * lag <= maxLag; ++k, ++used)
            sum += std::max ({ ac[size_t (k * lag - 1)], ac[size_t (k * lag)], ac[size_t (k * lag + 1)] });
        return used > 0 ? sum / used : 0.0;
    };
    int bestLag = -1;
    double bestScore = 0.0;
    std::vector<double> scores (static_cast<size_t> (maxLag) + 1, 0.0);
    for (int lag = minLag; lag <= maxLag; ++lag)
    {
        scores[size_t (lag)] = combScore (lag);
        if (scores[size_t (lag)] > bestScore)
        {
            bestScore = scores[size_t (lag)];
            bestLag = lag;
        }
    }
    if (bestLag < 0 || bestScore < 0.15)
        return result;
    for (int lag = minLag; lag < bestLag; ++lag)
        if (scores[size_t (lag)] >= 0.85 * bestScore && ac[size_t (lag)] >= ac[size_t (lag - 1)] && ac[size_t (lag)] >= ac[size_t (lag + 1)])
        {
            bestLag = lag;
            break;
        }

    // Sub-frame period: parabola through the autocorrelation peak, then averaged over the
    // multiples that are still inside the excerpt (each multiple k pins the period to 1/k).
    double periodFrames = bestLag;
    {
        double weighted = 0.0, weights = 0.0;
        for (int k = 1; k * bestLag + 1 <= maxLag; ++k)
        {
            int peak = k * bestLag;
            for (const int d : { -1, 1 })
                if (ac[size_t (k * bestLag + d)] > ac[size_t (peak)])
                    peak = k * bestLag + d;
            if (peak <= 1 || peak + 1 > maxLag + 1)
                continue;
            const double l = ac[size_t (peak - 1)], c = ac[size_t (peak)], r = ac[size_t (peak + 1)];
            const double denominator = l - 2.0 * c + r;
            const double offset = denominator < -1.0e-12 ? std::clamp (0.5 * (l - r) / denominator, -0.5, 0.5) : 0.0;
            weighted += (peak + offset) * k * c; // the k-th multiple carries k times the precision
            weights += double (k * k) * c;
        }
        if (weights > 0.0)
            periodFrames = weighted / weights;
    }

    result.changePeriod = periodFrames / fps;
    result.bpm = fold (240.0 / result.changePeriod);   // one change per bar (or 2, 4 bars)
    result.confidence = std::clamp ((bestScore - 0.15) / 0.5, 0.0, 1.0);
    return result;
}
