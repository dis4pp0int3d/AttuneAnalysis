// SPDX-License-Identifier: MIT
#pragma once

#include "AnalysisCancellation.h"

#include <vector>

// Tempo of material without attacks (pads, strings, drones with chord changes).
//
// A pad has no onsets to count, but its chords change on bar lines. This finds how often the
// harmony changes (a chroma novelty curve, then its periodicity) and reads that period as one,
// two or four 4/4 bars. Which of those it is cannot be heard, so the tempo is reported in
// foldLow..2*foldLow; the half/double alternative is as plausible.
struct HarmonicTempoEstimate
{
    double bpm = 0.0;           // 0 when the harmony does not change regularly
    double changePeriod = 0.0;  // seconds between chord changes
    double confidence = 0.0;    // 0..1 heuristic: how regular the changes are
};

namespace HarmonicTempo
{
    constexpr double foldLow = 85.0; // reported tempi lie in [85, 170)

    // Mono audio at 22050 Hz (22.05 kHz).
    // tuningCents: where the music sits relative to A4 = 440 Hz (from the key detector).
    HarmonicTempoEstimate analyze (const std::vector<float>& mono22050, double sampleRate,
                                   const AnalysisCancellation& cancelled = {}, double tuningCents = 0.0);

    // Moves a tempo by octaves into [foldLow, 2 * foldLow).
    double fold (double bpm);
}
