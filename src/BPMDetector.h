// SPDX-License-Identifier: MIT
#pragma once

#include "AnalysisCancellation.h"

#include <vector>

struct TempoEstimate
{
    double bpm = 0.0;            // 0 when no tempo was found
    double confidence = 0.0;     // 0..1 heuristic, not a calibrated probability
    double alternativeBpm = 0.0; // strongest competing reading (usually half or double time)
    // Share of the onset envelope carried by its strongest 10 % of frames. Hits (drums,
    // plucks, chords struck together) make a spiky envelope (high share); the slow beating of
    // sustained notes in pads and drones makes a smooth wave (low share).
    double onsetSpikiness = 0.0;
};

// Tempo detector for mono audio at 22050 Hz (22.05 kHz).
//
// 1. Onset envelope: half-wave rectified log-magnitude spectral flux, local mean removed.
// 2. Every tempo from 60 to 200 BPM is scored by two independent views of that envelope:
//      - its spectrum at the beat, eighth and sixteenth rates (is there a pulse at exactly
//        this rate?), which rejects the 2/3 and 4/3 "ghost" tempi autocorrelation alone
//        cannot tell apart, and
//      - its autocorrelation one to four beats apart (does the pattern repeat every beat?),
//        which rejects tempi that are too fast.
//    A broad prior around 120 BPM settles the octave when the audio supports both.
// 3. The winner is refined against the whole envelope to a few hundredths of a BPM.
class BPMDetector
{
public:
    static TempoEstimate analyze (const std::vector<float>& mono22050, double sampleRate,
                                  const AnalysisCancellation& cancelled = {});
};
