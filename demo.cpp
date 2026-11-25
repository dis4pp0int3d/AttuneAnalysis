// SPDX-License-Identifier: MIT
#include "BPMDetector.h"
#include "KeyDetector.h"
#include <array>
#include <cmath>
#include <iostream>

// Generate four bars of a C-major progression at 120 BPM.
int main()
{
    constexpr int rate = 22050;
    constexpr double tau = 6.283185307179586;
    const std::array<std::array<int, 3>, 4> chords {{{60,64,67}, {65,69,72}, {67,71,74}, {60,64,67}}};
    std::vector<float> audio(rate * 8);
    for (std::size_t i = 0; i < audio.size(); ++i)
    {
        const double t = double(i) / rate;
        const double beatTime = std::fmod(t, 0.5);
        double chord = 0;
        for (int midi : chords[i / (2 * rate)])
        {
            const double hz = 440.0 * std::exp2((midi - 69) / 12.0);
            chord += std::sin(tau * hz * t);
        }
        const double click = std::sin(tau * 1500.0 * beatTime) * std::exp(-150.0 * beatTime);
        audio[i] = static_cast<float>(0.2 * chord * std::exp(-4.0 * beatTime) + 0.3 * click);
    }

    const auto tempo = BPMDetector::analyze(audio, rate);
    const auto key = KeyDetector::configured(KeyDetector::defaultStyle).detectKey(audio, rate);
    std::cout << tempo.bpm << " BPM, " << (key ? key->fullName : "unknown key") << '\n';
    return std::abs(tempo.bpm - 120.0) < 1.0 && key && key->shortName == "C" ? 0 : 1;
}
