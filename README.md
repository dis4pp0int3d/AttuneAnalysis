# Attune Analysis

A small C++17 library for estimating tempo and musical key from mono audio. It also reports a Camelot code and tuning offset in cents. It builds with the C++ standard library alone.

## Try it

From this directory, run these commands with CMake 3.20 or later and a C++17 compiler.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release -V
```

The demo generates four bars in C major at 120 BPM and prints the detected tempo and key.

## Use it in your project

Add this directory with `add_subdirectory` and link your target to `attune_analysis`. Decode your audio, mix it to mono and resample it to 22,050 Hz. Pass the samples as a `std::vector<float>`, typically in the range −1 to +1.

```cpp
#include "BPMDetector.h"
#include "KeyDetector.h"

// mono holds the decoded samples at 22,050 Hz.
auto tempo = BPMDetector::analyze(mono, 22050);
auto key = KeyDetector::configured(KeyDetector::defaultStyle).detectKey(mono, 22050);
```

Read `tempo.bpm` for the main estimate and `tempo.alternativeBpm` for another possible reading. A BPM of zero means no tempo was found. For pads with regular chord changes, try `HarmonicTempo::analyze` from `HarmonicTempo.h`.

When `key` has a value, read `key->fullName`, `key->camelot` and `key->tuningCents`. Short clips and half-time or double-time rhythms can give ambiguous results. The full Attune plug-in also uses loop length and file names to refine its estimates.

The code uses the MIT license. KeyFinder model weights and tonal profiles use CC BY-SA 4.0. [Credits and license details](ATTRIBUTION.md)
