# Credits and licenses

The tempo detectors, standalone FFT, demo and Attune code changes are distributed under the [MIT license](LICENSE).

## Key detection

IFEELVOID developed the original key detector and model for [KeyFinder](https://github.com/ifeelvoid/keyfinder). `KeyDetector.cpp` and `KeyDetector.h` include code from KeyFinderVST. Attune adapted tuning estimation, confidence scoring and analysis cancellation. This standalone edition also replaces the FFT implementation. [Original MIT notice](licenses/KeyFinder-MIT.txt)

## Model weights and tonal profiles

The numerical weights in `KeyModelWeights.h` and the tonal profiles in `KeyProfiles.h` are distributed under [CC BY-SA 4.0 International](https://creativecommons.org/licenses/by-sa/4.0/). The weights match the KeyFinder C++ release. The tonal profiles were moved to a separate header with their values preserved.

The weights derive from annotations in the datasets below and track labels from the IFEELVOID catalog. The archive contains the derived numerical parameters. [KeyFinder model provenance](https://github.com/ifeelvoid/keyfinder/blob/main/docs/briefs/0010-genre-key-models/model-provenance.json)

- Ángel Faraldo and Daniel G. Camhi. [GiantSteps+ EDM Key Dataset, v2](https://doi.org/10.5281/zenodo.4153506).
- Eduard Mas Marín and Ángel Faraldo. [GiantSteps-MTG Key Dataset](https://doi.org/10.5281/zenodo.1101082).

When sharing these data, retain the credits, source links and license. Identify any changes and share adapted data under the same terms.
