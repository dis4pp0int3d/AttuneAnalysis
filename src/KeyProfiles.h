// SPDX-License-Identifier: CC-BY-SA-4.0
// KeyFinder learned50 profiles. See ATTRIBUTION.md.
#pragma once
#include <array>

namespace keyfinder_profiles
{
    // Identical to KeyDetector.swift: development-only GiantSteps+ v2 profiles,
    // 50% empirical median / 50% normalized Krumhansl-Schmuckler. Provenance:
    // docs/briefs/0004-public-key-challenge/profile-provenance.json.
    constexpr std::array<double, 12> majorProfile {
        0.1275930110697475, 0.0623281238962653, 0.0856139951706759, 0.062434445057663025,
        0.10059235106369269, 0.08946810513660708, 0.06685370853481933, 0.11450725102532114,
        0.06293247226287722, 0.08759937643779363, 0.061632033039242085, 0.07844512730529522 };
    constexpr std::array<double, 12> minorProfile {
        0.12295831518055632, 0.06713257770803258, 0.08267739330362738, 0.10456108651196024,
        0.06793500512433348, 0.08110914848416328, 0.06509663723126663, 0.10267883579369759,
        0.08369819182992237, 0.06756749780032058, 0.0806585670113963, 0.07392674402072341 };
}
