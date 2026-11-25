// SPDX-License-Identifier: MIT
#pragma once
#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace attune_detail
{
inline constexpr double pi = 3.14159265358979323846;

// Radix-2 forward FFT. Each detector reuses one plan across its analysis frames.
class MagnitudeFFT
{
public:
    explicit MagnitudeFFT(int order)
        : size(std::size_t{1} << order), reversed(size), roots(size / 2), work(size)
    {
        for (std::size_t i = 1; i < size; ++i)
            reversed[i] = (reversed[i / 2] >> 1) | ((i & 1) * (size / 2));
        for (std::size_t i = 0; i < roots.size(); ++i)
            roots[i] = std::polar(1.0, -2.0 * pi * double(i) / double(size));
    }

    // Replace the first N/2+1 samples with unnormalised magnitudes from DC to Nyquist.
    void transform(float* samples)
    {
        for (std::size_t i = 0; i < size; ++i)
            work[reversed[i]] = samples[i];
        for (std::size_t width = 2; width <= size; width *= 2)
        {
            const auto half = width / 2;
            const auto stride = size / width;
            for (std::size_t base = 0; base < size; base += width)
                for (std::size_t j = 0; j < half; ++j)
                {
                    const auto even = work[base + j];
                    const auto odd = roots[j * stride] * work[base + j + half];
                    work[base + j] = even + odd;
                    work[base + j + half] = even - odd;
                }
        }
        for (std::size_t i = 0; i <= size / 2; ++i)
            samples[i] = static_cast<float>(std::abs(work[i]));
    }

private:
    std::size_t size;
    std::vector<std::size_t> reversed;
    std::vector<std::complex<double>> roots, work;
};
}
