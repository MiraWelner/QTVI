#include "getBinBreaksAndCount.hpp"

// MATLAB: getBinBreaksAndCount.m (inline version used in AnnealSegments.m)
// Computes bin_breaks as a 1D vector of 1-based inclusive bin endpoints,
// and bin_count based on remainder logic.

BinBreaksResult getBinBreaksAndCount(
    uint64_t total_len,
    uint64_t bin_size_samples,
    double ppgSR,
    double min_bin_size_mins)
{
    BinBreaksResult result;

    // MATLAB: remainder = mod(length(ppg), bin_size_samples)
    double remainder_mins = (double)(total_len % bin_size_samples) / ppgSR / 60.0;

    // MATLAB: if remainder < min_bin_size_mins -> floor, else ceil
    if (remainder_mins < min_bin_size_mins)
        result.bin_count = (int)std::floor((double)total_len / bin_size_samples);
    else
        result.bin_count = (int)std::ceil((double)total_len / bin_size_samples);

    // MATLAB: bin_breaks = (bin_size_samples + 1 : bin_size_samples : length(ppg))
    for (uint64_t b = bin_size_samples + 1; b <= total_len; b += bin_size_samples)
        result.bin_breaks.push_back(b);

    // MATLAB: pad or adjust last break
    if ((int)result.bin_breaks.size() < result.bin_count)
        result.bin_breaks.push_back(total_len);
    else if (!result.bin_breaks.empty())
        result.bin_breaks.back() = total_len;

    return result;
}
