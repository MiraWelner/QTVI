#pragma once
#pragma once
/**
 * @file   data_bin_header.hpp
 * @brief  The 592-byte header of a file_to_bin v1 data .bin, and where its
 *         trailing sleep-stage block starts. Read only.
 *
 * Layout (the same one anneal_handler.cpp read_data_bin reads):
 *   uint32 version (= 1), uint32 n_channels (= 36), uint32 sleep_state_len,
 *   uint32 sizes_up[36], uint32 sizes_raw[36], float native_rates[36],
 *   float up_rates[36], uint32 sleep_size
 * then per channel { sizes_up doubles, 2 * sizes_raw doubles }, then
 * sleep_size doubles of sleep stage.
 *
 * Its own file because two unrelated readers need different fields of it:
 * beat_times needs the ECG1 upsample rate and length, record_sleep needs the
 * epoch length and the stage block. Neither should have to know the other.
 */

#include <array>
#include <cstdint>
#include <fstream>
#include <string>

namespace data_bin_header {

    inline constexpr int kNumChannels = 36;
    inline constexpr int kSlotEcg1 = 1;
    inline constexpr uint64_t kHeaderBytes = (4 + 4 * kNumChannels) * 4;   // 592

    struct Header {
        bool ok = false;   // read, and version / channel count understood
        uint32_t sleep_state_len = 0;
        std::array<uint32_t, kNumChannels> sizes_up{};
        std::array<uint32_t, kNumChannels> sizes_raw{};
        std::array<float, kNumChannels> native_rates{};
        std::array<float, kNumChannels> up_rates{};
        uint32_t sleep_size = 0;
        uint64_t file_size = 0;

        // The annealer's primarySR: it cuts bins from the UPSAMPLED ECG1.
        double ecgRateHz() const { return static_cast<double>(up_rates[kSlotEcg1]); }
        uint64_t ecgLength() const { return sizes_up[kSlotEcg1]; }

        uint64_t sleepBlockOffset() const {
            uint64_t off = kHeaderBytes;
            for (int i = 0; i < kNumChannels; ++i)
                off += 8ull * sizes_up[i] + 16ull * sizes_raw[i];
            return off;
        }
    };

    inline Header read(const std::string& path) {
        Header h;
        std::ifstream f(path, std::ios::binary);
        if (!f) return h;
        f.seekg(0, std::ios::end);
        h.file_size = static_cast<uint64_t>(f.tellg());
        f.seekg(0, std::ios::beg);
        if (h.file_size < kHeaderBytes) return h;

        uint32_t version = 0, nCh = 0;
        f.read(reinterpret_cast<char*>(&version), 4);
        f.read(reinterpret_cast<char*>(&nCh), 4);
        // Same refusal as the annealer: a file it would not read is not read here.
        if (version != 1 || nCh != static_cast<uint32_t>(kNumChannels)) return h;
        f.read(reinterpret_cast<char*>(&h.sleep_state_len), 4);
        f.read(reinterpret_cast<char*>(h.sizes_up.data()), kNumChannels * 4);
        f.read(reinterpret_cast<char*>(h.sizes_raw.data()), kNumChannels * 4);
        f.read(reinterpret_cast<char*>(h.native_rates.data()), kNumChannels * 4);
        f.read(reinterpret_cast<char*>(h.up_rates.data()), kNumChannels * 4);
        f.read(reinterpret_cast<char*>(&h.sleep_size), 4);
        h.ok = static_cast<bool>(f);
        return h;
    }

}   // namespace data_bin_header