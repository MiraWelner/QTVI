#pragma once
/**
 * @file   record_sleep.hpp
 * @brief  A record's sleep staging, and the stage at any time in it.
 *
 * From the ORIGINAL data .bin, not the annealed one. The annealed file keeps
 * only the epochs whose END falls inside a good stretch, as a bare list with
 * no epoch numbers, so a time near a splice can sit in an epoch that list
 * dropped. The original holds every epoch.
 *
 * Epoch e covers [e * E, (e + 1) * E) seconds from the start of the record --
 * the hypnogram's convention. Codes 0..4 are Wake, NREM1, NREM2, NREM3, REM.
 * stageAt returns -1 for anything else (unscored, movement, artifact, NaN, a
 * non-integer code) and for a time before the record or past the last epoch.
 *
 * Nothing here knows about beats or templates. Read only.
 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "annealing/data_bin_header.hpp"

namespace record_sleep {

    inline constexpr int kNumStages = 5;   // Wake, NREM1, NREM2, NREM3, REM
    inline constexpr const char* kStageColumn[kNumStages] = {
        "pct_wake", "pct_nrem1", "pct_nrem2", "pct_nrem3", "pct_rem"
    };

    struct RecordSleep {
        double epoch_sec = 0.0;
        std::vector<double> stages;   // one per epoch, record order

        // file_to_bin writes a single -1.0 for "no staging" (chart_utils.hpp,
        // sleep_data_present). Demand at least one real code, so Bittium and
        // CHAOS read as absent rather than as a record that is never asleep.
        bool present() const {
            if (!(epoch_sec > 0.0) || stages.empty()) return false;
            for (double v : stages)
                if (std::isfinite(v) && v >= 0.0 && v <= 4.0) return true;
            return false;
        }

        int stageAt(double t) const {
            if (!std::isfinite(t) || t < 0.0 || !(epoch_sec > 0.0)) return -1;
            const double e = std::floor(t / epoch_sec);
            if (e >= static_cast<double>(stages.size())) return -1;
            const double v = stages[static_cast<size_t>(e)];
            if (!std::isfinite(v)) return -1;
            const int code = static_cast<int>(v);
            return (code >= 0 && code <= 4 && v == static_cast<double>(code)) ? code : -1;
        }
    };

    // Empty (present() == false) when the header is unreadable or there is no
    // stage block. A truncated file yields the epochs that are there.
    inline RecordSleep read(const std::string& path, const data_bin_header::Header& h) {
        RecordSleep rs;
        if (!h.ok) return rs;
        rs.epoch_sec = static_cast<double>(h.sleep_state_len);
        const uint64_t off = h.sleepBlockOffset();
        if (h.sleep_size == 0 || off >= h.file_size) return rs;
        const uint64_t n = std::min<uint64_t>(h.sleep_size, (h.file_size - off) / 8);
        std::ifstream f(path, std::ios::binary);
        if (!f) return rs;
        rs.stages.resize(static_cast<size_t>(n));
        f.seekg(static_cast<std::streamoff>(off), std::ios::beg);
        f.read(reinterpret_cast<char*>(rs.stages.data()), static_cast<std::streamsize>(n * 8));
        if (!f) rs.stages.clear();
        return rs;
    }

}   // namespace record_sleep