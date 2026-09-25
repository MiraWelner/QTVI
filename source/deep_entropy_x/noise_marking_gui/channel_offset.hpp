#pragma once
/**
 * @file   channel_output.hpp
 * @brief  ChannelOutput: one channel's view of one bin's partition.
 *
 */

#include <cstdint>
#include <vector>

#include "template_bank.hpp"
#include "seed_pool.hpp"
#include "pvc_filter.hpp"

namespace tbank {

    struct ChannelOutput {
        TemplateBank               bank;
        std::vector<int32_t>       assignment;   // per beat
        seed_pool::SeedSelection   seed;
        std::vector<CapRaiseEvent> cap_raises;
        BinCounts                  counts;
        std::vector<BeatFlags>     flags;
        pvc_filter::FilterResult   pvc;
    };

}  // namespace tbank