#pragma once
#include "gtest/gtest.h"

// GTEST_SKIP arrived in Google Test 1.10; the VS template ships 1.8.1.
// Fallback: report the reason and return. The test shows PASSED, not SKIPPED.
#ifndef GTEST_SKIP
#define GTEST_SKIP() \
    return GTEST_MESSAGE_("Skipped: ", ::testing::TestPartResult::kSuccess)
#endif