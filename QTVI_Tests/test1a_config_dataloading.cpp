// ============================================================================
// test1a_config_dataloading.cpp -- QTVI_Tests (Google Test)
//
// REQUIREMENT: Load config.csv for all datasets. Verify each channel's rate
// pair, absent channels report zero, and the use_consensus_rpeak flag
// defaults to true when blank.
//
// FIXTURE: QTVI_Tests\data\config\config.csv -- a copy of the real config
// (mesa, bittium, chaos, shhs), frozen here, with three fixture edits:
//   original_file_path / output_folder = "input" / "output" for every dataset:
//       blank, load_config opens a folder-picker dialog and the test hangs;
//       real D:\ paths, it makes folders in the real output tree
//   use_consensus_rpeak: MESA BLANK (tests the default), BITTIUM FALSE (proves
//       the column is read -- TRUE alone cannot, it equals the default)
//
// KEEPING IT IN STEP: the fixture is frozen, so editing the REAL config never
// breaks this test. To test a newer config, copy it here again, redo the
// three edits, and update kExpected below to match.
//
// PROJECT SETTINGS this file needs (QTVI_Tests -> Properties):
//   C/C++ -> Preprocessor -> Preprocessor Definitions:
//       TESTS_DATA_DIR=R"($(ProjectDir)data)"
//   C/C++ -> General -> Additional Include Directories:
//       $(ProjectDir)..\source\deep_entropy_x   (the folder holding config_file_handling\)
//   Qt Project Settings -> Qt Modules: Widgets   (config_loader uses QFileDialog)
// ============================================================================
#include "pch.h"                                     // the template's; it includes gtest
#include "config_file_handling/config_loader.hpp"   // the code under test

// Without this definition every fixture path below fails to compile, with
// errors that do not name the cause. Stop with one that does.
#ifndef TESTS_DATA_DIR
#error "TESTS_DATA_DIR is not defined. QTVI_Tests -> Properties (All Configurations, All Platforms) -> C/C++ -> Preprocessor -> Preprocessor Definitions: add  TESTS_DATA_DIR=R\"($(ProjectDir)data)\""
#endif

#include <filesystem>
#include <fstream>
#include <map>
#include <string>

namespace fs = std::filesystem;

// ---- the expected answers ---------------------------------------------------
//
// One entry per dataset: the number load_config takes, the name it should
// report, the rate pairs config.csv gives it, and the consensus flag it should
// end up with. Any channel NOT listed here must come back as 0.
struct Expected {
    int choice;
    const char* dataset;
    std::map<std::string, std::pair<double, double>> rates;   // channel -> (raw, upsample)
    bool consensus;
};

static const Expected kExpected[] = {
    { 1, "MESA",
      { {"ecg", {256, 1000}}, {"ppg", {256, 500}}, {"eeg", {256, 1000}},
        {"eog_l", {256, 500}}, {"eog_r", {256, 500}}, {"emg", {256, 500}},
        {"pres", {32, 32}}, {"flow", {32, 32}}, {"snore", {32, 32}},
        {"thor", {32, 32}}, {"abdo", {32, 32}}, {"leg", {32, 32}},
        {"auxac", {32, 32}}, {"pos", {32, 32}}, {"therm", {32, 32}},
        {"oxstatus", {1, 1}}, {"spo2", {1, 1}}, {"hr", {1, 1}},
        {"dhr", {256, 500}} },
      true },   // blank -> must default to true
    { 2, "BITTIUM",
      { {"ecg", {500, 1000}}, {"ppg", {0, 500}}, {"accel", {25, 25}},
        {"temp", {1, 1}}, {"marker", {1, 1}}, {"pacemaker", {8, 8}} },
      false },   // "FALSE" -> false
    { 3, "CHAOS",
      { {"ecg", {500, 1000}}, {"ppg", {125, 500}}, {"cvp", {125, 125}},
        {"abp", {125, 500}}, {"art", {125, 500}}, {"art_pulm", {125, 500}},
        {"resp", {62.5, 500}} },
      true },   // "TRUE" -> true
    { 4, "SHHS",
      { {"ecg", {125, 1000}}, {"eeg", {125, 500}}, {"eog_l", {50, 500}},
        {"eog_r", {50, 500}}, {"emg", {125, 500}}, {"flow", {10, 32}},
        {"thor", {10, 32}}, {"abdo", {10, 32}}, {"pos", {1, 1}},
        {"oxstatus", {1, 1}}, {"spo2", {1, 1}}, {"hr", {1, 1}} },
      true },   // "TRUE" -> true
};

// ---- every rate pair config_entry has ---------------------------------------
//
// "Absent channels report zero" means checking ALL of them, not just the ones
// a dataset uses -- so the test needs the full list. A pointer-to-member
// (double config_entry::*) lets one loop read any of them by name.
struct RatePair { const char* name; double config_entry::* raw; double config_entry::* up; };

static const RatePair kAllPairs[] = {
    { "ecg",       &config_entry::ecg_raw_rate,       &config_entry::ecg_upsample_rate },
    { "ppg",       &config_entry::ppg_raw_rate,       &config_entry::ppg_upsample_rate },
    { "cvp",       &config_entry::cvp_raw_rate,       &config_entry::cvp_upsample_rate },
    { "pres",      &config_entry::pres_raw_rate,      &config_entry::pres_upsample_rate },
    { "abp",       &config_entry::abp_raw_rate,       &config_entry::abp_upsample_rate },
    { "art",       &config_entry::art_raw_rate,       &config_entry::art_upsample_rate },
    { "art_pulm",  &config_entry::art_pulm_raw_rate,  &config_entry::art_pulm_upsample_rate },
    { "accel",     &config_entry::accel_raw_rate,     &config_entry::accel_upsample_rate },
    { "temp",      &config_entry::temp_raw_rate,      &config_entry::temp_upsample_rate },
    { "marker",    &config_entry::marker_raw_rate,    &config_entry::marker_upsample_rate },
    { "resp",      &config_entry::resp_raw_rate,      &config_entry::resp_upsample_rate },
    { "pacemaker", &config_entry::pacemaker_raw_rate, &config_entry::pacemaker_upsample_rate },
    { "eeg",       &config_entry::eeg_raw_rate,       &config_entry::eeg_upsample_rate },
    { "eeg_4",     &config_entry::eeg_4_raw_rate,     &config_entry::eeg_4_upsample_rate },
    { "flow",      &config_entry::flow_raw_rate,      &config_entry::flow_upsample_rate },
    { "snore",     &config_entry::snore_raw_rate,     &config_entry::snore_upsample_rate },
    { "thor",      &config_entry::thor_raw_rate,      &config_entry::thor_upsample_rate },
    { "abdo",      &config_entry::abdo_raw_rate,      &config_entry::abdo_upsample_rate },
    { "leg",       &config_entry::leg_raw_rate,       &config_entry::leg_upsample_rate },
    { "auxac",     &config_entry::auxac_raw_rate,     &config_entry::auxac_upsample_rate },
    { "therm",     &config_entry::therm_raw_rate,     &config_entry::therm_upsample_rate },
    { "pos",       &config_entry::pos_raw_rate,       &config_entry::pos_upsample_rate },
    { "oxstatus",  &config_entry::oxstatus_raw_rate,  &config_entry::oxstatus_upsample_rate },
    { "spo2",      &config_entry::spo2_raw_rate,      &config_entry::spo2_upsample_rate },
    { "hr",        &config_entry::hr_raw_rate,        &config_entry::hr_upsample_rate },
    { "dhr",       &config_entry::dhr_raw_rate,       &config_entry::dhr_upsample_rate },
    { "eog_l",     &config_entry::eog_l_raw_rate,     &config_entry::eog_l_upsample_rate },
    { "eog_r",     &config_entry::eog_r_raw_rate,     &config_entry::eog_r_upsample_rate },
    { "emg",       &config_entry::emg_raw_rate,       &config_entry::emg_upsample_rate },
};

// ---- arrange + act, shared by every test below ------------------------------
//
// load_config reads "config.csv" from the CURRENT folder and creates output
// folders next to it. So: copy the fixture into a scratch folder, switch into
// it, load, and switch back -- the fixture in data\ is never touched.
// THE FIXTURE MUST NOT MAKE load_config OPEN A WINDOW. A blank
// original_file_path or output_folder makes load_config open a folder dialog,
// and a test program has no running Qt application -- so Qt stops the whole
// program ("Must construct a QApplication before a QWidget") and every test
// after it is lost. Check those two rows first, and fail with the reason.
static bool fixturePathsFilled(const fs::path& fixture) {
    std::ifstream f(fixture);
    if (!f) { ADD_FAILURE() << "cannot open the fixture: " << fixture.string(); return false; }
    std::string line;
    bool ok = true;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string key = line.substr(0, line.find(','));
        if (key != "original_file_path" && key != "output_folder") continue;
        // every cell after the first must be non-blank
        size_t start = line.find(',');
        while (start != std::string::npos) {
            const size_t end = line.find(',', start + 1);
            const std::string cell = line.substr(start + 1,
                end == std::string::npos ? std::string::npos : end - start - 1);
            if (cell.find_first_not_of(" \t") == std::string::npos) {
                ADD_FAILURE() << "fixture " << fixture.string() << ": '" << key
                    << "' has a blank cell -- load_config would open a folder dialog "
                    "and crash the test program. Fill every cell in that row (e.g. "
                    "'input' / 'output').";
                ok = false;
                break;
            }
            start = end;
        }
    }
    return ok;
}

static bool loadFromFixture(int choice, config_entry& out) {
    const fs::path fixture = fs::path(TESTS_DATA_DIR) / "config" / "config.csv";
    if (!fixturePathsFilled(fixture)) return false;
    const fs::path scratch = fs::temp_directory_path() / "qtvi_test_config";
    fs::create_directories(scratch);
    fs::copy_file(fixture, scratch / "config.csv", fs::copy_options::overwrite_existing);

    const fs::path original = fs::current_path();
    fs::current_path(scratch);
    const bool ok = load_config(choice, out);
    fs::current_path(original);
    return ok;
}

// ---- the tests: one requirement each ----------------------------------------
//
// TEST(Group, Name): Test Explorer lists them as
// test1a_config_dataloading > AllDatasetsLoad, ...
// EXPECT_*  records a failure and carries on (the CHECK of before);
// ASSERT_*  records a failure and stops this test (the REQUIRE of before).
// SCOPED_TRACE adds "MESA" etc. to any failure inside the loop, so you can
// tell WHICH dataset failed.

TEST(test1a_config_dataloading, AllDatasetsLoad) {
    for (const Expected& e : kExpected) {
        SCOPED_TRACE(e.dataset);
        config_entry cfg;
        EXPECT_TRUE(loadFromFixture(e.choice, cfg));
        EXPECT_EQ(cfg.dataset_type, std::string(e.dataset));
    }
}

TEST(test1a_config_dataloading, EachChannelsRatePairIsRead) {
    for (const Expected& e : kExpected) {
        SCOPED_TRACE(e.dataset);
        config_entry cfg;
        ASSERT_TRUE(loadFromFixture(e.choice, cfg));   // no point checking rates if it did not load
        for (const auto& [name, rate] : e.rates) {
            SCOPED_TRACE(name);
            bool known = false;
            for (const RatePair& p : kAllPairs) {
                if (name != p.name) continue;
                known = true;
                EXPECT_EQ(cfg.*(p.raw), rate.first);
                EXPECT_EQ(cfg.*(p.up), rate.second);
            }
            // A name kAllPairs does not have (a typo in kExpected) would
            // otherwise be skipped, and the test would pass checking nothing.
            EXPECT_TRUE(known);
        }
    }
}

TEST(test1a_config_dataloading, AbsentChannelsReportZero) {
    for (const Expected& e : kExpected) {
        SCOPED_TRACE(e.dataset);
        config_entry cfg;
        ASSERT_TRUE(loadFromFixture(e.choice, cfg));
        for (const RatePair& p : kAllPairs) {
            if (e.rates.count(p.name)) continue;   // this dataset has it -- checked above
            SCOPED_TRACE(p.name);
            EXPECT_EQ(cfg.*(p.raw), 0.0);
            EXPECT_EQ(cfg.*(p.up), 0.0);
        }
    }
}

TEST(test1a_config_dataloading, UseConsensusRpeakDefaultsToTrueWhenBlank) {
    for (const Expected& e : kExpected) {
        SCOPED_TRACE(e.dataset);
        config_entry cfg;
        ASSERT_TRUE(loadFromFixture(e.choice, cfg));
        EXPECT_EQ(cfg.use_consensus_rpeak, e.consensus);
    }
}