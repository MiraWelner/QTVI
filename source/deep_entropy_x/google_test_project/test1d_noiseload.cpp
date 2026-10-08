// ============================================================================
// test1a_noiseload.cpp -- QTVI_Tests (Google Test)
//
// REQUIREMENT: a subject's noise markings are read properly.
//
// WHAT "READ" MEANS. The noise-marking GUI writes two files on Finish:
//   <stem>_noise_markings.bin   read back by EVERYTHING -- the GUI's reload,
//                               the anneal step, the template build -- through
//                               noise_markings::loadRows
//   <stem>_noise_markings.csv   the same rows as text, for people; never read
// So the test reads the .bin with loadRows and checks it against its CSV: the
// CSV is the answer key, written by the same writer at the same moment.
//
// FIXTURE: QTVI_Tests\data\noise\ holds each subject as a PAIR:
//   3010023_20110817_noise_markings.bin
//   3010023_20110817_noise_markings.csv
// Every CSV there is tested; its .bin must sit next to it. To add a subject,
// copy both files from noise_marking_output into data\noise\.
//
// PROJECT SETTINGS: the same as test1a_config_dataloading (TESTS_DATA_DIR,
// include directories). Additional Include Directories must also reach
// noise_marking_gui\, and the project needs user_annotation_handler.cpp only
// if a test WRITES markings -- these only read, so it does not.
// ============================================================================
#include "pch.h"
#include "../noise_marking_gui/user_annotation_handler.h"   // noise_markings::loadRows
#include "../noise_marking_gui/annotation_types.hpp"        // label -> code

#ifndef TESTS_DATA_DIR
#error "TESTS_DATA_DIR is not defined. QTVI_Tests -> Properties (All Configurations, All Platforms) -> C/C++ -> Preprocessor -> Preprocessor Definitions: add  TESTS_DATA_DIR=R\"($(ProjectDir)data)\""
#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
namespace nm = noise_markings;

// ---- known answers, independent of both files ------------------------------
//
// What this subject's noise markings ARE, written out by hand. The .bin-vs-CSV
// test cannot catch a writer bug that puts the same wrong value in both files;
// this can. One entry per subject you have checked by eye.
struct KnownMark { const char* channel; const char* type; double startSec, endSec; double thr, blk; };
static const std::map<std::string, std::vector<KnownMark>> kKnown = {
    { "3010023_20110817", {
        { "ECG1", "1) R Peak Noise",     4.830, 16.412, NAN,   NAN   },
        { "ECG1", "3) Blank.+Thresh.",  19.872, 25.302, 0.700, 900.0 },
        { "PPG",  "2) Minor Noise",     13.168, 19.584, NAN,   NAN   },
        { "PPG",  "4) PVC",              6.080,  8.866, NAN,   NAN   },
        { "PPG",  "4) PVC",             23.500, 26.240, NAN,   NAN   },
    } },
};

// ---- the CSV, as the answer key ---------------------------------------------
struct CsvRow {
    double startSample, endSample, startSec, endSec;
    std::string channel, type;
    double thr, blk;   // NaN where the CSV cell is blank
};

static std::vector<CsvRow> readCsv(const fs::path& p) {
    std::vector<CsvRow> out;
    std::ifstream f(p);
    std::string line;
    std::getline(f, line);                                   // header
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> c;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')) c.push_back(cell);
        while (c.size() < 8) c.push_back("");
        auto num = [](const std::string& s) { return s.empty() ? std::nan("") : std::stod(s); };
        out.push_back({ num(c[0]), num(c[1]), num(c[2]), num(c[3]), c[4], c[5], num(c[6]), num(c[7]) });
    }
    return out;
}

// Every CSV in data\noise\, with the .bin that must sit next to it.
struct Pair { std::string stem; fs::path csv, bin; };
static std::vector<Pair> fixturePairs() {
    std::vector<Pair> out;
    const std::string tail = "_noise_markings.csv";
    for (const auto& e : fs::directory_iterator(fs::path(TESTS_DATA_DIR) / "noise")) {
        const std::string name = e.path().filename().string();
        if (name.size() <= tail.size() || name.compare(name.size() - tail.size(), tail.size(), tail) != 0)
            continue;
        const std::string stem = name.substr(0, name.size() - tail.size());
        out.push_back({ stem, e.path(), e.path().parent_path() / (stem + "_noise_markings.bin") });
    }
    return out;
}

// The CSV prints seconds, threshold and blanking to 3 decimals; the .bin has
// them exact. Half the last printed digit is the most they can differ.
static constexpr double kCsvRounding = 0.0005;

// NaN in both counts as equal: a blank threshold is NaN in the .bin.
static bool sameValue(double a, double b, double tol) {
    return (std::isnan(a) && std::isnan(b)) || std::abs(a - b) <= tol;
}

// ---- the tests ----------------------------------------------------------------

TEST(test1a_noiseload, EverySubjectHasItsBin) {
    const std::vector<Pair> pairs = fixturePairs();
    ASSERT_FALSE(pairs.empty()) << "no *_noise_markings.csv in " << (fs::path(TESTS_DATA_DIR) / "noise").string();
    for (const Pair& p : pairs) {
        SCOPED_TRACE(p.stem);
        EXPECT_TRUE(fs::exists(p.bin)) << "missing " << p.bin.string()
            << " -- copy it from noise_marking_output next to its CSV";
    }
}

TEST(test1a_noiseload, BinLoads) {
    for (const Pair& p : fixturePairs()) {
        SCOPED_TRACE(p.stem);
        if (!fs::exists(p.bin)) continue;                    // reported by EverySubjectHasItsBin
        const nm::RowsResult r = nm::loadRows(p.bin.string());
        EXPECT_TRUE(r.read) << r.error;
        EXPECT_TRUE(r.error.empty()) << r.error;
        EXPECT_FALSE(r.legacy) << "written in the old pre-versioned format";
    }
}

TEST(test1a_noiseload, BinMatchesItsCsvRowForRow) {
    for (const Pair& p : fixturePairs()) {
        SCOPED_TRACE(p.stem);
        if (!fs::exists(p.bin)) continue;
        const nm::RowsResult r = nm::loadRows(p.bin.string());
        ASSERT_TRUE(r.read) << r.error;
        const std::vector<CsvRow> csv = readCsv(p.csv);
        ASSERT_EQ(r.rows.size(), csv.size()) << "the .bin and the CSV hold different numbers of marks";
        // Same writer, same loop: the two files list the marks in the same order.
        for (size_t i = 0; i < csv.size(); ++i) {
            SCOPED_TRACE("row " + std::to_string(i + 1) + ": " + csv[i].channel + " " + csv[i].type);
            const nm::Row& b = r.rows[i];
            EXPECT_EQ(int(b.channel_code), int(nm::code_for_channel(csv[i].channel)));
            EXPECT_EQ(int(b.annotation_code), int(annotation_types::code_for_label(csv[i].type.c_str())));
            EXPECT_EQ(b.start_sample, csv[i].startSample);
            EXPECT_EQ(b.end_sample, csv[i].endSample);
            EXPECT_NEAR(b.start_sec, csv[i].startSec, kCsvRounding);
            EXPECT_NEAR(b.end_sec, csv[i].endSec, kCsvRounding);
            EXPECT_TRUE(sameValue(b.threshold, csv[i].thr, kCsvRounding))
                << "threshold: .bin " << b.threshold << ", CSV " << csv[i].thr;
            EXPECT_TRUE(sameValue(b.blanking_ms, csv[i].blk, kCsvRounding))
                << "blanking_ms: .bin " << b.blanking_ms << ", CSV " << csv[i].blk;
        }
    }
}

TEST(test1a_noiseload, KnownSubjectsHaveTheirKnownMarks) {
    for (const auto& [stem, marks] : kKnown) {
        SCOPED_TRACE(stem);
        const fs::path bin = fs::path(TESTS_DATA_DIR) / "noise" / (stem + "_noise_markings.bin");
        ASSERT_TRUE(fs::exists(bin)) << "kKnown lists " << stem << " but " << bin.string() << " is missing";
        const nm::RowsResult r = nm::loadRows(bin.string());
        ASSERT_TRUE(r.read) << r.error;
        ASSERT_EQ(r.rows.size(), marks.size());
        for (size_t i = 0; i < marks.size(); ++i) {
            SCOPED_TRACE("mark " + std::to_string(i + 1) + ": " + marks[i].channel + " " + marks[i].type);
            const nm::Row& b = r.rows[i];
            EXPECT_EQ(int(b.channel_code), int(nm::code_for_channel(marks[i].channel)));
            EXPECT_EQ(int(b.annotation_code), int(annotation_types::code_for_label(marks[i].type)));
            EXPECT_NEAR(b.start_sec, marks[i].startSec, kCsvRounding);
            EXPECT_NEAR(b.end_sec, marks[i].endSec, kCsvRounding);
            EXPECT_TRUE(sameValue(b.threshold, marks[i].thr, kCsvRounding)) << "threshold " << b.threshold;
            EXPECT_TRUE(sameValue(b.blanking_ms, marks[i].blk, kCsvRounding)) << "blanking_ms " << b.blanking_ms;
        }
    }
}