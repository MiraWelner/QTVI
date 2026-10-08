// test2h_record_summary.cpp -- Task H (Section 4.8): the per-record diagnostic summary.
//
// Loads the synthetic record SYNTH_001 from TESTS_DATA_DIR/diagnostic/ with the
// pipeline's own readers -- the same files a real run produces:
//   SYNTH_001_beats.bin                     templates_io::readBeatsBin
//   SYNTH_001_templates.bin                 templates_io::readTemplatesBin
//   SYNTH_001_peak_locations_all_beats.bin  read_output_binfile(peaks, annealed)
//   SYNTH_001_annealed.bin
// (written by make_task_h_test_data.cpp), runs Task A's ECG SQI and the
// summary, and checks it.
//
// AllPanelsPopulated is the spec's acceptance test: "Each processed record
// produces the summary with all listed panels populated from the Task A to D
// outputs." The others check each panel holds the right numbers. Expected
// values come from the loaded files themselves (R peaks, PVC flags, templates),
// so the test needs nothing beyond the four .bin files.

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "logging/record_summary.hpp"
#include "logging/record_summary_plot.hpp"
#include "template_generation/template_io.hpp"
#include "peak_finding/peakfinding_io.hpp"

#ifndef TESTS_DATA_DIR
#error "TESTS_DATA_DIR is not defined. QTVI_Tests -> Properties (All Configurations, All Platforms) -> C/C++ -> Preprocessor -> Preprocessor Definitions: add  TESTS_DATA_DIR=R\"($(ProjectDir)data)\""
#endif

namespace {

    const std::string kId = "SYNTH_001";

    std::string dataPath(const std::string& suffix) {
        return std::string(TESTS_DATA_DIR) + "/diagnostic/" + kId + suffix;
    }

    // ECG sample rate from the annealed header (after magic, version, bin count, PPG rate).
    double annealedEcgRate(const std::string& path) {
        std::ifstream f(path, std::ios::binary);
        double sr = 0.0;
        f.seekg(8 + 4 + 8 + 8);
        f.read(reinterpret_cast<char*>(&sr), 8);
        return sr;
    }

    // The record, as the pipeline would hand it to buildSummary.
    struct Record {
        double fs = 0;
        std::vector<output_binfile_data> peaks;
        template_structs::TemplateFile tf;
        template_structs::BeatsFile bf;
        std::vector<std::array<std::vector<size_t>, 3>> sliceOfRow;
        // Per beat number (bins in order, slices in order): bin, slice, PVC.
        std::vector<int> binOf;
        std::vector<size_t> sliceOf;
        std::vector<bool> pvc;
    };

    Record load() {
        Record rec;
        rec.peaks = read_output_binfile(dataPath("_peak_locations_all_beats.bin"), dataPath("_annealed.bin"));
        rec.fs = annealedEcgRate(dataPath("_annealed.bin"));
        const std::size_t nBins = rec.peaks.size();
        rec.tf.bins.resize(nBins);
        rec.sliceOfRow.resize(nBins);
        auto& rows = rec.bf.per_channel_beats["CH1"];
        auto& rhythm = rec.bf.per_channel_rhythm["CH1"];
        rows.resize(nBins);
        rhythm.resize(nBins);

        std::vector<templates_io::BinBlock<templates_io::BeatRecord>> beats;
        if (!templates_io::readBeatsBin(dataPath("_beats.bin"), beats))
            throw std::runtime_error("readBeatsBin failed: " + dataPath("_beats.bin"));
        for (const auto& blk : beats) {
            if (blk.channel != "CH1") continue;
            std::vector<size_t> slice(nBins, 0);
            for (std::size_t k = 0; k < blk.records.size(); ++k) {
                const auto& r = blk.records[k];
                const size_t s = slice[r.bin]++;
                rec.tf.bins[r.bin].ch1_raw.r_col = r.r_col;
                rec.binOf.push_back(static_cast<int>(r.bin));
                rec.sliceOf.push_back(s);
                rec.pvc.push_back(r.premature != 0);
                if (!r.became_beat) continue;
                rows[r.bin].emplace_back(blk.samples.begin() + k * blk.width,
                    blk.samples.begin() + (k + 1) * blk.width);
                rhythm[r.bin].push_back(r.premature ? 1 : 0);
                rec.sliceOfRow[r.bin][0].push_back(s);
            }
        }

        std::vector<templates_io::BinBlock<templates_io::TemplateRecord>> tmpls;
        if (!templates_io::readTemplatesBin(dataPath("_templates.bin"), tmpls))
            throw std::runtime_error("readTemplatesBin failed: " + dataPath("_templates.bin"));
        for (const auto& blk : tmpls) {
            if (blk.channel != "CH1") continue;
            for (std::size_t k = 0; k < blk.records.size(); ++k) {
                const auto& r = blk.records[k];
                auto& bank = rec.tf.bins[r.bin].ecg_bank[0];
                if (static_cast<int>(bank.templates.size()) <= r.template_id)
                    bank.templates.resize(r.template_id + 1);
                auto& t = bank.templates[r.template_id];
                t.label_code = r.label_code;
                t.spawn_seq = blk.trailers[k].spawn_seq;
                t.subtype = blk.trailers[k].subtype;
                t.confirmed_by_operator = blk.trailers[k].confirmed_by_operator != 0;
                t.members = blk.members[k];
                t.members_clean = blk.members_clean[k];
            }
        }
        return rec;
    }

    class test2h_record_summary : public ::testing::Test {
    protected:
        Record rec = load();
        EcgSQIResult sqi;

        record_summary::DiagnosticSummary build() {
            const LeadPolarity pol{ { false, false, false } };
            sqi = scoreEcgSQI(rec.tf, rec.bf, rec.fs, pol);     // Task A, as commit() runs it
            record_summary::SummaryInputs in;
            in.recordID = kId;
            in.tmpl = &rec.tf;
            in.beats = &rec.bf;
            in.peakResults = &rec.peaks;
            in.ecgSliceOfRow = &rec.sliceOfRow;
            in.sqi = &sqi;
            in.ecgFs = rec.fs;
            in.ecgRecRateHz = rec.fs;
            return record_summary::buildSummary(in);
        }

        std::size_t nBeats() const { return rec.pvc.size(); }
        // RR from beat n to the next, from the R peaks in the file.
        double rrMs(std::size_t n) const {
            const auto& R = rec.peaks[rec.binOf[n]].ch1.raw;
            const size_t s = rec.sliceOf[n];
            return 1000.0 * (static_cast<double>(R[s + 1]) - static_cast<double>(R[s])) / rec.fs;
        }
        bool sameBin(std::size_t n) const { return n + 1 < nBeats() && rec.binOf[n] == rec.binOf[n + 1]; }
    };

    std::size_t countFinite(const std::vector<double>& v) {
        return static_cast<std::size_t>(std::count_if(v.begin(), v.end(),
            [](double x) { return std::isfinite(x); }));
    }

    std::string slurp(const std::filesystem::path& p) {
        std::ifstream f(p);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    }

}  // namespace

// ---- ACCEPTANCE -----------------------------------------------------------

TEST_F(test2h_record_summary, AllPanelsPopulated) {
    const auto ds = build();
    EXPECT_EQ(ds.recordID, kId);
    EXPECT_EQ(ds.n_beats, nBeats());
    EXPECT_GT(countFinite(ds.chiSq0_vs_beat), 0u) << "chi-sq_0 vs beat empty";
    EXPECT_GT(countFinite(ds.chiSqAbs_vs_beat), 0u) << "chi-sq_abs vs beat empty";
    EXPECT_FALSE(ds.rr_vs_pp.empty()) << "RR vs PP empty";
    EXPECT_FALSE(ds.rr_vs_qt.empty()) << "RR vs QT empty";
    EXPECT_FALSE(ds.poincare.empty()) << "Poincare empty";
    EXPECT_FALSE(ds.normalBeats.empty()) << "normal beat tracings empty";
    EXPECT_FALSE(ds.abnormalBeats.empty()) << "abnormal beat tracings empty";
}

// ---- chi-sq: Task A's own score for every beat ----------------------------

TEST_F(test2h_record_summary, ChiSqMatchesTaskAScores) {
    const auto ds = build();
    EXPECT_EQ(countFinite(ds.chiSq0_vs_beat), nBeats());   // every beat is in a template
    for (const EcgSQIBeat& r : sqi.beats) {
        const size_t slice = rec.sliceOfRow[r.bin][0][r.row];
        std::size_t n = 0;
        while (n < nBeats() && !(rec.binOf[n] == static_cast<int>(r.bin) && rec.sliceOf[n] == slice)) ++n;
        ASSERT_LT(n, nBeats());
        EXPECT_DOUBLE_EQ(ds.chiSq0_vs_beat[n], r.q.chiSq0);
        EXPECT_DOUBLE_EQ(ds.chiSqAbs_vs_beat[n], r.q.chiSqAbs);
    }
}

// ---- RR vs PP: P wave to P wave ---------------------------------------------

TEST_F(test2h_record_summary, RrVsPpPairsConsecutiveSinusBeats) {
    const auto ds = build();
    std::vector<double> rr;
    for (std::size_t n = 0; n < nBeats(); ++n)
        if (sameBin(n) && !rec.pvc[n] && !rec.pvc[n + 1]) rr.push_back(rrMs(n));
    ASSERT_EQ(ds.rr_vs_pp.size(), rr.size());
    for (std::size_t i = 0; i < rr.size(); ++i) {
        EXPECT_NEAR(ds.rr_vs_pp[i].first, rr[i], 1e-6);
        // PP = RR + PR(n+1) - PR(n); PR varies by < 30 ms beat to beat here.
        EXPECT_NEAR(ds.rr_vs_pp[i].second, rr[i], 30.0) << "pair " << i;
    }
    // Every sinus beat's P sits before its R; no PVC has one.
    for (std::size_t n = 0; n < nBeats(); ++n) {
        if (rec.pvc[n]) { EXPECT_FALSE(std::isfinite(ds.p_re_r_ms[n])) << "PVC beat " << n; continue; }
        ASSERT_TRUE(std::isfinite(ds.p_re_r_ms[n])) << "sinus beat " << n;
        EXPECT_LT(ds.p_re_r_ms[n], -30.0);
        EXPECT_GT(ds.p_re_r_ms[n], -300.0);
    }
}

// ---- Poincare -------------------------------------------------------------------

TEST_F(test2h_record_summary, PoincareIsConsecutiveRrInsideBins) {
    const auto ds = build();
    std::vector<std::pair<double, double>> expect;
    for (std::size_t n = 0; n < nBeats(); ++n)
        if (sameBin(n)) expect.push_back({ rrMs(n), rrMs(n + 1) });
    ASSERT_EQ(ds.poincare.size(), expect.size());
    for (std::size_t i = 0; i < expect.size(); ++i) {
        EXPECT_NEAR(ds.poincare[i].first, expect[i].first, 1e-6);
        EXPECT_NEAR(ds.poincare[i].second, expect[i].second, 1e-6);
    }
}

TEST_F(test2h_record_summary, SpliceGapBreaksTheIntervalAcrossIt) {
    const auto before = build();
    // Cut 1000 samples out of bin 0 between its 11th and 12th R peaks
    // (beats 10 and 11 are both sinus, so one PP point goes).
    ASSERT_FALSE(rec.pvc[10] || rec.pvc[11]);
    const auto& R = rec.peaks[0].ch1.raw;
    const uint64_t first = rec.peaks[0].ecg_bin_indexs[0].first;
    const uint64_t last = rec.peaks[0].ecg_bin_indexs[0].second;
    const uint64_t cut = first + R[10] + 20;
    rec.peaks[0].ecg_bin_indexs = { { first, cut }, { cut + 1 + 1000, last + 1000 } };
    const auto after = build();
    EXPECT_EQ(after.rr_vs_pp.size(), before.rr_vs_pp.size() - 1);
    EXPECT_EQ(after.poincare.size(), before.poincare.size() - 2);
}

// ---- RR vs QT: one point per template -----------------------------------------------

TEST_F(test2h_record_summary, RrVsQtIsPerTemplate) {
    const auto ds = build();
    std::size_t nTemplates = 0;
    for (const auto& b : rec.tf.bins) nTemplates += b.ecg_bank[0].templates.size();
    ASSERT_EQ(ds.qt_templates.size(), nTemplates);
    ASSERT_EQ(ds.rr_vs_qt.size(), nTemplates);
    std::size_t ectopic = 0;
    for (const auto& t : ds.qt_templates) {
        EXPECT_TRUE(std::isfinite(t.rr_ms));
        EXPECT_GT(t.qt_ms, 0.0);
        EXPECT_NE(t.qt_source, record_summary::QtSource::NONE);
        ectopic += t.abnormal ? 1 : 0;
    }
    EXPECT_EQ(ectopic, rec.tf.bins.size());   // one PVC template per bin
}

// ---- tracings -------------------------------------------------------------------------

TEST_F(test2h_record_summary, TracingsSplitNormalFromAbnormal) {
    const auto ds = build();
    const std::size_t nPvc = static_cast<std::size_t>(std::count(rec.pvc.begin(), rec.pvc.end(), true));
    EXPECT_EQ(ds.normalBeats.size(), std::min(nBeats() - nPvc, record_summary::kMaxTracings));
    EXPECT_EQ(ds.abnormalBeats.size(), std::min(nPvc, record_summary::kMaxTracings));
    for (std::size_t n : ds.abnormalBeatNumbers) EXPECT_TRUE(rec.pvc[n]);
    for (std::size_t n : ds.normalBeatNumbers) EXPECT_FALSE(rec.pvc[n]);
}

// ---- the output files ---------------------------------------------------------------------

TEST_F(test2h_record_summary, WritesCsvQtCsvAndSvg) {
    const auto ds = build();
    const std::filesystem::path dir = std::filesystem::temp_directory_path() / "qtvi_test2h";
    std::filesystem::create_directories(dir);
    ASSERT_TRUE(record_summary::writeSummaryCsv(dir.string(), kId, ds));
    ASSERT_TRUE(record_summary_plot::writeSummarySvg(dir.string(), kId, ds));
    const std::string c = slurp(dir / (kId + "_diagnostic_summary.csv"));
    for (const char* panel : { "\nchisq0,", "\nchisqabs,", "\nrr_vs_pp,", "\npoincare,",
                               "\nrr_vs_qt,", "\nnormal_beats,", "\nabnormal_beats," })
        EXPECT_NE(c.find(panel), std::string::npos) << "CSV has no rows for " << panel + 1;
    EXPECT_NE(slurp(dir / (kId + "_diagnostic_summary.svg")).find("</svg>"), std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(dir / (kId + "_diagnostic_summary_qt.csv")));
    std::filesystem::remove_all(dir);
}