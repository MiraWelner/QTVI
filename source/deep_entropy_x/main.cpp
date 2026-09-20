/**
 * @file   main.cpp
 * @brief  Entry point for the noise marking and template marking pipeline. Handles user input for dataset selection, initials, and file processing.
 */

#include "analysis_job.hpp"
#include "config_file_handling/config_loader.hpp"
#include "noise_marking_gui/gui_handler.h"
#include "noise_marking_gui/user_annotation_handler.h"
#include "template_marking_gui/parse_data_from_filename.hpp"
#include "logging/user_mark_log.hpp"
#include "template_viewer/template_viewer.hpp"

#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <QScreen>
#include <QEventLoop>
#include <QObject>
#include <QFileInfo>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <thread>
#include <QtWidgets/QMessageBox>
#include <vector>

#include <theme/theme.h>


 // Ask the user which dataset to load. The number returned here IS the dataType
 // that load_config() maps to a data_type row in config.csv, so this menu and the
 // mapping in config_loader.cpp have to be edited together. SHHS already
 // had channel labels in apply_dataset_specific_channel_labels() and rows in the
 // config, but this menu was still pinned to three datasets, so there was no way
 // to select them.
static int get_dataset_choice() {
    static constexpr int n_valid_datasets = 4;
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\n"
        "4: SHHS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return -1;
    while (choice < 1 || choice > n_valid_datasets) {
        std::cout << "Please enter 1-" << n_valid_datasets << ": ";
        if (!(std::cin >> choice)) return -1;
    }
    return choice;
}

static std::string get_initials() {
    std::cout << "Enter your initials so it can be identified who logged what: ";
    std::string raw;
    std::getline(std::cin >> std::ws, raw);
    std::string out;
    for (char c : raw)
        if (std::isalnum(static_cast<unsigned char>(c)))
            out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (out.empty()) out = "anon";
    return out;
}

static std::vector<std::filesystem::path> load_binfiles(const config_entry& cfg) {
    std::vector<std::filesystem::path> binFiles;
    for (const auto& entry :
        std::filesystem::directory_iterator(cfg.bin_file_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin")
            binFiles.push_back(entry.path());
    }
    std::sort(binFiles.begin(), binFiles.end());
    return binFiles;
}

static bool runNoiseMarking(const config_entry& cfg, const std::filesystem::path& binFs, QVector<GenExcStruct>& outAll,
    std::filesystem::path& outCurrent, beat_log& beatLog, bool& outEcg1Inverted, bool& outEcg2Inverted, bool& outEcg3Inverted) {
    /*Launch the GUI to do the noise marking. One important thing that takes place is that the gui object (a noise_marking_gui) has
    an invertedForSignal attribute for each channel. */
    auto gui = std::make_unique<noise_marking_gui>();
    gui->set_params_to_config_defaults(cfg);
    gui->setBeatLog(&beatLog);

    QScreen* screen = QGuiApplication::primaryScreen();
    if (screen) {
        const QRect available = screen->availableGeometry();
        gui->setMaximumHeight(available.height() - 40);
    }
    gui->setWindowTitle(
        QString::fromStdString("Marking: " + binFs.filename().string()));
    gui->setFileSource(QString::fromStdString(binFs.string()));

    if (gui->exec() != QDialog::Accepted)
        return false;

    outAll = gui->getAllMarkings();
    outCurrent = gui->getFilePath().toStdString();

    outEcg1Inverted = gui->invertedForSignal("ECG1");
    outEcg2Inverted = gui->invertedForSignal("ECG2");
    outEcg3Inverted = gui->invertedForSignal("ECG3");

    return true;
}

static void exportMarkings(const config_entry& cfg, const std::filesystem::path& binFile, const GenExcStruct* markings) {
    /*Export the recorded markings to <noise_data_path>/<stem>_noise_markings.{csv,bin}.
    The .bin is the file the anneal step (processOneFile) reads back in.*/
    annotation_handler nm;
    auto rateForLabel = [&](const QString& label) -> double {
        if (label == "PPG")      return cfg.ppg_upsample_rate;
        if (label == "ABP")      return cfg.abp_upsample_rate;
        if (label == "ACCEL")    return cfg.accel_upsample_rate;
        if (label == "ART")      return cfg.art_upsample_rate;
        if (label == "ART_PULM") return cfg.art_pulm_upsample_rate;
        return cfg.ecg_upsample_rate;   // ECG1/ECG2/ECG3
        };

    if (markings) {
        for (int i = 0; i < markings->noiseExc.size(); ++i) {
            const double sr = rateForLabel(markings->data_type[i]);
            nm.addSegment(
                static_cast<int>(markings->noiseExc[i].first * sr),
                static_cast<int>(markings->noiseExc[i].second * sr),
                markings->data_type[i].toStdString(),
                markings->marking_type[i].toStdString(), sr,
                markings->threshold[i], markings->blanking[i]);
        }
    }
    const std::filesystem::path base =
        std::filesystem::path(cfg.noise_data_path)
        / (binFile.stem().string() + "_noise_markings");
    nm.exportCSV(base.string() + ".csv");
    nm.exportBinary(base.string() + ".bin");
    std::cout << "Saved Noise Markings for " << binFile.filename().string() << "\n";
}

// ---------------------------------------------------------------------------
// Stage 3: template marking. Builds every anchor alignment, then opens ONE
// TemplateViewerWindow on the templates file and blocks (local event loop)
// until the viewer finishes. It used to open one window per alignment.
// ---------------------------------------------------------------------------
static void runTemplateMarking(const config_entry& cfg, std::shared_ptr<analysis_job::AnalysisJob> job, const QString& fileId,
    std::vector<analysis_job::BankSnapshot>& outBanks) {
    const QString displayId = fileId;

    // ---- EVERY ALIGNMENT BUILT ONCE, BEFORE THE WINDOW OPENS ----------
    // The anchor alignments are built in analysis_job::prepare now -- all four are
    // folded into job.tmpl before the provisional templates file is written,
    // while job.beats is still the pristine R-pass matrix. So there is nothing
    // to do here: the file the viewer is about to open already carries them.
    //
    // The buildAllAnchors call that used to sit here (on its own thread, behind
    // a progress dialog) is gone with the function; post_process.hpp does that
    // work upstream instead.

    {
        TemplateViewerWindow viewer;
        viewer.setBoundaryTrainingDir(QString::fromStdString(cfg.training_log));
        // Was never wired up: m_vcgOutputPath defaulted to an empty QString,
        // so QDir(m_vcgOutputPath) in save_bin_and_csv() resolved to the
        // process's current working directory instead of cfg.vcg_output --
        // <id>_vcg.csv (and vcg_basis, if it shares this directory) landed
        // wherever the executable was launched from.
        viewer.set_vcg_output_dir(QString::fromStdString(cfg.vcg_output));
        // Sections 5.2-5.4 normalization CSVs (<id>_feature_norm.csv,
        // <id>_cv_check.csv) land alongside the bin archive. Same
        // "setter must actually be called or the path stays empty and
        // nothing writes" lesson as setVcgOutputDir above.
        viewer.setNormOutputDir(QString::fromStdString(cfg.template_path));
        // (setAnchorPassCount / setAnchorStep / setAnchorLabel are gone with the
        //  cycle -- there is no pass, no step, and no single alignment to name
        //  in the title bar. The focus panel names the alignment it is showing.)

        QEventLoop loop;

        // Queued connection so that even if loadSubject() emits finished()
        // synchronously (its read-error / empty-bins path), the quit is still
        // delivered once exec() starts.
        QObject::connect(&viewer, &TemplateViewerWindow::finished,
            &loop, &QEventLoop::quit, Qt::QueuedConnection);
        // (requestQAlignReload is gone: one save writes every alignment, so
        //  there is nothing to reload and the button always reads "Finish".)

        viewer.show();
        // FROM MEMORY, NOT FROM A FILE. This passed
        // job->viewerTemplatePath, which meant analysis_job::prepare had to write
        // templates.bin before the squared/absval blocks existed -- a file
        // that looked complete and was not. The TemplateFile is already in
        // this process; the disk round trip produced nothing but a filename.
        //
        // templates.bin is now written once, at the end of the anchor cycle,
        // so its existence means complete and viewerTemplatePath is no longer
        // read here.
        viewer.loadSubject(job->tmpl,
            QString::fromStdString(cfg.template_path),
            QString::fromStdString(cfg.fiducial_marker_locations),
            displayId, cfg.ecg_upsample_rate,
            cfg.ppg_upsample_rate, cfg.abp_upsample_rate,
            cfg.art_upsample_rate, cfg.art_pulm_upsample_rate,
            cfg.notch_filter_hz);
        loop.exec();

        // ---- THE OPERATOR'S BANKS, OUT TO THE CALLER -------------------
        //
        // The confirmation flags live on the VIEWER's m_bins -- showPage sets
        // confirmed_by_operator as each panel is built, and a right-click sets
        // marked_invalid_template -- while job->tmpl is a separate TemplateFile
        // that nothing touches after prepare(). They have to be copied back or
        // the `confirmed` column reads "presumed" for the whole record no
        // matter how much marking was done.
        //
        // The copy-back AND the write used to happen right here, and both were
        // wrong in the same way: the finalize worker is very often STILL
        // RUNNING at this point (nothing joined it -- ensureWorkerDone was
        // built and never called), and mergeTemplatesSlow mutates job->tmpl.
        // So the write raced the worker over the same object, and which blocks
        // reached disk depended on the timing of a background thread.
        //
        // The banks are handed out instead, and the caller joins the worker
        // before calling analysis_job::commit. Copied, not referenced: the
        // viewer is destroyed at the closing brace below.
        outBanks.clear();
        outBanks.reserve(viewer.bins().size());
        for (const TemplateBin& b : viewer.bins())
            outBanks.push_back(analysis_job::BankSnapshot{ b.ecg_bank, b.ppg_bank });
        // viewer is destroyed here (window closes).
    }
}


int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    Theme::apply(app);
    const int dataset_choice = get_dataset_choice();
    config_entry cfg;
    if (!load_config(dataset_choice, cfg)) {
        std::cerr << "Error Loading config.csv\n";
        return 1;
    }
    // Per-reviewer log folder. Each reviewer only skips files THEY logged.
    const std::string initials = get_initials();

    cfg.log_path = cfg.output_path + "/log_" + initials;
    std::filesystem::create_directories(cfg.log_path);
    std::cout << "Logging to: " << cfg.log_path << "\n";

    std::filesystem::create_directories(cfg.annealed_data_path);
    std::filesystem::create_directories(cfg.r_peak_data_path);
    std::filesystem::create_directories(cfg.noise_data_path);
    std::filesystem::create_directories(cfg.template_path);
    std::filesystem::create_directories(cfg.fiducial_marker_locations);
    std::filesystem::create_directories(cfg.quality_metric);
    std::filesystem::create_directories(cfg.training_log);
    std::filesystem::create_directories(cfg.template_path);


    const std::vector<std::filesystem::path> binFiles = load_binfiles(cfg);
    if (binFiles.empty()) {
        std::cerr << "No .bin files in: " << cfg.bin_file_path << "\n";
        return 0;
    }

    // NO WORKER PARKING. There used to be a vector of outstanding threads here,
    // with a reap()/force-drain pair and a kMaxOutstanding cap, so the loop could
    // advance to the next file while a finalize worker was still running.
    //
    // That is no longer possible, and the reason is worth stating: _bins.bin is
    // written by analysis_job::commit, commit has to follow the join (it
    // serializes the very object finalize mutates), and it has to precede the
    // next file. A worker therefore cannot outlive its own iteration, so there
    // is nothing left to park, reap or cap.
    //
    // What that costs: one file's squared/absval compute no longer overlaps the
    // NEXT file's marking. It still overlaps its OWN marking, which is where the
    // time actually was -- the old cap comment said as much ("in normal use human
    // marking is slower than the compute, so this is rarely hit").

    auto reportFinalizeError = [](const analysis_job::AnalysisJob& job) {
        if (!job.error.empty())
            std::cerr << "  ERROR (squared/absval finalize) " << job.stem << ": "
            << job.error << "\n";
        // NO CLEANUP. This removed the _templates.partial.bin that prepare()
        // used to write for the viewer to open. There is no provisional file:
        // _bins.bin is written once, by commit(), so nothing is left over.
        };

    for (const std::filesystem::path& binFs : binFiles) {
        const std::string stem = binFs.stem().string();
        std::cout << "\n=== " << stem << " ===\n";

        // Skip files that already have a log: an existing log means this file
        // was marked on a prior run, so move on to the next one.
        const std::string logPath = cfg.log_path + "/" + stem + "_log.csv";
        if (std::filesystem::exists(logPath)) {
            std::cout << "  log already exists in " << cfg.log_path << "; skipping " << stem << "\n";
            continue;
        }

        // ---- Per-file beat log (placeholder rows seeded with config defaults) ----
        // The GUI fills this live while the dialog is open and flushes/writes
        // it every 30 s; here we just seed the blanking/threshold columns.
        beat_log beatLog;
        beatLog.setDefaultParams(cfg.blanking_period, cfg.threshold);

        // ---- Stage 1: noise marking -------------------------------------
        std::cout << "Noise marking: " << binFs.filename().string() << "\n";
        QVector<GenExcStruct> allMarkings;
        std::filesystem::path currentBinFile;
        bool ecg1Inverted = false, ecg2Inverted = false, ecg3Inverted = false;
        if (!runNoiseMarking(cfg, binFs, allMarkings, currentBinFile, beatLog,
            ecg1Inverted, ecg2Inverted, ecg3Inverted)) {
            std::cout << "  skipped by user; not processing/templating.\n";
            continue;
        }

        // The user may have used "Load" inside the noise GUI to switch to a
        // different file; template/log the file actually marked, not the loop's
        // binFs. (In the normal case getFilePath() == binFs, so this is a no-op.)
        const std::filesystem::path effBin =
            currentBinFile.empty() ? binFs : currentBinFile;
        const std::string effStem = effBin.stem().string();

        // Commit the final partial buffer (the last <30 s the timer didn't
        // reach) and write the populated log out.
        beatLog.flushPending();
        beatLog.writeCsv(cfg.log_path + "/" + effStem + "_log.csv");

        // ---- Export markings (input to the anneal step) -----------------
        if (allMarkings.isEmpty()) {
            exportMarkings(cfg, currentBinFile, nullptr);
        }
        else {
            for (const GenExcStruct& m : allMarkings)
                exportMarkings(cfg, std::filesystem::path(m.filePath.toStdString()), &m);
        }

        // ---- Stage 2 (fast) + Stage 3 (marking), with the squared/absval
        //      half of Stage 2 running in parallel --------------------------
        // analysis_job::prepare does anneal, the hardware lag, raw R-peaks and
        // the raw/unfiltered/PPG templates, and hands back in-memory state --
        // no provisional file. The squared/absval R-peak detection and
        // templating are deferred to a worker that runs while the operator
        // marks templates.
        auto jobOpt = analysis_job::prepare(cfg, effBin,
            ecg1Inverted, ecg2Inverted, ecg3Inverted);
        if (!jobOpt) {
            // prepare() returns nullopt when the recording is shorter than one
            // bin; it has already said so on stderr. Dereferencing it was
            // unconditional here, so that case was a crash.
            std::cout << "  prep failed or skipped; nothing to mark.\n";
            continue;
        }
        // shared_ptr so the worker lambda co-owns the job: the job data always
        // outlives the worker, whatever order they unwind in.
        auto job = std::make_shared<analysis_job::AnalysisJob>(std::move(*jobOpt));

        // The `done` atomic IS GONE with the parking machinery. It existed
        // purely so reap() could test a worker for completion without blocking;
        // with a single worker that is always joined below, std::thread::join
        // is the whole of what is needed.
        std::thread worker;
        if (job->needsFinalize) {
            worker = std::thread([job] {
                // Pure compute + file writes to canonical paths. No Qt here.
                analysis_job::finalize(*job);
                });
        }

        // ---- Stage 3: template marking --------------------------------
        std::vector<analysis_job::BankSnapshot> operatorBanks;
        runTemplateMarking(cfg, job, QString::fromStdString(job->stem), operatorBanks);

        // ---- Commit: join, THEN write -----------------------------------
        //
        // THE JOIN IS A CORRECTNESS REQUIREMENT, NOT A PERFORMANCE CHOICE.
        // commit() writes job->tmpl; finalize() mutates job->tmpl, because
        // mergeTemplatesSlow packs the squared/absval blocks into it. Running
        // the two concurrently is a data race on the object being serialized,
        // and it fails by producing a plausible file rather than by crashing.
        //
        // That race was live. The old inline write sat immediately after the
        // viewer's event loop with nothing joining the worker first --
        // ensureWorkerDone was constructed at this line and never called -- so
        // which blocks reached disk depended on the timing of a background
        // thread.
        //
        // In practice the operator has been marking for minutes and the worker
        // finished long ago, so this join returns immediately.
        if (worker.joinable()) worker.join();
        if (job->needsFinalize) reportFinalizeError(*job);
        analysis_job::commit(*job, operatorBanks);
    }

    std::cout << "\nAll files processed.\n";
    return 0;
}