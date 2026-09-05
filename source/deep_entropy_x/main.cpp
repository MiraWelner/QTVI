/**
 * @file   main.cpp
 * @brief  Entry point for the noise marking and template marking pipeline. Handles user input for dataset selection, initials, and file processing.
 */

#include "post_process.hpp"
#include "config_file_handling/config_loader.hpp"
#include "noise_marking_gui/gui_handler.h"
#include "noise_marking_gui/user_annotation_handler.h"
#include "template_marking_gui/parse_data_from_filename.hpp"
#include "logging/user_mark_log.hpp"
#include "template_marking_gui/TemplateViewerWindow.hpp"

#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <QScreen>
#include <QEventLoop>
#include <QObject>
#include <QFileInfo>

#include <algorithm>
#include <atomic>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <thread>
#include <functional>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QProgressDialog>
#include <vector>

#include <theme/theme.h>


 // Ask the user which dataset to load. The number returned here IS the dataType
 // that load_config() maps to a data_type row in config.csv, so this menu and the
 // mapping in config_loader.cpp have to be edited together. SHHS1/SHHS2 already
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
                markings->marking_type[i].toStdString(), sr);
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
static void runTemplateMarking(const config_entry& cfg, std::shared_ptr<post_process_detail::ViewerJob> job, const QString& fileId, std::function<void()> ensureWorkerDone) {
    const QString displayId = fileId;

    // ---- EVERY ALIGNMENT BUILT ONCE, BEFORE THE WINDOW OPENS ----------
    //
    // This was a LOOP over nAnchorPasses = anchorSequence().size() + 1: open a
    // viewer on the R templates, wait for "Finish and Next", tear the window
    // down, re-align on the next landmark, open a fresh viewer, repeat. Four
    // windows to mark one subject, with a rebuild between each.
    //
    // The viewer now loads all four alignments from one file and switches
    // between them per landmark (see anchor_view.hpp), so the alignments are
    // built up front and the window opens once. job->anchorStep is unused.
    //
    // Still off the GUI thread behind a modal busy dialog, for the same reason
    // as before: the align work is seconds, not milliseconds, and a frozen
    // window looks like a crash. And the finalize worker is still joined first
    // -- buildAllAnchors folds every block into job->tmpl and writes the
    // canonical file, which needs finalize's squared/absval done and job->tmpl
    // race-free. The old code could defer that join for the non-final steps
    // because they staged into job->anchorAccum; with one call there is no
    // staging, so the join is unconditional.
    {
        QProgressDialog busy(QStringLiteral("Building anchor alignments\xE2\x80\xA6"),
            QString(), 0, 0, nullptr);
        busy.setWindowModality(Qt::ApplicationModal);
        busy.setCancelButton(nullptr);
        busy.setMinimumDuration(0);
        busy.show();

        bool ok = false;
        QEventLoop waitLoop;
        std::thread build([&ok, job, &waitLoop, ensureWorkerDone]() {
            // The wait loop below blocks until this thread posts quit. Guard
            // the whole body so quit is ALWAYS posted -- if ensureWorkerDone()
            // or the build throws, we must still wake waitLoop or the GUI
            // thread hangs forever in waitLoop.exec().
            try {
                if (ensureWorkerDone) ensureWorkerDone();
                ok = post_process_detail::buildAllAnchors(*job);
            }
            catch (const std::exception& e) {
                ok = false;
                if (job->error.empty()) job->error = e.what();
            }
            catch (...) {
                ok = false;
                if (job->error.empty()) job->error = "unknown exception building anchor alignments";
            }
            QMetaObject::invokeMethod(&waitLoop, &QEventLoop::quit,
                Qt::QueuedConnection);   // guaranteed reached: no path skips it
            });
        waitLoop.exec();     // GUI stays alive; dialog animates
        build.join();
        busy.close();

        // NOT FATAL. job->viewerTemplatePath still points at a readable
        // templates file carrying the R base, and marking on R alone is what
        // this tool did before the anchor passes existed. The operator is told
        // because the consequence is invisible otherwise: TemplateBin::chFor
        // falls back to the R base for a missing alignment, so every bar's
        // close-up would silently show the R template and the per-alignment
        // CSV columns would be four copies of it.
        if (!ok) {
            QMessageBox::warning(nullptr, "Anchor alignments",
                QString("Could not build the anchor alignments:\n%1\n\n"
                    "Marking will continue on the R alignment only -- the "
                    "close-up panel will show the R template under every "
                    "landmark.")
                .arg(QString::fromStdString(job->error)));
        }
    }

    {
        TemplateViewerWindow viewer;
        viewer.setBoundaryTrainingDir(QString::fromStdString(cfg.training_log));
        // Was never wired up: m_vcgOutputPath defaulted to an empty QString,
        // so QDir(m_vcgOutputPath) in save_bin_and_csv() resolved to the
        // process's current working directory instead of cfg.vcg_output --
        // <id>_vcg.csv (and vcg_basis, if it shares this directory) landed
        // wherever the executable was launched from.
        viewer.setVcgOutputDir(QString::fromStdString(cfg.vcg_output));
        // Sections 5.2-5.4 normalization CSVs (<id>_feature_norm.csv,
        // <id>_cv_check.csv) land alongside the bin archive. Same
        // "setter must actually be called or the path stays empty and
        // nothing writes" lesson as setVcgOutputDir above.
        viewer.setNormOutputDir(QString::fromStdString(cfg.bin_archive_path));
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
        viewer.loadSubject(QString::fromStdString(job->viewerTemplatePath.string()),
            QString::fromStdString(cfg.qtvi_marker_path),
            displayId, cfg.ecg_upsample_rate,
            cfg.ppg_upsample_rate, cfg.abp_upsample_rate,
            cfg.art_upsample_rate, cfg.art_pulm_upsample_rate,
            cfg.notch_filter_hz);
        loop.exec();
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
    std::filesystem::create_directories(cfg.qtvi_marker_path);
    std::filesystem::create_directories(cfg.quality_metric);
    std::filesystem::create_directories(cfg.training_log);
    std::filesystem::create_directories(cfg.bin_archive_path);


    const std::vector<std::filesystem::path> binFiles = load_binfiles(cfg);
    if (binFiles.empty()) {
        std::cerr << "No .bin files in: " << cfg.bin_file_path << "\n";
        return 0;
    }

    // Background squared/absval finalize jobs that are still running. We do
    // NOT join a file's worker before moving to the next file -- that's what
    // caused the stall after the template file was saved. Instead the worker
    // is parked here, the loop advances immediately to the next file, and we
    // reap finished workers opportunistically (and drain the rest at exit).
    struct Outstanding {
        std::thread th;
        std::shared_ptr<post_process_detail::ViewerJob> job;
        std::shared_ptr<std::atomic<bool>> done;
    };
    std::vector<Outstanding> outstanding;

    auto finishJob = [](const std::shared_ptr<post_process_detail::ViewerJob>& job) {
        if (!job->error.empty()) {
            std::cerr << "  ERROR (squared/absval finalize) " << job->stem << ": "
                << job->error << "\n";
        }
        // Best-effort removal of that file's provisional viewer file.
        if (job->needsFinalize && job->provisionalPath != job->templatePath) {
            std::error_code ec;
            std::filesystem::remove(job->provisionalPath, ec);
        }
        };

    // Reap completed background jobs. force==true joins everything (used at
    // shutdown); force==false only joins workers that have already finished,
    // so it never blocks the main thread / the next file from loading.
    auto reap = [&](bool force) {
        for (size_t i = 0; i < outstanding.size();) {
            Outstanding& o = outstanding[i];
            if (force || o.done->load(std::memory_order_acquire)) {
                if (o.th.joinable()) o.th.join();
                finishJob(o.job);
                outstanding.erase(outstanding.begin() + i);
            }
            else {
                ++i;
            }
        }
        };

    // Safety valve: if marking races far ahead of compute, cap how many
    // background jobs (and their in-memory peak/template data) pile up by
    // blocking on the oldest. In normal use human marking is slower than
    // the compute, so this is rarely hit.
    constexpr size_t kMaxOutstanding = 4;

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
        // prepareViewerJob does anneal + raw R-peaks + raw/unfiltered/PPG
        // templates and writes a provisional file for the viewer. The
        // squared/absval R-peak detection and templating are deferred to a
        // worker thread that runs while the user marks templates.
        auto jobOpt = post_process_detail::prepareViewerJob(cfg, effBin,
            ecg1Inverted, ecg2Inverted, ecg3Inverted);
        // shared_ptr so the worker lambda safely co-owns the job; it is
        // joined later (reaped when finished, or drained at shutdown), so the
        // job data always outlives the worker.
        auto job = std::make_shared<post_process_detail::ViewerJob>(std::move(*jobOpt));

        auto done = std::make_shared<std::atomic<bool>>(false);
        std::thread worker;
        if (job->needsFinalize) {
            worker = std::thread([job, done] {
                // Pure compute + file writes to canonical paths. No Qt here.
                post_process_detail::finalizeViewerJob(*job);
                done->store(true, std::memory_order_release);
                });
        }

        // ---- Stage 3: template marking --------------------------------
        // buildAllAnchors (called at the top of runTemplateMarking, before the
        // window opens) needs the finalize worker joined: it folds every
        // alignment into job->tmpl and writes the canonical templates file, so
        // finalize's squared/absval must already be done and job->tmpl must be
        // race-free.
        //
        // THIS NOW HAPPENS ON EVERY FILE, not just the ones where the operator
        // reached the final anchor. The old cycle staged the non-final anchors
        // into job->anchorAccum precisely so the join could be deferred and the
        // R pass could open while the squared build was still running; with one
        // build up front there is nothing to stage, so the marking window waits
        // for finalize. The parking below therefore almost always takes the
        // second branch.
        auto ensureWorkerDone = [&worker]() { if (worker.joinable()) worker.join(); };
        runTemplateMarking(cfg, job, QString::fromStdString(job->stem), ensureWorkerDone);

        // The worker is normally already joined by the line above, so this
        // usually finishes the job here rather than parking it. Both branches
        // are kept: needsFinalize == false means no worker ever started, and a
        // build that threw before ensureWorkerDone() could leave one running.
        if (worker.joinable())
            outstanding.push_back(Outstanding{ std::move(worker), job, done });
        else if (job->needsFinalize)
            finishJob(job);   // worker was joined by buildAllAnchors

        // Clean up any background jobs that have already finished (non-blocking).
        reap(/*force=*/false);

        // Bound how many in-flight jobs may accumulate.
        while (outstanding.size() > kMaxOutstanding) {
            Outstanding& o = outstanding.front();
            if (o.th.joinable()) o.th.join();
            finishJob(o.job);
            outstanding.erase(outstanding.begin());
        }
    }

    // Drain remaining background jobs before exiting.
    if (!outstanding.empty()) {
        reap(/*force=*/true);
    }

    std::cout << "\nAll files processed.\n";
    return 0;
}