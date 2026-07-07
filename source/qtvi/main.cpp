/**
 * @file   main.cpp
 * @brief  Entry point for the noise marking and template marking pipeline. Handles user input for dataset selection, initials, and file processing.
 */

#include "post_process.hpp"
#include "config_loader.hpp"
#include "noise_marking_gui\gui_handler.h"
#include "noise_marking_gui\user_annotation_handler.h"
#include "template_marking_gui\parse_data_from_filename.hpp"
#include "beat_log.hpp"
#include "template_marking_gui\TemplateViewerWindow.hpp"

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
#include <vector>

#include <theme/theme.h>


static int get_dataset_choice() {
    /* Ask the user to select a MESA, Bittium, or CHAOS dataset. */
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return -1;
    while (choice < 1 || choice > 3) {
        std::cout << "Invalid choice. Please enter 1, 2, or 3: ";
        if (!(std::cin >> choice)) return -1;
    }
    return choice;
}

static std::string get_initials() {
    std::cout << "Enter your initials: ";
    std::string raw;
    std::getline(std::cin >> std::ws, raw);
    std::string out;                       // keep alphanumerics only, lowercased
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

// ---------------------------------------------------------------------------
// Stage 1: noise marking. Returns false if the user skipped the dialog
// (rejected), in which case the rest of the pipeline is skipped for this file.
// On accept, fills `outAll` with every file's markings touched in the dialog
// (the dialog lets the user browse to other files) and `outCurrent` with the
// path that was actually loaded when the dialog closed.
// ---------------------------------------------------------------------------
static bool runNoiseMarking(const config_entry& cfg,
    const std::filesystem::path& binFs,
    QVector<GenExcStruct>& outAll,
    std::filesystem::path& outCurrent,
    beat_log& beatLog) {
    auto gui = std::make_unique<noise_marking_gui>();
    gui->set_params_to_config_defaults(cfg);
    gui->setBeatLog(&beatLog);
    // NOTE: no setPostProcessQueue() -- there is no queue anymore. The GUI's
    // m_postQueue stays nullptr and every `if (m_postQueue)` branch is skipped.

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
    return true;   // gui destroyed here, before the heavy pipeline runs
}

// ---------------------------------------------------------------------------
// Export the recorded markings to <noise_data_path>/<stem>_noise_markings.{csv,bin}.
// The .bin is the file the anneal step (processOneFile) reads back in.
// ---------------------------------------------------------------------------
static void exportMarkings(const config_entry& cfg,
    const std::filesystem::path& binFile,
    const GenExcStruct* markings) {
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
    std::cout << "  saved markings for " << binFile.filename().string() << "\n";
}

// ---------------------------------------------------------------------------
// Stage 3: template marking. Opens TemplateViewerWindow on the templates file
// produced by Stage 2 and blocks (local event loop) until the viewer finishes.
// ---------------------------------------------------------------------------
static void runTemplateMarking(const config_entry& cfg,
    const std::filesystem::path& templateFile) {
    // Build the same display id template_marking.cpp used, so the saved
    // "<id>_template_markings.bin" name is unchanged. Fall back to the raw
    // stem when the filename doesn't parse (e.g. unexpected token layout).
    QString nameStem = QString::fromStdString(templateFile.stem().string());
    const int suffixPos = nameStem.indexOf("_templates");
    if (suffixPos >= 0) nameStem = nameStem.left(suffixPos);

    QString displayId = nameStem;
    if (auto parsed = parseTemplateFileName(nameStem)) {
        displayId = parsed->date.isEmpty()
            ? parsed->subjectId
            : parsed->subjectId + "_" + parsed->date;
    }

    std::cout << "Template marking: " << displayId.toStdString() << "\n";

    TemplateViewerWindow viewer;

    QEventLoop loop;
    // Queued connection so that even if loadSubject() emits finished()
    // synchronously (its read-error / empty-bins path), the quit is still
    // delivered once exec() starts -- avoids the hang the standalone
    // template_marking.cpp had in that edge case.
    QObject::connect(&viewer, &TemplateViewerWindow::finished,
        &loop, &QEventLoop::quit, Qt::QueuedConnection);

    viewer.show();
    viewer.loadSubject(QString::fromStdString(templateFile.string()),
        QString::fromStdString(cfg.qtvi_marker_path),
        displayId, cfg.ecg_upsample_rate);
    loop.exec();
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
        if (!runNoiseMarking(cfg, binFs, allMarkings, currentBinFile, beatLog)) {
            std::cout << "  skipped by user; not processing/templating.\n";
            continue;
        }

        // Commit the final partial buffer (the last <30 s the timer didn't
        // reach) and write the populated log out.
        beatLog.flushPending();
        beatLog.writeCsv(cfg.log_path + "/" + stem + "_log.csv");   // was output_path + "/logs/"

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
        std::cout << "Processing (fast: anneal / raw r-peaks / raw templates): "
            << stem << "\n";
        auto jobOpt = post_process_detail::prepareViewerJob(cfg, binFs);
        if (!jobOpt) {
            std::cout << "  prep failed or skipped; not templating.\n";
            continue;
        }
        // shared_ptr so the worker lambda safely co-owns the job; it is
        // joined later (reaped when finished, or drained at shutdown), so the
        // job data always outlives the worker.
        auto job = std::make_shared<post_process_detail::ViewerJob>(std::move(*jobOpt));

        if (!std::filesystem::exists(job->viewerTemplatePath)) {
            std::cerr << "  no template file produced for " << stem
                << "; skipping template marking.\n";
            continue;
        }

        auto done = std::make_shared<std::atomic<bool>>(false);
        std::thread worker;
        if (job->needsFinalize) {
            std::cout << "  squared/absval r-peaks + templates finalizing in background\n";
            worker = std::thread([job, done] {
                // Pure compute + file writes to canonical paths. No Qt here.
                post_process_detail::finalizeViewerJob(*job);
                done->store(true, std::memory_order_release);
                });
        }

        // ---- Stage 3: template marking (runs concurrently with worker) ----
        runTemplateMarking(cfg, job->viewerTemplatePath);

        // Do NOT join here. Park the worker and advance to the next file so
        // the next file loads while this file's squared/absval work finishes.
        if (worker.joinable())
            outstanding.push_back(Outstanding{ std::move(worker), job, done });

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
        std::cout << "\nFinishing " << outstanding.size()
            << " background squared/absval job(s)...\n";
        reap(/*force=*/true);
    }

    std::cout << "\nAll files marked, processed, and templated.\n";
    return 0;
}