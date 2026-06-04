/**
 * @file   noise_marking_gui.cpp
 * @brief  Combined entry point for the (formerly two) GUIs.
 *
 *         For every .bin in the dataset's bin folder, this runs the full
 *         pipeline LINEARLY on the main thread, one file at a time:
 *
 *           1. Noise marking      (noise_marking_gui dialog)
 *           2. Processing         (post_process_detail::processOneFile:
 *                                   anneal -> r-peaks -> templates)
 *           3. Template marking   (TemplateViewerWindow)
 *
 *         The old background PostProcessQueue is gone -- step 2 is now a
 *         direct synchronous call between the two dialogs, so a file is
 *         fully finished (marked, processed, templated) before the next
 *         one is opened.
 *
 *         This file replaces BOTH former mains:
 *           - the old noise_marking_gui.cpp main (queue-based), and
 *           - template_marking.cpp main (separate program).
 *         template_marking.cpp should be removed from the build.
 *
 *         The merged target uses the noise-marking config_entry.hpp /
 *         config_loader.* (the superset); the template-marking copies of
 *         those two files are dropped.
 */

#include "post_process.hpp"
#include "config_loader.hpp"
#include "noise_marking_gui\gui_handler.h"
#include "noise_marking_gui\user_annotation_handler.h"
#include "template_marking_gui\parse_data_from_filename.hpp"
#include "template_marking_gui\TemplateViewerWindow.hpp"

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
    std::filesystem::path& outCurrent) {
    auto gui = std::make_unique<noise_marking_gui>();
    gui->setConfig(cfg);
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
    annotation_handler nm(cfg.finalSamplingRate);
    if (markings) {
        for (int i = 0; i < markings->noiseExc.size(); ++i) {
            nm.addSegment(
                static_cast<int>(markings->noiseExc[i].first * cfg.finalSamplingRate),
                static_cast<int>(markings->noiseExc[i].second * cfg.finalSamplingRate),
                markings->data_type[i].toStdString(),
                markings->marking_type[i].toStdString());
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
        displayId, cfg.finalSamplingRate);
    loop.exec();
}


int main(int argc, char* argv[]) {
    QApplication app(argc, argv);
    Theme::apply(app);

    const int dataset_choice = get_dataset_choice();
    auto cfgOpt = load_config(dataset_choice);
    if (!cfgOpt) {
        std::cerr << "Error Loading config.csv\n";
        return 1;
    }
    const config_entry& cfg = *cfgOpt;

    const std::vector<std::filesystem::path> binFiles = load_binfiles(cfg);
    if (binFiles.empty()) {
        std::cerr << "No .bin files in: " << cfg.bin_file_path << "\n";
        return 0;
    }

    for (const std::filesystem::path& binFs : binFiles) {
        const std::string stem = binFs.stem().string();
        std::cout << "\n=== " << stem << " ===\n";

        // ---- Stage 1: noise marking -------------------------------------
        std::cout << "Noise marking: " << binFs.filename().string() << "\n";
        QVector<GenExcStruct> allMarkings;
        std::filesystem::path currentBinFile;
        if (!runNoiseMarking(cfg, binFs, allMarkings, currentBinFile)) {
            std::cout << "  skipped by user; not processing/templating.\n";
            continue;
        }

        // ---- Export markings (input to the anneal step) -----------------
        if (allMarkings.isEmpty()) {
            exportMarkings(cfg, currentBinFile, nullptr);
        }
        else {
            for (const GenExcStruct& m : allMarkings)
                exportMarkings(cfg, std::filesystem::path(m.filePath.toStdString()), &m);
        }

        // ---- Stage 2: process (anneal -> r-peaks -> templates) ----------
        // Synchronous, on this thread. Was previously the background queue.
        std::cout << "Processing (anneal / r-peaks / templates): " << stem << "\n";
        post_process_detail::processOneFile(cfg, binFs);

        // Template file name matches what processOneFile wrote:
        //   <template_path>/<stem>_<binMinutes>_templates.bin
        const std::filesystem::path templateFile =
            std::filesystem::path(cfg.template_path) /
            (stem + "_" +
                std::to_string(static_cast<int>(cfg.bin_length_minutes)) +
                "_templates.bin");

        if (!std::filesystem::exists(templateFile)) {
            std::cerr << "  no template file produced for " << stem
                << " (expected " << templateFile.filename().string()
                << "); skipping template marking.\n";
            continue;
        }

        // ---- Stage 3: template marking ----------------------------------
        runTemplateMarking(cfg, templateFile);
    }

    std::cout << "\nAll files marked, processed, and templated.\n";
    return 0;
}