/**
 * @file   noise_marking_gui.cpp
 * @brief  Entry point. For each source file in cfg.originalFilePath:
 *           1. Convert (or pull from cache) into a .bin in cfg.binFilePath
 *           2. Open the marking dialog on that .bin
 *           3. Export marking results (CSV + .bin) into cfg.noiseDataPath
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-18
 */

#include "gui_handler.hpp"
#include "post_process.hpp"
#include "user_annotation_handler.h"
#include "config_loader.hpp"
#include "file_to_bin.hpp"
#include "post_process_queue.hpp"


#include <QtWidgets/QApplication>
#include <QDir>
#include <QFileInfo>
#include <chrono>
#include <thread>
#include <iostream>
#include <memory>

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 1;

    config_entry cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: dataset " << choice << " not in config.csv\n";
        return 1;
    }
    if (!promptForMissingPaths(cfg)) return 1;

    QString markingFolder = QString::fromStdString(cfg.noise_data_path);
    QDir().mkpath(markingFolder);

    QStringList srcFiles = discoverSourceFiles(cfg);
    if (srcFiles.isEmpty()) {
        std::cerr << "No " << cfg.mainExt << " files in: "
            << cfg.input_path << "\n";
        return 0;
    }

    PostProcessQueue postQueue;

    for (const QString& srcPath : srcFiles) {
        std::cout << "Loading file for noise marking: "
            << QFileInfo(srcPath).fileName().toStdString() << "\n";

        auto binFs = make_binfile(srcPath.toStdString(), cfg);
        if (binFs.empty()) {
            std::cerr << "  conversion failed; skipping\n";
            continue;
        }

        // Sequential loop usually skips this; covers the case where a
        // previous session queued this file and it's still in flight.
        if (postQueue.isLocked(binFs)) {
            std::cout << "  waiting for background processing of "
                << binFs.filename().string() << " to finish...\n";
            while (postQueue.isLocked(binFs)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
            }
        }

        QString binPath = QString::fromStdString(binFs.string());

        auto gui = std::make_unique<noise_marking_gui>();
        gui->setConfig(cfg);
        gui->setPostProcessQueue(&postQueue);
        gui->setWindowTitle("Marking: " + QFileInfo(binPath).fileName());
        gui->setFileSource(binPath);

        if (gui->exec() != QDialog::Accepted) continue;

        QVector<GenExcStruct> allMarkings = gui->getAllMarkings();
        QString currentBinFile = gui->getFilePath();

        auto exportOne = [&](const QString& binFile,
            const GenExcStruct* markings) {
                annotation_handler nm(cfg.finalSamplingRate);
                if (markings) {
                    for (int i = 0; i < markings->noiseExc.size(); ++i) {
                        nm.addSegment(
                            (size_t)(markings->noiseExc[i].first * cfg.finalSamplingRate),
                            (size_t)(markings->noiseExc[i].second * cfg.finalSamplingRate),
                            markings->data_type[i].toStdString(),
                            markings->marking_type[i].toStdString());
                    }
                }
                QFileInfo info(binFile);
                QString base = QDir(markingFolder).filePath(
                    info.baseName() + "_noise_markings");
                nm.exportCSV(base.toStdString() + ".csv");
                nm.exportBinary(base.toStdString() + ".bin");
                std::cout << "Saved markings for "
                    << info.fileName().toStdString() << "\n";
            };

        if (allMarkings.isEmpty()) {
            exportOne(currentBinFile, nullptr);
        }
        else {
            for (const GenExcStruct& m : allMarkings) {
                exportOne(m.filePath, &m);
            }
        }

        // Kick off background post-processing. Locks binFs until done;
        // the next iteration's isLocked() check above honors the lock.
        postQueue.enqueue(binFs, cfg);
    }

    std::cout << "All files marked. Waiting for background "
        "post-processing to finish...\n";
    postQueue.drain();
    std::cout << "All done.\n";
    return 0;
}