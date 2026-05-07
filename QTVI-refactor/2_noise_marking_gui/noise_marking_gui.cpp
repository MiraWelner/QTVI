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
#include "noise_manager.h"
#include "bin_handler.hpp"

#include <QtWidgets/QApplication>
#include <QDir>
#include <QFileInfo>

#include <iostream>
#include <memory>

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 1;

    config_entry cfg;
    if (!loadConfig(choice, cfg)) {
        std::cerr << "Error: dataset " << choice << " not in config.csv\n";
        return 1;
    }

    QString markingFolder = QString::fromStdString(cfg.noise_data_path);
    QDir().mkpath(markingFolder);

    QStringList srcFiles = discoverSourceFiles(cfg);
    if (srcFiles.isEmpty()) {
        std::cerr << "No " << cfg.mainExt << " files in: "
            << cfg.original_file_path << "\n";
        return 0;
    }

    // Single pass over every source file: convert, mark, export.
    // (The previous version had an outer "while (madeProgress)" loop, but
    // the source-file list is fixed at startup so a second pass would do
    // nothing -- removed.)
    for (const QString& srcPath : srcFiles) {
        std::cout << "Loading file for noise marking: "
            << QFileInfo(srcPath).fileName().toStdString() << "\n";
        auto binFs = convertToBin(srcPath.toStdString(), cfg);
        if (binFs.empty()) {
            std::cerr << "  conversion failed; skipping\n";
            continue;
        }
        QString binPath = QString::fromStdString(binFs.string());

        
        
        auto gui = std::make_unique<noise_marking_gui>();

        gui->setConfig(cfg);
        gui->setWindowTitle("Marking: " + QFileInfo(binPath).fileName());
        gui->setFileSource(binPath);

        if (gui->exec() != QDialog::Accepted) {
            // User skipped or cancelled.
            continue;
        }

        QVector<GenExcStruct> allMarkings = gui->getAllMarkings();
        QString currentBinFile = gui->getFilePath();

        // exportOne: writes one (CSV, .bin) pair to noiseDataPath. If
        // `markings` is non-null, its segments are added; otherwise the
        // exported file just records "no markings" for this recording.
        auto exportOne = [&](const QString& binFile,
            const GenExcStruct* markings) {
                NoiseManager nm(cfg.finalSamplingRate);
                if (markings) {
                    for (int i = 0; i < markings->noiseExc.size(); ++i) {
                        // Each (start, end) is in seconds; convert to samples
                        // at the marking rate. NB: this uses cfg.finalSamplingRate
                        // for every channel. If per-channel rates ever diverge,
                        // this needs to look up the channel's rate separately.
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
    }

    std::cout << "Starting post-processing...\n";
    int n = processDataset(cfg);
    std::cout << "Post-processing complete. Processed " << n << " files.\n";
    return 0;
}