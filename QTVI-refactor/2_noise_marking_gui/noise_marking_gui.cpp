/**
 * @file   noise_marking_gui.cpp
 * @brief  Entry point � iterates .bin files, opens the marking GUI for each,
 *         and exports annotations as CSV + binary.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-18
 */
#include "gui_handler.hpp"
#include "NoiseManager.h"

#include <QtWidgets/QApplication>
#include <QDir>
#include <QDirIterator>
#include <QFile>
#include <QFileInfo>
#include <QSet>

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <memory>

static const std::string CONFIG_PATH = "config.csv";
static const int         SAMPLING_RATE = 1000;

// ============================================================================
// Config helpers
// ============================================================================

struct ConfigEntry {
    std::string outputPath;
    std::string markingPath;
};

static std::vector<std::string> parseCsvRow(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (char c : line) {
        if (c == ',') { fields.push_back(cur); cur.clear(); }
        else          cur += c;
    }
    fields.push_back(cur);

    for (auto& f : fields) {
        size_t first = f.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) { f.clear(); continue; }
        size_t last = f.find_last_not_of(" \t\r\n");
        f = f.substr(first, last - first + 1);
    }
    return fields;
}

static bool loadConfig(int dataType, ConfigEntry& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) return false;

    std::string target = (dataType == 1) ? "MESA"
        : (dataType == 2) ? "BITTIUM"
        : (dataType == 3) ? "CHAOS" : "";

    std::string line;
    std::getline(file, line);  // skip header

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        auto row = parseCsvRow(line);
        if (row.size() < 6) continue;

        std::string rowType = row[0];
        std::transform(rowType.begin(), rowType.end(), rowType.begin(), ::toupper);

        if (rowType == target) {
            out.outputPath = row[4];
            out.markingPath = row[5];
            return true;
        }
    }
    return false;
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    std::cout << "NOISE MARKING\n"
        << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 1;

    ConfigEntry cfg;
    if (!loadConfig(choice, cfg)) {
        std::cerr << "Error: Could not find configuration for selection "
            << choice << " in config.csv\n";
        return 1;
    }

    QString binFolder = QString::fromStdString(cfg.outputPath);
    QString markingFolder = QString::fromStdString(cfg.markingPath);
    QDir().mkpath(markingFolder);

    // Collect .bin files (skip noise-marking outputs)
    QDirIterator it(binFolder, { "*.bin" }, QDir::Files, QDirIterator::Subdirectories);
    QStringList binFiles;
    while (it.hasNext()) {
        QString p = it.next();
        if (!p.contains("_noise_markings.bin"))
            binFiles << p;
    }
    binFiles.sort();

    if (binFiles.isEmpty()) {
        std::cerr << "No .bin files found in: " << cfg.outputPath << "\n";
        return 0;
    }

    QSet<QString> savedFiles;
    QSet<QString> skippedFiles;   // explicitly skipped via Skip button

    bool madeProgress = true;
    while (madeProgress) {
        madeProgress = false;

        for (const QString& binPath : binFiles) {
            if (savedFiles.contains(binPath) || skippedFiles.contains(binPath))
                continue;

            QFileInfo info(binPath);

            auto gui = std::make_unique<noise_marking_gui>();
            gui->setWindowTitle("Marking: " + info.fileName());
            gui->setFileSource(binPath);

            if (gui->exec() != QDialog::Accepted) {
                skippedFiles.insert(binPath);
                madeProgress = true;
                continue;
            }

            QVector<GenExcStruct> allMarkings = gui->getAllMarkings();
            QString currentFile = gui->getFilePath();

            // If the user browsed to a different file, the loop's file
            // was never actually worked on — don't mark it as done.
            // Only mark files that were explicitly exported.

            if (allMarkings.isEmpty()) {
                // User hit Save without marking anything — export empty
                // file for whatever is currently loaded.
                NoiseManager noiseHandler(SAMPLING_RATE);
                QFileInfo curInfo(currentFile);
                QString outBase = QDir(markingFolder).filePath(
                    curInfo.baseName() + "_noise_markings");
                noiseHandler.exportCSV(outBase.toStdString() + ".csv");
                noiseHandler.exportBinary(outBase.toStdString() + ".bin");
                std::cout << "Saved markings for "
                    << curInfo.fileName().toStdString()
                    << std::endl;
                savedFiles.insert(currentFile);
            }
            else {
                for (const GenExcStruct& markings : allMarkings) {
                    NoiseManager noiseHandler(SAMPLING_RATE);
                    for (int i = 0; i < markings.noiseExc.size(); ++i) {
                        noiseHandler.addSegment(
                            static_cast<size_t>(markings.noiseExc[i].first * SAMPLING_RATE),
                            static_cast<size_t>(markings.noiseExc[i].second * SAMPLING_RATE),
                            markings.data_type[i].toStdString(),
                            markings.marking_type[i].toStdString()
                        );
                    }

                    QFileInfo markedFileInfo(markings.filePath);
                    QString outBase = QDir(markingFolder).filePath(
                        markedFileInfo.baseName() + "_noise_markings");
                    noiseHandler.exportCSV(outBase.toStdString() + ".csv");
                    noiseHandler.exportBinary(outBase.toStdString() + ".bin");
                    std::cout << "Saved markings for "
                        << markedFileInfo.fileName().toStdString()
                        << std::endl;

                    savedFiles.insert(markings.filePath);
                }
            }

            madeProgress = true;
        }
    }

    std::cout << "Marking Complete.\n";
    return 0;
}