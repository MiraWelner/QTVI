/**
 * @file   noise_marking_gui.cpp
 * @brief  Entry point of the noise marking program. Makes a GUI where you can mark the noise you see in the
 *         ECG/PPG, and outputs a csv and a .bin file
 *
 * Output .bin format:
 *   Header:
 *     [uint64]  count — number of annotation segments
 *
 *   Followed by count rows of 6 doubles each (48 bytes per row):
 *     [0] startSample   — first sample index of the annotation
 *     [1] endSample     — last sample index of the annotation
 *     [2] startSec      — start time in seconds (startSample / sampleRate)
 *     [3] endSec        — end time in seconds (endSample / sampleRate)
 *     [4] labelId       — signal type: 0=unknown, 1=PPG, 2=ECG, 3=BOTH
 *     [5] typeId        — marking type:
 *                            0=unknown, 1=Noise/Artifact, 2=AF, 3=SVT, 4=VT,
 *                            5=PVC, 6=PAC, 7=Benign Arrhythmia, 8=Significant Arrhythmia
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 * @date   2026-03-18
 */
#include "gui_handler.h"
#include "NoiseManager.h"
#include <QtWidgets/QApplication>
#include <QDir>
#include <QDirIterator>
#include <QFile>
#include <QFileInfo>
#include <iostream>
#include <memory>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>

static const std::string CONFIG_PATH = "config.csv";
static const int SAMPLING_RATE = 2000;
static const int SLEEP_EPOCH_SIZE = 30;

struct input_bin_header {
    //the intput bin file is loaded and contains data about how 
    //long the signals are
    uint64_t ecg1_signal_length;
    uint64_t ecg2_signal_length;
    uint64_t ecg3_signal_length;
    uint64_t ppg_signal_length;
    uint64_t sleep_signal_length;
};

struct config_csv_data {
    std::string outputPath;
    std::string markingPath;
};


static std::vector<std::string> parse_csv_row(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (size_t i = 0; i < line.length(); ++i) {
        if (line[i] == ',') { fields.push_back(cur); cur = ""; }
        else cur += line[i];
    }
    fields.push_back(cur);
    for (auto& f : fields) {
        size_t first = f.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) f = "";
        else {
            size_t last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
    }
    return fields;
}

static bool load_config(int data_type, config_csv_data& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) return false;

    std::string target = (data_type == 1) ? "MESA" : (data_type == 2) ? "BITTIUM" : (data_type == 3) ? "CHAOS" : "";
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        if (row.size() < 6) continue;

        std::string rType = row[0];
        std::transform(rType.begin(), rType.end(), rType.begin(), ::toupper);

        if (rType == target) {
            out.outputPath = row[4];
            out.markingPath = row[5];
            return true;
        }
    }
    return false;
}

static input_bin_header read_bin_header(const QString& path) {
    /**
     * @brief  Reads the 64-byte binary header from a .bin file produced by file_to_bin.
     *
     * @param  path  Path to the .bin file
     * @return       Populated BinHeader struct (zeroed if the file can't be opened)
     */
    input_bin_header header;
    std::memset(&header, 0, sizeof(header));
    QFile file(path);
    if (!file.open(QIODevice::ReadOnly)) return header;
    file.read(reinterpret_cast<char*>(&header), sizeof(header));
    return header;
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    std::cout << "NOISE MARKING\n";
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 1;

    config_csv_data cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: Could not find configuration for selection " << choice << " in config.csv" << std::endl;
        return 1;
    }

    QString binFolder = QString::fromStdString(cfg.outputPath);
    QString markingFolder = QString::fromStdString(cfg.markingPath);

    QDir().mkpath(markingFolder);

    // Find all .bin files (excluding noise marking output files)
    QDirIterator it(binFolder, QStringList() << "*.bin", QDir::Files, QDirIterator::Subdirectories);
    QStringList binFiles;
    while (it.hasNext()) {
        QString p = it.next();
        if (!p.contains("_noise_markings.bin")) binFiles << p;
    }
    binFiles.sort();

    if (binFiles.isEmpty()) {
        std::cerr << "No .bin files found in: " << cfg.outputPath << std::endl;
        return 0;
    }

    for (const QString& binPath : binFiles) {
        QFileInfo info(binPath);
        input_bin_header header = read_bin_header(binPath);

        auto gui = std::make_unique<noise_marking_gui>();
        gui->setWindowTitle("Marking: " + info.fileName());
        gui->setFileSource(binPath);

        if (gui->exec() != QDialog::Accepted) continue;

        const GenExcStruct markings = gui->getMarkings();
        NoiseManager noiseHandler(SAMPLING_RATE);

        for (int i = 0; i < markings.noiseExc.size(); ++i) {
            noiseHandler.addSegment(
                static_cast<size_t>(markings.noiseExc[i].first * SAMPLING_RATE),
                static_cast<size_t>(markings.noiseExc[i].second * SAMPLING_RATE),
                markings.data_type[i].toStdString(),
                markings.marking_type[i].toStdString()
            );
        }

        QString outBase = QDir(markingFolder).filePath(info.baseName() + "_noise_markings");
        noiseHandler.exportCSV(outBase.toStdString() + ".csv");
        noiseHandler.exportBinary(outBase.toStdString() + ".bin");
        std::cout << "Saved markings for " << info.fileName().toStdString() << std::endl;
    }

    std::cout << "Marking Complete." << std::endl;
    return 0;
}