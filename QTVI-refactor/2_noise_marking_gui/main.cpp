#include "noise_marking_gui.h"
#include "NoiseManager.hpp"
#include <QtWidgets/QApplication>
#include <QDir>
#include <QDirIterator>
#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QDebug>
#include <QStringList>
#include <QVector>
#include <iostream> 
#include <memory>
#include <cstring>

// Structure to match the 64-byte Universal Header from edf_to_bin.cpp
struct BinHeader {
    double ecgRate;      // Sampling rate for ECG signals
    double ppgRate;      // Sampling rate for PPG signal
    double epochSize;    // Epoch size (typically 30.0)
    uint64_t size1;      // Number of samples in ECG signal 1
    uint64_t size2;      // Number of samples in ECG signal 2
    uint64_t size3;      // Number of samples in ECG signal 3
    uint64_t sizeP;      // Number of samples in PPG signal
    uint64_t sizeS;      // Number of sleep stage values
};

// Structure to match the General config.csv row
struct DatasetConfig {
    QString binFolder;      // Column 4 (Index 4)
    QString markingFolder;  // Column 5 (Index 5)
};

// --- HELPERS ---

// Helper to handle commas inside quotes in CSV
QStringList parseCsvLine(const QString& line) {
    QStringList fields;
    QString cur;
    bool inQuotes = false;
    for (int i = 0; i < line.length(); ++i) {
        QChar c = line[i];
        if (c == '\"') inQuotes = !inQuotes;
        else if (c == ',' && !inQuotes) {
            fields << cur.trimmed().remove('\"');
            cur = "";
        }
        else cur += c;
    }
    fields << cur.trimmed().remove('\"');
    return fields;
}

DatasetConfig loadDatasetConfig(int choice) {
    DatasetConfig config;
    QFile file("config.csv");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        qCritical() << "Error: Could not find config.csv in the application directory.";
        return config;
    }

    QString target = (choice == 1) ? "MESA" : (choice == 2) ? "BITTIUM" : (choice == 3) ? "CHAOS" : "";
    QTextStream in(&file);
    in.readLine(); // Skip header row

    while (!in.atEnd()) {
        QString line = in.readLine();
        if (line.trimmed().isEmpty()) continue;

        QStringList fields = parseCsvLine(line);
        if (fields.size() > 5 && fields[0].toUpper() == target.toUpper()) {
            config.binFolder = fields[4];      // Column 5: Where the BINs are
            config.markingFolder = fields[5];  // Column 6: Where markings go
            break;
        }
    }
    return config;
}

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    // Step 1: User Selection
    std::cout << "Select dataset for Noise Marking:\n";
    std::cout << "1: MESA\n2: Bittium\n3: Chaos\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return 0;

    // Step 2: Load Paths from the same config.csv used by the converter
    DatasetConfig cfg = loadDatasetConfig(choice);

    if (cfg.binFolder.isEmpty() || cfg.markingFolder.isEmpty()) {
        qCritical() << "Error: Could not find paths for selection " << choice << " in config.csv";
        return -1;
    }

    QDir().mkpath(cfg.markingFolder);

    // Step 3: Find BIN files
    QDirIterator it(cfg.binFolder, QStringList() << "*.bin", QDir::Files, QDirIterator::Subdirectories);
    QStringList binFiles;
    while (it.hasNext()) {
        QString p = it.next();
        if (!p.contains("_noise_markings.bin")) binFiles << p;
    }
    binFiles.sort();

    if (binFiles.isEmpty()) {
        qWarning() << "No .bin files found in: " << cfg.binFolder;
        return 0;
    }

    // Step 4: Process Files
    for (const QString& binPath : binFiles) {
        QFileInfo info(binPath);

        // Read the Universal Header (64 bytes)
        BinHeader header;
        // Initialize with default values
        std::memset(&header, 0, sizeof(BinHeader));
        header.ecgRate = 256.0;
        header.ppgRate = 256.0;
        header.epochSize = 30.0;

        QFile file(binPath);
        if (file.open(QIODevice::ReadOnly)) {
            qint64 bytesRead = file.read(reinterpret_cast<char*>(&header), sizeof(BinHeader));
            if (bytesRead != sizeof(BinHeader)) {
                qWarning() << "Warning: Could not read complete header from" << binPath;
                qWarning() << "Expected" << sizeof(BinHeader) << "bytes, got" << bytesRead;
            }
            file.close();

            // Validate header values
            if (header.ecgRate <= 0.0 || header.ecgRate > 10000.0) {
                qWarning() << "Warning: Invalid ECG rate in header:" << header.ecgRate << "- using default 256.0";
                header.ecgRate = 256.0;
            }
            if (header.ppgRate <= 0.0 || header.ppgRate > 10000.0) {
                qWarning() << "Warning: Invalid PPG rate in header:" << header.ppgRate << "- using default 256.0";
                header.ppgRate = 256.0;
            }

            // Debug output
            qDebug() << "File:" << info.fileName();
            qDebug() << "  ECG Rate:" << header.ecgRate;
            qDebug() << "  PPG Rate:" << header.ppgRate;
            qDebug() << "  Epoch Size:" << header.epochSize;
            qDebug() << "  Signal 1 samples:" << header.size1;
            qDebug() << "  Signal 2 samples:" << header.size2;
            qDebug() << "  Signal 3 samples:" << header.size3;
            qDebug() << "  PPG samples:" << header.sizeP;
            qDebug() << "  Sleep stages:" << header.sizeS;
        }
        else {
            qWarning() << "Warning: Could not open file to read header:" << binPath;
        }


        auto gui = std::make_unique<noise_marking_gui>();
        gui->setWindowTitle("Marking: " + info.fileName());
        gui->setFileSource(binPath);

        if (gui->exec() == QDialog::Accepted) {
            const GenExcStruct markings = gui->getMarkings();
            NoiseManager noiseHandler(header.ecgRate);

            for (int i = 0; i < markings.noiseExc.size(); ++i) {
                qDebug() << "1";
                const size_t start = static_cast<size_t>(markings.noiseExc[i].first * header.ecgRate);
                qDebug() << "1";
                const size_t end = static_cast<size_t>(markings.noiseExc[i].second * header.ecgRate);
                noiseHandler.addSegment(start, end,
                    markings.data_type[i].toStdString(),
                    markings.marking_type[i].toStdString());
            }

            // Save results to the folder specified in Column 6 (Index 5)
            QString outBase = QDir(cfg.markingFolder).filePath(info.baseName() + "_noise_markings");
            noiseHandler.exportCSV(outBase.toStdString() + ".csv");
            noiseHandler.exportBinary(outBase.toStdString() + ".bin");

            qInfo() << "Saved markings for" << info.fileName() << "to" << cfg.markingFolder;
        }
    }

    return 0;
}
