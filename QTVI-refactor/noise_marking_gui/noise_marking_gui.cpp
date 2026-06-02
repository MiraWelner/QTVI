#include "gui_handler.h"
#include "post_process.hpp"
#include "user_annotation_handler.h"
#include "config_loader.hpp"
#include "post_process_queue.hpp"

#include <QtWidgets/QApplication>
#include <QDir>
#include <QFileInfo>
#include <QDirIterator>
#include <filesystem>
#include <chrono>
#include <thread>
#include <iostream>
#include <memory>
#include <theme/theme.h>


int get_dataset_choice() {
    /*
        Ask the user to select a MESA, Bittium, or CHAOS dataset
    */
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return -1;
    while (choice < 1 || choice > 3) {
        std::cout << "Invalid choice. Please enter 1, 2, or 3: ";
        if (!(std::cin >> choice)) return -1;
	}
    return choice;
}


std::vector<std::filesystem::path> load_binfiles(config_entry cfg) {
    std::vector<std::filesystem::path> binFiles;
    for (const auto& entry : std::filesystem::directory_iterator(cfg.bin_file_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin")
            binFiles.push_back(entry.path());
    }
    std::sort(binFiles.begin(), binFiles.end());

	return binFiles;
}


int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    Theme::apply(app);

    int dataset_choice = get_dataset_choice();
    auto cfgOpt = load_config(dataset_choice);
    if (!cfgOpt) {
        std::cerr << "Error Loading config.csv";
        return 1;
    }
    const config_entry& cfg = *cfgOpt;

    PostProcessQueue postQueue;

    std::vector<std::filesystem::path> binFiles = load_binfiles(cfg);

    if (!binFiles.size()) {
        std::cerr << "No .bin files in: " << cfg.bin_file_path << "\n";
        return 0;
    }

    for (const std::filesystem::path& binFs : binFiles) {
        std::cout << "Loading file for noise marking: "
            << binFs.filename().string() << "\n";

        if (postQueue.isLocked(binFs)) {
            std::cout << "  waiting for background processing of "
                << binFs.filename().string() << " to finish...\n";
            while (postQueue.isLocked(binFs)) {
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
            }
        }

        auto gui = std::make_unique<noise_marking_gui>();
        gui->setConfig(cfg);
        gui->setPostProcessQueue(&postQueue);
        QScreen* screen = QGuiApplication::primaryScreen();
        QRect available = screen->availableGeometry();
        gui->setMaximumHeight(available.height() - 40);
        gui->setWindowTitle(QString::fromStdString(
            "Marking: " + binFs.filename().string()));
        gui->setFileSource(QString::fromStdString(binFs.string()));

        if (gui->exec() != QDialog::Accepted) continue;

        QVector<GenExcStruct> allMarkings = gui->getAllMarkings();
        std::filesystem::path currentBinFile = gui->getFilePath().toStdString();

        auto exportOne = [&](const std::filesystem::path& binFile,
            const GenExcStruct* markings) {
                annotation_handler nm(cfg.finalSamplingRate);
                if (markings) {
                    for (int i = 0; i < markings->noiseExc.size(); ++i) {
                        nm.addSegment(
                            (int)(markings->noiseExc[i].first * cfg.finalSamplingRate),
                            (int)(markings->noiseExc[i].second * cfg.finalSamplingRate),
                            markings->data_type[i].toStdString(),
                            markings->marking_type[i].toStdString());
                    }
                }
                std::filesystem::path base = std::filesystem::path(cfg.noise_data_path)
                    / (binFile.stem().string() + "_noise_markings");
                nm.exportCSV(base.string() + ".csv");
                nm.exportBinary(base.string() + ".bin");
                std::cout << "Saved markings for "
                    << binFile.filename().string() << "\n";
            };

        if (allMarkings.isEmpty()) {
            exportOne(currentBinFile, nullptr);
        }
        else {
            for (const GenExcStruct& m : allMarkings) {
                exportOne(std::filesystem::path(m.filePath.toStdString()), &m);
            }
        }

        postQueue.enqueue(binFs, cfg);
    }


    std::cout << "All files marked. Waiting for background "
        "post-processing to finish...\n";
    postQueue.drain();
    std::cout << "All done.\n";
    return 0;
}