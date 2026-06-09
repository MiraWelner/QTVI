/**
 * @file   main.cpp
 * @brief  Standalone utility for converting EDF / DAT recordings into the
 *         40-channel uniform .bin format used by the QTVi pipeline.
 *
 *         Walks the input directory configured for the chosen dataset,
 *         calls make_binfile on each source file, and writes the result
 *         to the bin_file_path subfolder of the output tree. Existing
 *         bin files are reused (make_binfile checks for them).
 *
 *         Run with no arguments and follow the prompts:
 *           1) MESA   2) BITTIUM   3) CHAOS
 *
 *         The config.csv must define input_path / output_path for the
 *         dataset, or the prompt will ask for them.
 */

#include "config_loader.hpp"
#include "config_entry.hpp"
#include "file_to_bin.hpp"

#include <QApplication>
#include <QDir>
#include <QDirIterator>
#include <QFileInfo>
#include <QString>
#include <QStringList>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    // QApplication is needed because promptForMissingPaths uses
    // QFileDialog. If you ever want a fully-headless mode, the
    // dataset config + paths would have to be passed via command-line
    // args instead.
    QApplication app(argc, argv);

    // ---- Dataset prompt --------------------------------------------------
    std::cout << "\nConvert to uniform .bin\n"
        << "  1) MESA\n"
        << "  2) BITTIUM\n"
        << "  3) CHAOS\n"
        << "Enter number (1-3): " << std::flush;

    int choice = 0;
    if (!(std::cin >> choice) || choice < 1 || choice > 3) {
        std::cerr << "Invalid choice.\n";
        return 1;
    }

    // ---- Load config -----------------------------------------------------
    config_entry cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: dataset " << choice << " not found in config.csv\n";
        return 1;
    }
    if (!promptForMissingPaths(cfg)) {
        std::cerr << "Error: input/output path not provided\n";
        return 1;
    }

    // ---- Discover source files ------------------------------------------
    QStringList files = discoverSourceFiles(cfg);
    if (files.isEmpty()) {
        std::cerr << "No source files found in " << cfg.input_path
            << " with extension " << cfg.dataset_type << "\n";
        return 1;
    }
    std::cout << "Found " << files.size() << " source file(s).\n\n";

    // ---- Convert one at a time ------------------------------------------
    for (int i = 0; i < files.size(); ++i) {
        const fs::path srcPath = files[i].toStdString();
        std::cout << "[" << (i + 1) << "/" << files.size() << "] "
            << srcPath.filename().string() << "\n";
        std::cout.flush();
        const fs::path outPath = make_binfile(srcPath, cfg);
    }
}