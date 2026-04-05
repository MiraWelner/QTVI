#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QEventLoop>
#include <iostream>
#include <string>
#include "ConfigReader.hpp"
#include "TemplateViewerWindow.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    // Find config.csv next to the executable or in working dir
    QString configPath = QCoreApplication::applicationDirPath() + "/config.csv";
    if (!QFileInfo::exists(configPath)) {
        configPath = "config.csv";
    }

    auto configs = readConfig(configPath);
    if (configs.empty()) {
        std::cerr << "Error: Could not read config.csv (tried: "
            << configPath.toStdString() << ")\n";
        return 1;
    }

    // Terminal dataset selection
    std::cout << "\n=== Select Dataset ===\n";
    for (size_t i = 0; i < configs.size(); ++i) {
        std::cout << "  " << (i + 1) << ") " << configs[i].dataType.toStdString() << "\n";
    }
    std::cout << "\nEnter number (1-" << configs.size() << "): ";
    std::cout.flush();

    int choice = 0;
    std::string line;
    if (!std::getline(std::cin, line)) return 0;

    try { choice = std::stoi(line); }
    catch (...) { choice = 0; }

    if (choice < 1 || choice >(int)configs.size()) {
        std::cerr << "Invalid selection.\n";
        return 1;
    }

    const auto& cfg = configs[choice - 1];
    std::cout << "\nTemplate path: " << cfg.templatePath.toStdString()
        << "\nMarking path:  " << cfg.markingPath.toStdString() << "\n\n";

    // Find all template_info.bin files
    QDir templateDir(cfg.templatePath);
    QStringList binFiles = templateDir.entryList({ "*_template_info.bin" }, QDir::Files);

    if (binFiles.isEmpty()) {
        std::cerr << "No *_template_info.bin files in: "
            << cfg.templatePath.toStdString() << "\n";
        return 1;
    }

    // Ensure marking output dir exists
    QDir().mkpath(cfg.markingPath);

    // Create viewer window
    TemplateViewerWindow viewer;
    viewer.show();

    int loaded = 0;
    int total = binFiles.size();
    int skipped = 0;

    for (int fi = 0; fi < total; ++fi) {
        const QString& binFile = binFiles[fi];

        // Extract subject ID
        QString stem = QFileInfo(binFile).baseName();
        int pos = stem.indexOf("_template_info");
        if (pos < 0) continue;
        QString subjectId = stem.left(pos);

        ++loaded;
        std::cout << "Showing " << subjectId.toStdString()
            << " | " << (fi + 1) << " of " << total
            << " | Remaining: " << (total - fi - 1) << "\n";
        std::cout.flush();

        // Load subject into viewer
        QString templateFile = cfg.templatePath + "/" + binFile;
        viewer.loadSubject(templateFile, cfg.markingPath, subjectId);

        // Wait for user to click Finish
        QEventLoop loop;
        QObject::connect(&viewer, &TemplateViewerWindow::finished,
            &loop, &QEventLoop::quit);
        loop.exec();
    }
    std::cout << "All subjects proccessed\n";

    return 0;
}