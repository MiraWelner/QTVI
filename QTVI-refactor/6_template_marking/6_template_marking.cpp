#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QEventLoop>
#include <iostream>
#include <string>
#include "ConfigReader.hpp"
#include "TemplateFileName.hpp"
#include "TemplateViewerWindow.h"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    QString configPath = QCoreApplication::applicationDirPath() + "/config.csv";
    if (!QFileInfo::exists(configPath)) configPath = "config.csv";

    auto configs = readConfig(configPath);
    if (configs.empty()) {
        std::cerr << "Error: Could not read config.csv (tried: "
            << configPath.toStdString() << ")\n";
        return 1;
    }

    // Terminal dataset selection
    std::cout << "\n=== Select Dataset ===\n";
    for (size_t i = 0; i < configs.size(); ++i)
        std::cout << "  " << (i + 1) << ") " << configs[i].dataType.toStdString() << "\n";
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
        << "\nMarking path:  " << cfg.markingPath.toStdString()
        << "\nExpected rate: " << cfg.upsampledRate << " Hz\n\n";

    // Discover candidates
    QDir templateDir(cfg.templatePath);

    if (!templateDir.exists()) {
        std::cerr << "Template directory does not exist: "
            << cfg.templatePath.toStdString() << "\n";
        return 1;
    }

    QStringList allFiles = templateDir.entryList({ "*_templates.bin" }, QDir::Files);

    if (allFiles.isEmpty()) {
        std::cerr << "No *_templates.bin files in: "
            << cfg.templatePath.toStdString() << "\n";
        // Show what IS there so the mismatch is obvious
        QStringList everything = templateDir.entryList(QDir::Files);
        if (!everything.isEmpty()) {
            std::cerr << "Directory contains " << everything.size()
                << " other file(s); first few:\n";
            for (int i = 0; i < std::min(5, (int)everything.size()); ++i)
                std::cerr << "  " << everything[i].toStdString() << "\n";
        }
        return 1;
    }

    // Filter by filename-encoded rate / bin length
    struct Subject { QString file; QString id; };
    std::vector<Subject> subjects;
    int skippedRate = 0, skippedBin = 0, skippedParse = 0;

    const int expectedRate = static_cast<int>(cfg.upsampledRate);
    const int expectedBin = cfg.binMinutes;

    for (const QString& f : allFiles) {
        QString stem = QFileInfo(f).completeBaseName();   // strip ".bin"
        int suffixPos = stem.indexOf("_templates");
        if (suffixPos < 0) { ++skippedParse; continue; }
        QString nameStem = stem.left(suffixPos);

        auto parsed = parseTemplateFileName(nameStem);
        if (!parsed) { ++skippedParse; continue; }

        if (expectedRate > 0 && parsed->rateHz != expectedRate) {
            std::cerr << "  skip (rate " << parsed->rateHz
                << " != " << expectedRate << "): "
                << f.toStdString() << "\n";
            ++skippedRate;
            continue;
        }

        if (expectedBin > 0 && parsed->binMinutes != expectedBin) {
            std::cerr << "  skip (bin " << parsed->binMinutes
                << " != " << expectedBin << "): "
                << f.toStdString() << "\n";
            ++skippedBin;
            continue;
        }

        subjects.push_back({ f, parsed->subjectId });
    }

    if (skippedParse + skippedRate + skippedBin > 0) {
        std::cout << "\nSkipped " << skippedParse << " unparseable, "
            << skippedRate << " rate-mismatched, "
            << skippedBin << " bin-mismatched files.\n\n";
    }

    if (subjects.empty()) {
        std::cerr << "No matching files after filtering.\n";
        return 1;
    }

    QDir().mkpath(cfg.markingPath);

    TemplateViewerWindow viewer;
    viewer.show();

    int total = static_cast<int>(subjects.size());
    for (int fi = 0; fi < total; ++fi) {
        const auto& s = subjects[fi];

        std::cout << "Showing " << s.id.toStdString()
            << " | " << (fi + 1) << " of " << total
            << " | Remaining: " << (total - fi - 1) << "\n";
        std::cout.flush();

        QString templateFile = cfg.templatePath + "/" + s.file;
        viewer.loadSubject(templateFile, cfg.markingPath, s.id);

        QEventLoop loop;
        QObject::connect(&viewer, &TemplateViewerWindow::finished,
            &loop, &QEventLoop::quit);
        loop.exec();
    }
    std::cout << "All subjects processed\n";

    return 0;
}