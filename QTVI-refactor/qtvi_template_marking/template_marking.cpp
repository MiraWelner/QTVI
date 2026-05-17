#include <QApplication>
#include <QDir>
#include <QFileInfo>
#include <QEventLoop>
#include <QObject>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "config_loader.hpp"
#include "config_entry.hpp"
#include "TemplateFileName.hpp"
#include "TemplateViewerWindow.hpp"

int main(int argc, char* argv[]) {
    QApplication app(argc, argv);

    // Dataset prompt -- matches the noise-marking GUI for consistency.
    std::cout << "\n=== Select Dataset ===\n"
        << "  1) MESA\n  2) BITTIUM\n  3) CHAOS\n"
        << "Enter number (1-3): " << std::flush;

    int choice = 0;
    if (!(std::cin >> choice)) return 1;

    config_entry cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: dataset " << choice << " not in config.csv\n";
        return 1;
    }
    if (!promptForMissingPaths(cfg)) return 1;

    const QString templatePath = QString::fromStdString(cfg.template_path);
    const QString markingPath = QString::fromStdString(cfg.qtvi_marker_path);
    const int     expectedRate = static_cast<int>(cfg.finalSamplingRate);
    const int     expectedBin = static_cast<int>(cfg.bin_length_minutes);

    std::cout << "\nTemplate path: " << templatePath.toStdString()
        << "\nMarking path:  " << markingPath.toStdString()
        << "\nExpected rate: " << expectedRate << " Hz"
        << "\nExpected bin:  " << expectedBin << " min\n\n";

    // Discover candidates.
    QDir templateDir(templatePath);
    if (!templateDir.exists()) {
        std::cerr << "Template directory does not exist: "
            << templatePath.toStdString() << "\n";
        return 1;
    }

    QStringList allFiles = templateDir.entryList({ "*_templates.bin" }, QDir::Files);
    if (allFiles.isEmpty()) {
        std::cerr << "No *_templates.bin files in: "
            << templatePath.toStdString() << "\n";
        QStringList everything = templateDir.entryList(QDir::Files);
        if (!everything.isEmpty()) {
            std::cerr << "Directory contains " << everything.size()
                << " other file(s); first few:\n";
            for (int i = 0; i < std::min(5, (int)everything.size()); ++i)
                std::cerr << "  " << everything[i].toStdString() << "\n";
        }
        return 1;
    }

    // Filter by filename-encoded rate / bin length.
    struct Subject { QString file; QString id; };
    std::vector<Subject> subjects;
    int skippedRate = 0, skippedBin = 0, skippedParse = 0;

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

    // template_marking owns its output tree -- create it before opening the
    // viewer so the first save doesn't fail on a missing directory.
    QDir().mkpath(markingPath);

    TemplateViewerWindow viewer;
    viewer.show();

    const int total = static_cast<int>(subjects.size());
    for (int fi = 0; fi < total; ++fi) {
        const auto& s = subjects[fi];

        std::cout << "Showing " << s.id.toStdString()
            << " | " << (fi + 1) << " of " << total
            << " | Remaining: " << (total - fi - 1) << "\n";
        std::cout.flush();

        QString templateFile = templatePath + "/" + s.file;
        viewer.loadSubject(templateFile, markingPath, s.id);

        QEventLoop loop;
        QObject::connect(&viewer, &TemplateViewerWindow::finished,
            &loop, &QEventLoop::quit);
        loop.exec();
    }
    std::cout << "All subjects processed\n";

    return 0;
}