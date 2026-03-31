#include <QApplication>
#include <QInputDialog>
#include <QDir>
#include <QFileInfo>
#include <QMessageBox>
#include <QEventLoop>
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
        QMessageBox::critical(nullptr, "Error",
            "Could not read config.csv\nTried: " + configPath);
        return 1;
    }

    // Dataset chooser dialog
    QStringList names;
    for (const auto& c : configs) names << c.dataType;

    bool ok = false;
    QString chosen = QInputDialog::getItem(nullptr, "Select Dataset",
        "Dataset:", names, 0, false, &ok);
    if (!ok) return 0;

    auto cfg = findConfig(configs, chosen);
    if (!cfg) {
        QMessageBox::critical(nullptr, "Error", "Dataset not found: " + chosen);
        return 1;
    }

    QMessageBox::information(nullptr, "Config",
        "Template path: " + cfg->templatePath + "\n"
        "Marking path: " + cfg->markingPath);

    // Find all template_info.bin files
    QDir templateDir(cfg->templatePath);
    QStringList binFiles = templateDir.entryList({ "*_template_info.bin" }, QDir::Files);

    if (binFiles.isEmpty()) {
        QMessageBox::warning(nullptr, "No Files",
            "No *_template_info.bin files in:\n" + cfg->templatePath);
        return 1;
    }

    // Ensure marking output dir exists
    QDir().mkpath(cfg->markingPath);

    // Create viewer window
    TemplateViewerWindow viewer;
    viewer.show();

    int loaded = 0;
    for (const QString& binFile : binFiles) {
        // Extract subject ID
        QString stem = QFileInfo(binFile).baseName();
        int pos = stem.indexOf("_template_info");
        if (pos < 0) continue;
        QString subjectId = stem.left(pos);

        // Skip if markings already exist
        QString markingFile = cfg->markingPath + "/" + subjectId + "_template_markings.bin";
        if (QFileInfo::exists(markingFile)) continue;

        // Load subject into viewer
        QString templateFile = cfg->templatePath + "/" + binFile;
        viewer.loadSubject(templateFile, cfg->markingPath, subjectId);
        loaded++;

        // Wait for user to click Finish
        QEventLoop loop;
        QObject::connect(&viewer, &TemplateViewerWindow::finished,
            &loop, &QEventLoop::quit);
        loop.exec();
    }

    if (loaded == 0) {
        QMessageBox::information(nullptr, "Done",
            QString("All %1 subjects already have markings.").arg(binFiles.size()));
    }
    else {
        QMessageBox::information(nullptr, "Done", "All subjects reviewed!");
    }

    return 0;
}