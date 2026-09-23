/**
 * @file   main.cpp
 * @brief  Entry point for the electronic record analysis and templating pipeline
 */

#include "analysis_job.hpp"
#include "config_file_handling/config_loader.hpp"
#include "noise_marking_gui/gui_handler.h"
#include "noise_marking_gui/user_annotation_handler.h"
#include "template_marking_gui/parse_data_from_filename.hpp"
#include "logging/user_mark_log.hpp"
#include "template_viewer/template_viewer.hpp"

#include <QtWidgets/QApplication>
#include <QGuiApplication>
#include <QScreen>
#include <QEventLoop>
#include <QObject>
#include <QFileInfo>

#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <thread>
#include <QtWidgets/QMessageBox>
#include <vector>

#include <theme/theme.h>

static int get_dataset_choice() {
	// Prompt the user to select a MESA, Bittium, Chaos, or SHHS dataset
    static constexpr int n_valid_datasets = 4;
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\n 4: SHHS\nChoice: ";
    int choice;
    if (!(std::cin >> choice)) return -1;
    while (choice < 1 || choice > n_valid_datasets) {
        std::cout << "Please enter 1-" << n_valid_datasets << ": ";
        if (!(std::cin >> choice)) return -1;
    }
    return choice;
}

static std::string get_initials() {
    //ask the user's initials such that the log file is saved in a folder with their initials
    std::cout << "Enter your initials: ";
    std::string raw;
    std::getline(std::cin >> std::ws, raw);
    std::string out;
    for (char c : raw)
        if (std::isalnum(static_cast<unsigned char>(c)))
            out += static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    if (out.empty()) out = "";
    return out;
}

static std::vector<std::filesystem::path> load_binfiles(const config_entry& cfg) {
    std::vector<std::filesystem::path> binFiles;
    for (const auto& entry :
        std::filesystem::directory_iterator(cfg.bin_file_path)) {
        if (entry.is_regular_file() && entry.path().extension() == ".bin")
            binFiles.push_back(entry.path());
    }
    std::sort(binFiles.begin(), binFiles.end());
    return binFiles;
}

static bool runNoiseMarking(const config_entry& cfg, const std::filesystem::path& binFs, QVector<GenExcStruct>& outAll, std::filesystem::path& outCurrent, beat_log& beatLog, bool& out_ecg1_inv, bool& out_ecg2_inv, bool& outEcg3Inverted) {
    // Launch the GUI to do the noise marking. One important thing that takes place is that the gui object (a noise_marking_gui) has
    // an invertedForSignal attribute for each channel
    auto gui = std::make_unique<noise_marking_gui>();
    gui->set_params_to_config_defaults(cfg);
    gui->setBeatLog(&beatLog);

    QScreen* screen = QGuiApplication::primaryScreen();
    if (screen) {
        const QRect available = screen->availableGeometry();
        gui->setMaximumHeight(available.height() - 40);
    }
    gui->setWindowTitle(
        QString::fromStdString("Marking: " + binFs.filename().string()));
    gui->setFileSource(QString::fromStdString(binFs.string()));

    if (gui->exec() != QDialog::Accepted)
        return false;

    outAll = gui->getAllMarkings();
    outCurrent = gui->getFilePath().toStdString();

    out_ecg1_inv = gui->invertedForSignal("ECG1");
    out_ecg2_inv = gui->invertedForSignal("ECG2");
    outEcg3Inverted = gui->invertedForSignal("ECG3");

    return true;
}

static void exportMarkings(const config_entry& cfg, const std::filesystem::path& binFile, const GenExcStruct* markings) {
    // Export the recorded markings to <noise_data_path>/<stem>_noise_markings.{csv,bin}. The .bin is the file the anneal step (processOneFile) reads back in.
    annotation_handler nm;
    auto rateForLabel = [&](const QString& label) -> double {
        if (label == "PPG")      return cfg.ppg_upsample_rate;
        if (label == "ABP")      return cfg.abp_upsample_rate;
        if (label == "ACCEL")    return cfg.accel_upsample_rate;
        if (label == "ART")      return cfg.art_upsample_rate;
        if (label == "ART_PULM") return cfg.art_pulm_upsample_rate;
        return cfg.ecg_upsample_rate;   // ECG1/ECG2/ECG3
        };

    if (markings) {
        for (int i = 0; i < markings->noiseExc.size(); ++i) {
            const double sr = rateForLabel(markings->data_type[i]);
            nm.addSegment(
                static_cast<int>(markings->noiseExc[i].first * sr),
                static_cast<int>(markings->noiseExc[i].second * sr),
                markings->data_type[i].toStdString(),
                markings->marking_type[i].toStdString(), sr,
                markings->threshold[i], markings->blanking[i]);
        }
    }
    const std::filesystem::path base =  std::filesystem::path(cfg.noise_data_path) / (binFile.stem().string() + "_noise_markings");
    nm.exportCSV(base.string() + ".csv");
    nm.exportBinary(base.string() + ".bin");
    std::cout << "Saved Noise Markings \n";
}

static void runTemplateMarking(const config_entry& cfg, std::shared_ptr<analysis_job::AnalysisJob> job, const QString& fileId, std::vector<analysis_job::BankSnapshot>& outBanks) {
        //Launch the template marking GUI
        TemplateViewerWindow viewer;
        viewer.setBoundaryTrainingDir(QString::fromStdString(cfg.training_log));
        viewer.set_vcg_output_dir(QString::fromStdString(cfg.vcg_output));
        viewer.setNormOutputDir(QString::fromStdString(cfg.template_path));
        QEventLoop loop;
        QObject::connect(&viewer, &TemplateViewerWindow::finished,  &loop, &QEventLoop::quit, Qt::QueuedConnection);

        viewer.show();
        viewer.loadSubject(job->tmpl,
            QString::fromStdString(cfg.template_path),
            QString::fromStdString(cfg.fiducial_marker_locations),
            fileId, cfg.ecg_upsample_rate,
            cfg.ppg_upsample_rate, cfg.abp_upsample_rate,
            cfg.art_upsample_rate, cfg.art_pulm_upsample_rate,
            cfg.notch_filter_hz);
        loop.exec();

        outBanks.clear();
        outBanks.reserve(viewer.bins().size());
        for (const TemplateBin& b : viewer.bins()) {
            outBanks.push_back(analysis_job::BankSnapshot{ b.ecg_bank, b.ppg_bank });
        }
}


int main(int argc, char* argv[]) {
    //set up gui
    QApplication app(argc, argv);
    Theme::apply(app);

    //get data type and initials for logging
    const int dataset_choice = get_dataset_choice();
    config_entry cfg;
    if (!load_config(dataset_choice, cfg)) {
        std::cerr << "Error Loading config.csv\n";
        return 1;
    }
    const std::string initials = get_initials();
    cfg.log_path = cfg.output_path + "/log_" + initials;
    std::filesystem::create_directories(cfg.log_path);

	//load the bin files
    const std::vector<std::filesystem::path> binFiles = load_binfiles(cfg);
    if (binFiles.empty()) {
        std::cerr << "No .bin files in: " << cfg.bin_file_path << "\n";
        return 0;
    }

    for (const std::filesystem::path& binFs : binFiles) {
        const std::string stem = binFs.stem().string();
        std::cout << "\nProcessing file:" << stem << "\n";

        // Skip files that already have a log - if they don't have a log, make one
        const std::string logPath = cfg.log_path + "/" + stem + "_log.csv";
        if (std::filesystem::exists(logPath)) {
            std::cout << "  log already exists in " << cfg.log_path << "; skipping " << stem << "\n";
            continue;
        }
        beat_log beatLog;
        beatLog.setDefaultParams(cfg.blanking_period, cfg.threshold);

		// Run the noise marking GUI
        std::cout << "Noise marking: " << binFs.filename().string() << "\n";
        QVector<GenExcStruct> allMarkings;
        std::filesystem::path currentBinFile;
        bool ecg1Inverted = false, ecg2Inverted = false, ecg3Inverted = false;

        // this is normally a no-op if the user loaded a bin file then this switches to the loaded file
        const std::filesystem::path effBin = currentBinFile.empty() ? binFs : currentBinFile;
        const std::string effStem = effBin.stem().string();

        // write CSV
        beatLog.flushPending();
        beatLog.writeCsv(cfg.log_path + "/" + effStem + "_log.csv");

        //export noise markings and load them in the template
        if (allMarkings.isEmpty()) {
            exportMarkings(cfg, currentBinFile, nullptr);
        }
        else {
            for (const GenExcStruct& m : allMarkings) {
                exportMarkings(cfg, std::filesystem::path(m.filePath.toStdString()), &m);
            }
        }

        auto jobOpt = analysis_job::prepare(cfg, effBin, ecg1Inverted, ecg2Inverted, ecg3Inverted);
        if (!jobOpt) {
            std::cout << "Bin Not Loaded\n";
            continue;
        }

        auto job = std::make_shared<analysis_job::AnalysisJob>(std::move(*jobOpt));

        std::thread worker([job] {
            analysis_job::finalize(*job);
            });
        std::vector<analysis_job::BankSnapshot> operatorBanks;
        runTemplateMarking(cfg, job, QString::fromStdString(job->stem), operatorBanks);
        worker.join();
        analysis_job::commit(*job, operatorBanks);
    }

    std::cout << "\nAll files processed.\n";
    return 0;
}