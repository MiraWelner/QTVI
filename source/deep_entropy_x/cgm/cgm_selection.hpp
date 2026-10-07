#pragma once
//
// cgm_selection.hpp
//
// BITTIUM ONLY, once per patient. For the Bittium .bin being processed, when
// config.csv has a cgm_folder path, open a Qt file dialog in that folder and
// let the operator click that patient's CGM files (one or more). Each picked
// file is processed on its own and written as
//
//     <output_folder>/cgm_output/<patientID>_<CGM file name>.bin
//     e.g. 024-02-03-25-49_ECG201V3_glucose_9-22-2026.bin
//
// (format: cgm_bin_io.hpp). Picking the same CGM file again overwrites its .bin.
//
// THE PATIENT ID comes from the Bittium .bin, NOT from the CGM file: the CGM
// names (ECG201V3_glucose_9-22-2026) and the Bittium IDs (024-02-03-25-49)
// share nothing. It is the .bin stem up to its first underscore, as
// parseTemplateFileName reads Bittium names:
//
//     024-02-03-25-49.bin            ->  024-02-03-25-49
//     024-02-03-25-49_1000_005.bin   ->  024-02-03-25-49
//
// Every way out is non-fatal: a blank cgm_folder, another dataset, a missing
// folder, a cancelled dialog or an unreadable file skips CGM (or that one
// file) and the ECG processing of this patient continues.
//
// Only selectAndProcessCgmFiles touches Qt; the naming and per-file helpers
// above it are plain C++ and are what the unit tests call.
//

#include "config_file_handling/config.hpp"
#include "cgm/cgm_io.hpp"
#include "cgm/cgm_pipeline.hpp"

#include <QFileDialog>
#include <QString>
#include <QStringList>

#include <filesystem>
#include <iostream>
#include <string>
#include <system_error>

namespace cgm_selection {

    inline std::string output_name(const std::string& patientId, const std::filesystem::path& cgmFile) {
        // <patientID>_<CGM file stem>.bin
        return patientId + "_" + cgmFile.stem().string() + ".bin";
    }

    struct CgmOutcome {
        bool        written = false;
        int         nReadings = 0;
        int         nValidBins = 0;
        cgm_pipeline::CgmPipelineResult result;   // full result, for the summary line
        cgm_bin_io::CgmBinFile file;               // what was written
        std::filesystem::path output;
        std::string error;
    };

    /// Read one CGM export, process it, write <patientID>_<name>.bin.
    inline CgmOutcome processCgmFile(const std::string& patientId,
        const std::filesystem::path& cgmFile,
        const std::filesystem::path& outDir,
        const cgm_pipeline::CgmParams& params = cgm_pipeline::CgmParams{}) {
        CgmOutcome R;
        const auto readings = cgm_pipeline::readLibreCsv(cgmFile.string());
        R.nReadings = static_cast<int>(readings.size());
        if (readings.empty()) { R.error = "no historic glucose rows"; return R; }

        R.result = cgm_pipeline::runCgmPipeline(readings, params);
        for (const auto& b : R.result.bins) R.nValidBins += b.raw.valid ? 1 : 0;
        R.file = cgm_bin_io::makeCgmBinFile(R.result, params);

        R.output = outDir / output_name(patientId, cgmFile);
        if (!cgm_bin_io::writeCgmBin(R.output.string(), R.file)) {
            R.error = "cannot write " + R.output.string();
            return R;
        }
        R.written = true;
        return R;
    }

    /// binStem: stem of the Bittium .bin for this patient. Returns the number
    /// of CGM .bin files written.
    inline int selectAndProcessCgmFiles(const config_entry& cfg, const std::string& binStem) {
        if (cfg.dataset_type != "BITTIUM" || cfg.cgm_folder.empty()) return 0;

        std::error_code ec;
        if (!std::filesystem::is_directory(cfg.cgm_folder, ec)) {
            std::cerr << "WARNING: cgm_folder '" << cfg.cgm_folder
                << "' is not a folder; skipping CGM\n";
            return 0;
        }
        if (cfg.cgm_output_path.empty()) {
            std::cerr << "WARNING: no cgm_output folder; skipping CGM\n";
            return 0;
        }

        const QStringList picked = QFileDialog::getOpenFileNames(
            nullptr, QString::fromStdString("Select CGM files for " + binStem),
            QString::fromStdString(cfg.cgm_folder),
            "CGM exports (*.csv);;All files (*)");
        if (picked.isEmpty()) {
            std::cout << "  [cgm] no CGM files selected for " << binStem << "\n";
            return 0;
        }

        cgm_pipeline::CgmParams params;           // spec defaults, plus:
        params.binMinutes = cfg.cgm_bin_minutes;  // config.csv cgm_bin_minutes

        int written = 0;
        for (const QString& q : picked) {
            const std::filesystem::path in = q.toStdString();
            const auto R = processCgmFile(id, in, cfg.cgm_output_path, params);
            const auto& res = R.result;
            ++written;
        }
        return written;
    }

}  // namespace cgm_selection