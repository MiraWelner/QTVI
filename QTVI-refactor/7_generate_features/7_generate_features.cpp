// All logic lives in header-only .hpp files.
// This is the sole translation unit.
#include "ppg_features.hpp"
#include "ppg_utils.hpp"
#include "find_foot_pulseox.hpp"
#include "dicrotic_sqi.hpp"
#include "beat_features.hpp"
#include "bin_io.hpp"
#include "flatten_and_generate.hpp"

#include <chrono>
#include <filesystem>

int main() {
    namespace fs = std::filesystem;
    using Clock = std::chrono::steady_clock;


    // ── Load config.csv ─────────────────────────────────────────────────
    auto configs = ppg::readConfigCsv("config.csv");
    if (configs.empty()) {
        std::cerr << "Error: Could not read config.csv or no entries found.\n";
        return 1;
    }

    // ── Dataset selection ───────────────────────────────────────────────
    std::cout << "\n=== Select Dataset ===\n";
    for (size_t i = 0; i < configs.size(); i++)
        std::cout << "  " << (i + 1) << ") " << configs[i].dataType << "\n";
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

    std::cout << "\nDataset:       " << cfg.dataType
        << "\nAnnealed path: " << cfg.annealedPath
        << "\nWave path:     " << cfg.wavePath
        << "\nTemplate path: " << cfg.templatePath
        << "\nMarking path:  " << cfg.markingPath
        << "\nOutput path:   " << cfg.featureOutputPath << "\n\n";

    bool use_templates = fs::exists(cfg.templatePath) &&
        !fs::is_empty(cfg.templatePath);
    bool skip_existing = true;
    std::cout << "Templates: " << (use_templates ? "enabled" : "disabled") << "\n\n";

    // ── Ensure output directory exists ───────────────────────────────────
    if (!fs::exists(cfg.featureOutputPath))
        fs::create_directories(cfg.featureOutputPath);

    // ── Discover matching files ─────────────────────────────────────────
    auto analysisFiles = ppg::find_analysis_files(
        cfg.annealedPath,
        cfg.wavePath,
        cfg.templatePath,
        cfg.markingPath,
        use_templates);

    std::cout << "Found " << analysisFiles.size() << " subjects to process.\n";
    std::cout << "*********************************************************************\n";

    double total_time = 0;
    int success = 0, fail = 0;

    for (size_t i = 0; i < analysisFiles.size(); i++) {
        const auto& af = analysisFiles[i];
        std::string out_file = (fs::path(cfg.featureOutputPath)
            / (af.uuid + "_feature_output.bin")).string();

        if (skip_existing && fs::exists(out_file)) {
            std::cout << af.uuid << "_feature_output.bin exists, skipping.\n";
            continue;
        }

        auto t_start = Clock::now();
        std::cout << "Beginning analysis of " << af.uuid
            << " | " << (i + 1) << " of " << analysisFiles.size() << "\n";

        double avg = (i > 0) ? total_time / i : 0;
        std::cout << "Avg Time (s): " << avg
            << "  Est finish (min): "
            << (avg * (analysisFiles.size() - i)) / 60.0 << "\n";

        int r = ppg::GenerateFeatures(
            af.anneal_path,
            af.wave_path,
            af.template_path,
            af.marking_path,
            cfg.featureOutputPath);

        if (r) success++; else fail++;

        double elapsed = std::chrono::duration<double>(Clock::now() - t_start).count();
        total_time += elapsed;
        std::cout << "Elapsed: " << elapsed << " s\n"
            << "____________________________________________________________________\n\n";
    }

    std::cout << "Done. Success: " << success << ", Fail: " << fail << "\n";
    return 0;
}