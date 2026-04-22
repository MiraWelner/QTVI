/**
 * @file   file_to_bin.cpp
 * @brief  Take in a MESA, Bittium, or CHAOS file and converts it to a .bin file of uniform format
 *
 * The format is the following:
 *  * 180-byte header (45 x uint32):
 *   Offset  0: signal_rate           (uint32)  -Sampling rate for continuous signals (1000 Hz after upsampling)
 *   Offset  4: boolean_rate          (uint32)  -Sampling rate for booleans that were originally 1 Hz
 *   Offset  8: pacemaker_event_rate  (uint32)  -Pacemaker event epoch duration (8 Hz, not upsampled)
 *   Offset 12: sleep_state_rate      (uint32)  -Sleep stage epoch duration in seconds (30s for MESA)
 *   Offset 16: size_ecg_1            (uint32)
 *   Offset 20: size_ecg_2            (uint32)
 *   Offset 24: size_ecg_3            (uint32)
 *   Offset 28: size_ppg              (uint32)
 *   Offset 32: size_accel_x          (uint32)
 *   Offset 36: size_accel_y          (uint32)
 *   Offset 40: size_accel_z          (uint32)
 *   Offset 44: size_marker           (uint32)
 *   Offset 48: size_temp             (uint32)
 *   Offset 52: size_pacemaker        (uint32)
 *   Offset 56: size_eog_l            (uint32)
 *   Offset 60: size_eog_r            (uint32)
 *   Offset 64: size_emg              (uint32)
 *   Offset 68: size_eeg_1            (uint32)
 *   Offset 72: size_eeg_2            (uint32)
 *   Offset 76: size_eeg_3            (uint32)
 *   Offset 80: size_pres             (uint32)   -MESA nasal pressure / Bittium+CHAOS CVP
 *   Offset 84: size_flow             (uint32)
 *   Offset 88: size_thor             (uint32)
 *   Offset 92: size_abdo             (uint32)
 *   Offset 96: size_leg              (uint32)
 *   Offset 100: size_therm           (uint32)
 *   Offset 104: size_pos             (uint32)
 *   Offset 108: size_ekg_off         (uint32)
 *   Offset 112: size_eog_l_off       (uint32)
 *   Offset 116: size_eog_r_off       (uint32)
 *   Offset 120: size_emg_off         (uint32)
 *   Offset 124: size_eeg1_off        (uint32)
 *   Offset 128: size_eeg2_off        (uint32)
 *   Offset 132: size_eeg3_off        (uint32)
 *   Offset 136: size_oxstatus        (uint32)
 *   Offset 140: size_spo2            (uint32)
 *   Offset 144: size_HR              (uint32)
 *   Offset 148: size_DHR             (uint32)
 *   Offset 152: size_resp            (uint32)
 *   Offset 156: size_abp             (uint32)   -Bittium+CHAOS arterial blood pressure (NLS_NOM_PRESS_BLD_ART_ABP)
 *   Offset 160: size_eeg_4           (uint32)   -Bittium+CHAOS 4th EEG channel
 *   Offset 164: size_art             (uint32)   -CHAOS systemic arterial (NLS_NOM_PRESS_BLD_ART)
 *   Offset 168: size_art_pulm        (uint32)   -CHAOS pulmonary arterial (NLS_NOM_PRESS_BLD_ART_PULM)
 *   Offset 172: size_sleep           (uint32)
 *   (total header: 4 rate + 40 channel-sizes + 1 sleep-size = 45 x uint32 = 180 bytes)
 *
 * Signal data (contiguous doubles, immediately after header):
 *   Written in the same order as the size fields above.
 *   Channels at >1 Hz are upsampled to 1 kHz.
 *   Channels already at 1 Hz are written raw (no upsampling).
 *   Missing channels are stored as a single -1.0 with size = 1.
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-03-18
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <filesystem>
#include <set>
#include <cmath>

extern "C" {
#include "edflib.h"
}
#include "pugixml.hpp"
#include "resample.hpp"

// 45 uint32 fields = 180 bytes (4 rate + 40 chan-size + 1 sleep-size)
static const int NUM_HEADER_FIELDS = 45;
static const std::streamoff HEADER_SIZE = NUM_HEADER_FIELDS * sizeof(uint32_t);
static const double SLEEP_STATE_LENGTH = 30.0;
static const std::string CONFIG_PATH = "config.csv";
static const double final_sampling_rate = 1000.0;
static const double BOOLEAN_RATE = 1.0;  // 1 Hz channels -do not upsample

struct config_csv_data {
    std::string dataType, mainExt, sleepExt, inputPath, outputPath;
    std::string ecg1Label, ecg2Label, ecg3Label, ppgLabel;
    double ecgRate, ppgRate;
};

std::vector<std::string> parse_csv_row(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    for (size_t i = 0; i < line.length(); ++i) {
        if (line[i] == ',') { fields.push_back(cur); cur = ""; }
        else cur += line[i];
    }
    fields.push_back(cur);
    for (auto& f : fields) {
        size_t first = f.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) f = "";
        else {
            size_t last = f.find_last_not_of(" \t\r\n");
            f = f.substr(first, last - first + 1);
        }
    }
    return fields;
}

bool contains(std::string search_string, std::string substring) {
    if (substring.empty()) return false;
    std::transform(search_string.begin(), search_string.end(), search_string.begin(), ::toupper);
    std::transform(substring.begin(), substring.end(), substring.begin(), ::toupper);
    return search_string.find(substring) != std::string::npos;
}

// Case-insensitive equality check on trimmed strings.
// Use this instead of contains() when channel names share a common prefix
// (e.g. NLS_NOM_PRESS_BLD_ART vs NLS_NOM_PRESS_BLD_ART_PULM vs NLS_NOM_PRESS_BLD_ART_ABP).
static bool equals_ci(std::string a, std::string b) {
    std::transform(a.begin(), a.end(), a.begin(), ::toupper);
    std::transform(b.begin(), b.end(), b.begin(), ::toupper);
    return a == b;
}

// ============================================================================
// Locate the real CSV header row in a CHAOS/Bittium .dat file. Some files have
// a metadata line and a confidentiality NOTICE before the real CSV header, so
// "first non-empty line" isn't reliable. We treat a line as the header iff it
// contains either a known channel-name token (NLS_NOM_, NLS_EEG_) or the words
// "Index" and "TimeStamp" on the same line.
//
// Leaves `in` positioned at the line AFTER the header (first data row).
// Returns the parsed header cells, or an empty vector if no header was found.
// ============================================================================
static std::vector<std::string> find_real_header(std::istream& in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        // Cheap pre-check: the header row always contains at least one known token.
        bool looksLikeHeader =
            contains(line, "NLS_NOM_") ||
            contains(line, "NLS_EEG_") ||
            (contains(line, "Index") && contains(line, "TimeStamp"));
        if (!looksLikeHeader) continue;
        return parse_csv_row(line);
    }
    return {};
}


// ============================================================================
// Infer the file's row rate (Hz) from the difference between the first two
// timestamps. In CHAOS .dat files, every row has a monitor timestamp even if
// most data columns on that row are empty, so this gives a single authoritative
// sample-slot rate for the whole file.
//
// Returns 0.0 on failure.
// ============================================================================
static double infer_row_rate(const std::filesystem::path& path,
    const std::string& tsColumnName)
{
    std::ifstream in(path);
    if (!in) return 0.0;

    std::vector<std::string> hdrs = find_real_header(in);
    if (hdrs.empty()) return 0.0;

    int tsCol = -1;
    for (int i = 0; i < (int)hdrs.size(); ++i) {
        if (contains(hdrs[i], tsColumnName)) { tsCol = i; break; }
    }
    if (tsCol < 0) return 0.0;

    auto parseMs = [](const std::string& s) -> long long {
        size_t sp = s.find_last_of(' ');
        std::string t = (sp == std::string::npos) ? s : s.substr(sp + 1);
        int hh = 0, mm = 0, ss = 0, ms = 0;
        if (std::sscanf(t.c_str(), "%d:%d:%d.%d", &hh, &mm, &ss, &ms) < 3)
            return -1;
        return ((long long)hh * 3600 + mm * 60 + ss) * 1000LL + ms;
        };

    long long firstMs = -1, secondMs = -1, lastMs = -1;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> cells = parse_csv_row(line);
        if (tsCol >= (int)cells.size()) continue;
        const std::string& ts = cells[tsCol];
        if (ts.empty()) continue;
        long long ms = parseMs(ts);
        if (ms < 0) continue;

        // Track the last-seen timestamp for the duration sanity check below.
        lastMs = ms;

        if (firstMs < 0) { firstMs = ms; continue; }
        if (secondMs < 0) {
            if (ms == firstMs) continue;   // same-ms duplicate; wait for the stride to appear
            secondMs = ms;
            // don't break -- keep reading so lastMs reaches end of file
        }
    }
    if (firstMs < 0 || secondMs < 0 || secondMs <= firstMs) return 0.0;

    // Sanity check: warn if the file's actual duration is wildly off expected.
    if (lastMs > firstMs) {
        double secs = (lastMs - firstMs) / 1000.0;
        if (secs < 14 * 60 || secs > 16 * 60) {
            std::cout << "Unexpected file length: " << secs << " seconds\n";
        }
    }

    double dtSec = (secondMs - firstMs) / 1000.0;
    if (dtSec <= 0.0) return 0.0;
    return 1.0 / dtSec;
}


// ============================================================================
// Read one column from a CHAOS .dat file into a dense vector of length N
// (= total number of data rows).  Empty cells are filled by linear
// interpolation between the surrounding populated cells.  Leading/trailing
// empty cells are filled with the first/last populated value, respectively.
//
// Returns an empty vector if the column is either absent or entirely empty.
//
// - If `exactMatch` is true, header names must match `label` exactly
//   (case-insensitive). Use this when multiple columns share a prefix
//   (e.g. NLS_NOM_PRESS_BLD_ART vs NLS_NOM_PRESS_BLD_ART_ABP).
// - If false, `contains()` matching is used.
// ============================================================================
static std::vector<double> read_dat_column_interpolated(
    const std::filesystem::path& path,
    const std::string& label,
    bool exactMatch)
{
    std::ifstream in(path);
    if (!in || label.empty()) return {};

    std::vector<std::string> hdrs = find_real_header(in);
    if (hdrs.empty()) return {};

    int colIdx = -1;
    for (int i = 0; i < (int)hdrs.size(); ++i) {
        bool match = exactMatch ? equals_ci(hdrs[i], label)
            : contains(hdrs[i], label);
        if (match) { colIdx = i; break; }
    }
    if (colIdx < 0) return {};

    // First pass: read all rows. `present[i]` flags whether row i had a value.
    std::vector<double> values;
    std::vector<bool>   present;
    values.reserve(1 << 18);
    present.reserve(1 << 18);

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        if (colIdx >= (int)row.size() || row[colIdx].empty()) {
            values.push_back(0.0);
            present.push_back(false);
            continue;
        }
        try {
            values.push_back(std::stod(row[colIdx]));
            present.push_back(true);
        }
        catch (...) {
            values.push_back(0.0);
            present.push_back(false);
        }
    }

    const size_t n = values.size();
    if (n == 0) return {};

    // Is there any populated cell at all?
    size_t firstPresent = n;
    size_t lastPresent = n;
    for (size_t i = 0; i < n; ++i)
        if (present[i]) { firstPresent = i; break; }
    if (firstPresent == n) return {};   // entirely empty column
    for (size_t i = n; i-- > 0; )
        if (present[i]) { lastPresent = i; break; }

    // Leading run of missings: fill with first known value.
    for (size_t i = 0; i < firstPresent; ++i) {
        values[i] = values[firstPresent];
        present[i] = true;
    }
    // Trailing run of missings: fill with last known value.
    for (size_t i = lastPresent + 1; i < n; ++i) {
        values[i] = values[lastPresent];
        present[i] = true;
    }

    // Interior gaps: walk forward. Whenever we enter a run of missings,
    // find the next present index and linearly interpolate across.
    size_t i = firstPresent;
    while (i < lastPresent) {
        if (present[i + 1]) { ++i; continue; }
        // Find the next present index j.
        size_t j = i + 2;
        while (j <= lastPresent && !present[j]) ++j;
        // i is present, j is present, (i, j) are missings.
        const double v0 = values[i];
        const double v1 = values[j];
        const double span = static_cast<double>(j - i);
        for (size_t k = i + 1; k < j; ++k) {
            const double f = static_cast<double>(k - i) / span;
            values[k] = v0 * (1.0 - f) + v1 * f;
            present[k] = true;
        }
        i = j;
    }

    return values;
}


// ============================================================================
// Write an EDF channel to binary output.
// If skip_resample is true, writes raw samples without upsampling.
// If channel index is invalid, writes a single -1.0 placeholder.
// ============================================================================
void edf_to_bin(int handle, int idx, long long n, double old_rate,
    std::ofstream& out, uint32_t& sizeOut, bool skip_resample = false) {
    if (idx < 0 || n <= 0) {
        double dummy = -1.0;
        out.write((char*)&dummy, 8);
        sizeOut = 1;
        return;
    }

    std::vector<double> buf(n);
    edfread_physical_samples(handle, idx, (int)n, buf.data());

    if (!skip_resample) {
        buf = upsample(buf, old_rate);
    }
    out.write((char*)buf.data(), buf.size() * 8);
    sizeOut = (uint32_t)buf.size();
}


// ============================================================================
// Write a placeholder (missing channel): single -1.0
// ============================================================================
static void write_missing(std::ofstream& out, uint32_t& sizeOut) {
    double v = -1.0;
    out.write((char*)&v, 8);
    sizeOut = 1;
}

// ============================================================================
// Config loader
// ============================================================================
static bool load_config(int data_type, config_csv_data& out) {
    std::ifstream file(CONFIG_PATH);
    if (!file.is_open()) return false;

    std::string target = (data_type == 1) ? "MESA" : (data_type == 2) ? "BITTIUM" : (data_type == 3) ? "CHAOS" : "";
    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        if (row.size() < 12) continue;

        std::string rType = row[0];
        std::transform(rType.begin(), rType.end(), rType.begin(), ::toupper);

        if (rType == target) {
            out.dataType = row[0];
            out.mainExt = row[1];
            out.sleepExt = row[2];
            out.inputPath = row[3];
            out.outputPath = row[4];
            out.ecg1Label = row[6];
            out.ecg2Label = row[7];
            out.ecg3Label = row[8];
            out.ppgLabel = row[9];
            try {
                out.ecgRate = (!row[10].empty()) ? std::stod(row[10]) : 0.0;
                out.ppgRate = (!row[11].empty()) ? std::stod(row[11]) : 0.0;
            }
            catch (...) {
                out.ecgRate = 256.0;
                out.ppgRate = 256.0;
            }
            return true;
        }
    }
    return false;
}

// ============================================================================
// Channel map -holds EDF signal index for each channel
// ============================================================================
struct ChannelMap {
    int ecg1 = -1, ecg2 = -1, ecg3 = -1, ppg = -1;
    int accel_x = -1, accel_y = -1, accel_z = -1;
    int marker = -1, temp = -1, pacemaker = -1;
    int eog_l = -1, eog_r = -1, emg = -1;
    int eeg1 = -1, eeg2 = -1, eeg3 = -1;
    int pres = -1, flow = -1, thor = -1, abdo = -1, leg = -1, therm = -1;
    int pos = -1;
    int ekg_off = -1, eog_l_off = -1, eog_r_off = -1, emg_off = -1;
    int eeg1_off = -1, eeg2_off = -1, eeg3_off = -1;
    int oxstatus = -1, spo2 = -1, hr = -1, dhr = -1, resp = -1;
    int abp = -1, eeg4 = -1;
    int art = -1, art_pulm = -1;
};

static ChannelMap build_edf_channel_map(const edf_hdr_struct* hdr, const config_csv_data& cfg) {
    ChannelMap cm;
    std::set<int> used;

    auto find = [&](const std::string& label) -> int {
        if (label.empty()) return -1;
        for (int i = 0; i < hdr->edfsignals; ++i) {
            if (!used.count(i) && contains(hdr->signalparam[i].label, label)) {
                used.insert(i);
                return i;
            }
        }
        return -1;
        };

    // Primary channels from config
    cm.ecg1 = find(cfg.ecg1Label);
    cm.ecg2 = find(cfg.ecg2Label);
    cm.ecg3 = find(cfg.ecg3Label);
    cm.ppg = find(cfg.ppgLabel);

    // Accelerometers (CHAOS)
    cm.accel_x = find("Accelerometer_X");
    cm.accel_y = find("Accelerometer_Y");
    cm.accel_z = find("Accelerometer_Z");

    // Boolean/low-rate channels (CHAOS)
    cm.marker = find("Marker");
    cm.temp = find("DEV_Temperature");
    cm.pacemaker = -1;

    // MESA polysomnography channels
    cm.eog_l = find("EOG-L");
    cm.eog_r = find("EOG-R");
    cm.emg = find("EMG");
    cm.eeg1 = find("EEG1");
    cm.eeg2 = find("EEG2");
    cm.eeg3 = find("EEG3");

    // MESA respiratory/other channels
    cm.pres = find("Pres");
    cm.flow = find("Flow");
    cm.thor = find("Thor");
    cm.abdo = find("Abdo");
    cm.leg = find("Leg");
    cm.therm = find("Therm");
    cm.pos = find("Pos");

    // MESA offset channels (1 Hz)
    cm.ekg_off = find("EKG_Off");
    cm.eog_l_off = find("EOG-L_Off");
    cm.eog_r_off = find("EOG-R_Off");
    cm.emg_off = find("EMG_Off");
    cm.eeg1_off = find("EEG1_Off");
    cm.eeg2_off = find("EEG2_Off");
    cm.eeg3_off = find("EEG3_Off");

    // MESA oximetry/heart rate
    cm.oxstatus = find("OxStatus");
    cm.spo2 = find("SpO2");
    if (cm.spo2 < 0) cm.spo2 = find("Sp02");
    cm.hr = find("HR");
    cm.dhr = find("DHR");

    // Respiratory (Bittium)
    cm.resp = find("NLS_NOM_RESP");
    if (cm.resp < 0) cm.resp = find("Resp");

    return cm;
}

// ============================================================================
// EDF helpers
// ============================================================================
static double edf_channel_rate(const edf_hdr_struct* hdr, int idx) {
    if (idx < 0) return 0.0;
    return (double)hdr->signalparam[idx].smp_in_datarecord /
        ((double)hdr->datarecord_duration / 10000000.0);
}

static long long edf_samples(const edf_hdr_struct* hdr, int idx) {
    if (idx < 0) return 0;
    return hdr->signalparam[idx].smp_in_file;
}

// ============================================================================
// Make binary file from EDF
// ============================================================================
static void make_binfile_edf(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_csv_data& cfg) {
    std::filesystem::path outPath = std::filesystem::path(cfg.outputPath) / (path.stem().string() + ".bin");
    std::ofstream out(outPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        return;
    }

    std::cout << "Processing EDF: " << path.filename().string() << std::endl;

    // Placeholder header
    std::vector<char> zeroes(HEADER_SIZE, 0);
    out.write(zeroes.data(), HEADER_SIZE);

    auto hdr = std::make_unique<edf_hdr_struct>();
    if (edfopen_file_readonly(path.string().c_str(), hdr.get(), EDFLIB_READ_ALL_ANNOTATIONS)) {
        std::cerr << "ERROR: cannot open EDF " << path << "\n";
        out.close();
        std::filesystem::remove(outPath);
        return;
    }

    ChannelMap cm = build_edf_channel_map(hdr.get(), cfg);

    uint32_t sizes[40] = {};

    // Helper: write one EDF channel.
    // Determines skip_resample automatically: if the channel's native rate
    // is <= BOOLEAN_RATE (1 Hz), write raw without upsampling.
    auto writeChannel = [&](int chIdx, int sizeIdx, double rateOverride = 0.0) {
        double rate = (rateOverride > 0.0) ? rateOverride : edf_channel_rate(hdr.get(), chIdx);
        bool skip = (chIdx >= 0) && (rate <= BOOLEAN_RATE);
        edf_to_bin(hdr->handle, chIdx, edf_samples(hdr.get(), chIdx), rate, out, sizes[sizeIdx], skip);
        };

    //  0-3: ECG 1-3, PPG  (upsample from ecgRate/ppgRate)
    writeChannel(cm.ecg1, 0, cfg.ecgRate);
    writeChannel(cm.ecg2, 1, cfg.ecgRate);
    writeChannel(cm.ecg3, 2, cfg.ecgRate);
    writeChannel(cm.ppg, 3, cfg.ppgRate);

    //  4-6: Accelerometers (upsample from native rate, e.g. 25 Hz in CHAOS)
    writeChannel(cm.accel_x, 4);
    writeChannel(cm.accel_y, 5);
    writeChannel(cm.accel_z, 6);

    //  7-9: Marker (1 Hz, no upsample), Temperature (1 Hz, no upsample), Pacemaker
    writeChannel(cm.marker, 7);
    writeChannel(cm.temp, 8);
    writeChannel(cm.pacemaker, 9);

    // 10-12: EOG-L, EOG-R, EMG (256 Hz in MESA -upsample)
    writeChannel(cm.eog_l, 10);
    writeChannel(cm.eog_r, 11);
    writeChannel(cm.emg, 12);

    // 13-15: EEG 1-3 (256 Hz in MESA -upsample)
    writeChannel(cm.eeg1, 13);
    writeChannel(cm.eeg2, 14);
    writeChannel(cm.eeg3, 15);

    // 16-22: Pres, Flow, Thor, Abdo, Leg, Therm, Pos (32 Hz in MESA -upsample)
    writeChannel(cm.pres, 16);
    writeChannel(cm.flow, 17);
    writeChannel(cm.thor, 18);
    writeChannel(cm.abdo, 19);
    writeChannel(cm.leg, 20);
    writeChannel(cm.therm, 21);
    writeChannel(cm.pos, 22);

    // 23-29: Offset channels (1 Hz -no upsample)
    writeChannel(cm.ekg_off, 23);
    writeChannel(cm.eog_l_off, 24);
    writeChannel(cm.eog_r_off, 25);
    writeChannel(cm.emg_off, 26);
    writeChannel(cm.eeg1_off, 27);
    writeChannel(cm.eeg2_off, 28);
    writeChannel(cm.eeg3_off, 29);

    // 30-31: OxStatus (1 Hz -no upsample), SpO2 (1 Hz -no upsample)
    writeChannel(cm.oxstatus, 30);
    writeChannel(cm.spo2, 31);

    // 32: HR (1 Hz -no upsample)
    writeChannel(cm.hr, 32);

    // 33: DHR (256 Hz in MESA -upsample)
    writeChannel(cm.dhr, 33);

    // 34: Resp
    writeChannel(cm.resp, 34);

    // 35: ABP (Bittium only -not in MESA EDF)
    writeChannel(cm.abp, 35);

    // 36: EEG4 (Bittium only -not in MESA EDF)
    writeChannel(cm.eeg4, 36);

    // 37: ART systemic (CHAOS only -not in MESA EDF)
    writeChannel(cm.art, 37);

    // 38: ART_PULM pulmonary (CHAOS only -not in MESA EDF)
    writeChannel(cm.art_pulm, 38);

    edfclose_file(hdr->handle);

    // Sleep stages from XML
    std::vector<double> stages;
    if (!cfg.sleepExt.empty() && !xmlPath.empty() && std::filesystem::exists(xmlPath)) {
        pugi::xml_document doc;
        if (doc.load_file(xmlPath.string().c_str())) {
            for (auto node : doc.select_nodes("//SleepStage")) {
                double v = node.node().text().as_double();
                stages.push_back(v == 5.0 ? 4.0 : v);
            }
        }
    }
    if (stages.empty()) {
        stages.push_back(-1.0);
    }
    uint32_t ss = (uint32_t)stages.size();
    out.write(reinterpret_cast<const char*>(stages.data()), ss * sizeof(double));

    // Write header
    out.seekp(0);

    uint32_t sig_rate = (uint32_t)final_sampling_rate;
    uint32_t bool_rate = (uint32_t)BOOLEAN_RATE;
    uint32_t pace_rate = 8;
    uint32_t sleep_rate = (uint32_t)SLEEP_STATE_LENGTH;

    out.write((char*)&sig_rate, 4);
    out.write((char*)&bool_rate, 4);
    out.write((char*)&pace_rate, 4);
    out.write((char*)&sleep_rate, 4);

    for (int i = 0; i < 40; ++i)
        out.write((char*)&sizes[i], 4);

    out.write((char*)&ss, 4);

    out.close();
}

// ============================================================================
// Make binary file from .dat (Bittium CSV)
// ============================================================================
static void make_binfile_dat(const std::filesystem::path& path,
    const config_csv_data& cfg) {
    std::filesystem::path outPath = std::filesystem::path(cfg.outputPath) / (path.stem().string() + ".bin");
    std::ofstream out(outPath, std::ios::binary | std::ios::trunc);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot create " << outPath << "\n";
        return;
    }

    std::cout << "Processing DAT: " << path.filename().string() << std::endl;

    std::vector<char> zeroes(HEADER_SIZE, 0);
    out.write(zeroes.data(), HEADER_SIZE);

    uint32_t sizes[40] = {};

    // --- Row rate from the first two monitor timestamps ---
    // Every row in a CHAOS .dat has a Monitor TimeStamp even when most data
    // cells on that row are blank, so the row rate is one authoritative
    // "sample slot rate" that applies to every column in the file.
    double row_rate = infer_row_rate(path, "Monitor TimeStamp");
    if (row_rate <= 0.0) row_rate = infer_row_rate(path, "System TimeStamp UTC");
    if (row_rate <= 0.0) {
        // Last-resort fallback: guess from the config's ECG rate.
        row_rate = (cfg.ecgRate > 0.0) ? cfg.ecgRate : 500.0;
        std::cout << "  [warn] couldn't infer row rate; using fallback "
            << row_rate << " Hz\n";
    }
    else {
        std::cout << "  [info] row rate = " << row_rate << " Hz\n";
    }

    // Helper: read a column, interpolate across gaps, upsample from row_rate
    // to the target rate, and write to `out`. Writes a single -1.0 placeholder
    // (sizeOut=1) if the column is absent or completely empty.
    auto writeCol = [&](const std::string& label, bool exact,
        uint32_t& sizeOut) {
            std::vector<double> samples =
                read_dat_column_interpolated(path, label, exact);
            if (samples.empty()) {
                double v = -1.0;
                out.write((char*)&v, 8);
                sizeOut = 1;
                return;
            }
            samples = upsample(samples, row_rate);
            out.write((char*)samples.data(), samples.size() * 8);
            sizeOut = (uint32_t)samples.size();
        };

    // 0-3: ECG 1-3, PPG
    writeCol(cfg.ecg1Label, false, sizes[0]);
    writeCol(cfg.ecg2Label, false, sizes[1]);
    writeCol(cfg.ecg3Label, false, sizes[2]);
    writeCol(cfg.ppgLabel, false, sizes[3]);

    // 4-6: Accelerometers (not present in CHAOS .dat)
    write_missing(out, sizes[4]);
    write_missing(out, sizes[5]);
    write_missing(out, sizes[6]);

    // 7-9: Marker, Temp, Pacemaker (not present)
    write_missing(out, sizes[7]);
    write_missing(out, sizes[8]);
    write_missing(out, sizes[9]);

    // 10-12: EOG, EMG (not present)
    write_missing(out, sizes[10]);
    write_missing(out, sizes[11]);
    write_missing(out, sizes[12]);

    // 13-15: EEG 1-3 (contains-match; the real columns have _LBL suffix)
    writeCol("NLS_EEG_NAMES_EEG_CHAN1", false, sizes[13]);
    writeCol("NLS_EEG_NAMES_EEG_CHAN2", false, sizes[14]);
    writeCol("NLS_EEG_NAMES_EEG_CHAN3", false, sizes[15]);

    // 16: Pres / CVP
    writeCol("NLS_NOM_PRESS_BLD_VEN_CENT", false, sizes[16]);

    // 17-22: Flow, Thor, Abdo, Leg, Therm, Pos (not present)
    for (int i = 17; i <= 22; ++i) write_missing(out, sizes[i]);

    // 23-29: MESA offset channels (not present)
    for (int i = 23; i <= 29; ++i) write_missing(out, sizes[i]);

    // 30-33: OxStatus, SpO2, HR, DHR (not present)
    for (int i = 30; i <= 33; ++i) write_missing(out, sizes[i]);

    // 34: Resp
    writeCol("NLS_NOM_RESP", false, sizes[34]);
    // 35: ABP (exact match -- distinguishes from ART / ART_PULM)
    writeCol("NLS_NOM_PRESS_BLD_ART_ABP", true, sizes[35]);
    // 36: EEG4 (contains-match)
    writeCol("NLS_EEG_NAMES_EEG_CHAN4", false, sizes[36]);
    // 37: ART systemic arterial (exact match)
    writeCol("NLS_NOM_PRESS_BLD_ART", true, sizes[37]);
    // 38: ART_PULM pulmonary arterial (exact match)
    writeCol("NLS_NOM_PRESS_BLD_ART_PULM", true, sizes[38]);

    // No sleep data for Bittium or CHAOS CSV
    std::vector<double> stages = { -1.0 };
    uint32_t ss = 1;
    out.write(reinterpret_cast<const char*>(stages.data()), ss * sizeof(double));

    // Write header
    out.seekp(0);

    uint32_t sig_rate = (uint32_t)final_sampling_rate;
    uint32_t bool_rate = (uint32_t)BOOLEAN_RATE;
    uint32_t pace_rate = 8;
    uint32_t sleep_rate = (uint32_t)SLEEP_STATE_LENGTH;

    out.write((char*)&sig_rate, 4);
    out.write((char*)&bool_rate, 4);
    out.write((char*)&pace_rate, 4);
    out.write((char*)&sleep_rate, 4);

    for (int i = 0; i < 40; ++i)
        out.write((char*)&sizes[i], 4);

    out.write((char*)&ss, 4);

    out.close();
}

// ============================================================================
// Dispatcher
// ============================================================================
static void make_binfile(const std::filesystem::path& path,
    const std::filesystem::path& xmlPath,
    const config_csv_data& cfg) {
    std::string ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

    if (ext == ".EDF") {
        make_binfile_edf(path, xmlPath, cfg);
    }
    else if (ext == ".DAT" || ext == ".CSV") {
        make_binfile_dat(path, cfg);
    }
    else {
        std::cerr << "ERROR: unsupported file type " << ext << " for " << path << "\n";
    }
}

int main(int argc, char* argv[]) {
    std::cout << "FILE TO BIN\n";
    std::cout << "Select Dataset:\n1: MESA\n2: Bittium\n3: CHAOS\nChoice: ";
    int choice;
    std::cin >> choice;

    config_csv_data cfg;
    if (!load_config(choice, cfg)) {
        std::cerr << "Error: Could not find configuration for selection " << choice << " in config.csv" << std::endl;
        return 1;
    }

    std::filesystem::create_directories(cfg.outputPath);
    std::string tExt = cfg.mainExt;
    std::transform(tExt.begin(), tExt.end(), tExt.begin(), ::toupper);

    for (const auto& entry : std::filesystem::recursive_directory_iterator(cfg.inputPath)) {
        if (!entry.is_regular_file()) continue;

        std::string fExt = entry.path().extension().string();
        std::transform(fExt.begin(), fExt.end(), fExt.begin(), ::toupper);

        if (fExt == tExt) {
            std::filesystem::path xml;
            if (!cfg.sleepExt.empty()) {
                std::string stem = entry.path().stem().string();
                for (const auto& f : std::filesystem::directory_iterator(entry.path().parent_path())) {
                    std::string cExt = f.path().extension().string();
                    std::transform(cExt.begin(), cExt.end(), cExt.begin(), ::toupper);
                    if (cExt == cfg.sleepExt && f.path().stem().string().find(stem) != std::string::npos) {
                        xml = f.path();
                        break;
                    }
                }
            }
            make_binfile(entry.path(), xml, cfg);
        }
    }
    std::cout << "Processing Complete." << std::endl;
    return 0;
}