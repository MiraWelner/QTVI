/**
 * @file   file_to_bin.cpp
 * @brief  Take in a MESA, Bittium, or CHAOS file and convert it to a .bin
 *         file of uniform format.
 *
 * ===== Bin format (v2: 41 channels, per-channel native rates) =====
 *
 *   512-byte header (128 x 32-bit fields):
 *
 *     Offset   0: signal_rate           (uint32)  common upsampled rate (1000 Hz)
 *     Offset   4: boolean_rate          (uint32)  shared rate for <=1 Hz channels (1 Hz)
 *     Offset   8: pacemaker_event_rate  (uint32)  pacemaker epoch rate (8 Hz)
 *     Offset  12: sleep_state_rate      (uint32)  sleep-stage epoch length in seconds (30 s)
 *
 *     Offset  16: upsampled sizes       41 x uint32  (size_<chan>)
 *     Offset 180: raw-pair sizes        41 x uint32  (size_<chan>_raw, counts PAIRS)
 *     Offset 344: native sampling rates 41 x float32 (Hz)
 *                 0.0 = channel absent; negative values never used.
 *     Offset 508: size_sleep            (uint32)
 *
 *   Channel index order (41 slots, identical across upsampled/raw/native-rate blocks):
 *      0: timestamp      (seconds from start of recording)
 *      1: ecg_1
 *      2: ecg_2
 *      3: ecg_3
 *      4: ppg
 *      5: accel_x
 *      6: accel_y
 *      7: accel_z
 *      8: marker
 *      9: temp
 *     10: pacemaker
 *     11: eog_l
 *     12: eog_r
 *     13: emg
 *     14: eeg_1
 *     15: eeg_2
 *     16: eeg_3
 *     17: pres               (MESA nasal pressure / CHAOS CVP)
 *     18: flow
 *     19: thor
 *     20: abdo
 *     21: leg
 *     22: therm
 *     23: pos
 *     24: ekg_off            (1 Hz offset channels below)
 *     25: eog_l_off
 *     26: eog_r_off
 *     27: emg_off
 *     28: eeg1_off
 *     29: eeg2_off
 *     30: eeg3_off
 *     31: oxstatus
 *     32: spo2
 *     33: HR
 *     34: DHR
 *     35: resp
 *     36: abp                (Bittium+CHAOS arterial, NLS_NOM_PRESS_BLD_ART_ABP)
 *     37: eeg_4              (Bittium+CHAOS 4th EEG channel)
 *     38: art                (CHAOS systemic arterial, NLS_NOM_PRESS_BLD_ART)
 *     39: art_pulm           (CHAOS pulmonary arterial, NLS_NOM_PRESS_BLD_ART_PULM)
 *     40: (reserved / future use -- always written as absent)
 *
 *   Header size check: 4 + 41 + 41 + 41 + 1 = 128 fields x 4 bytes = 512 bytes.
 *
 * ===== Signal data layout =====
 *
 *   Immediately after the header, per channel in index order:
 *     1. Upsampled samples    size_<chan> doubles (at signal_rate, 1 kHz)
 *     2. Raw (timestamp, value) pairs  size_<chan>_raw * 2 doubles, interleaved
 *        as t0, v0, t1, v1, ... where t_k is in seconds from the start of
 *        recording.
 *
 *   Channels at <= boolean_rate (1 Hz) skip upsampling: the "upsampled" block
 *   is a copy of the native samples.
 *   Missing channels write a single -1.0 in the upsampled slot (size = 1) and
 *   a single (-1.0, -1.0) sentinel pair in the raw slot (size = 1 pair).
 *   Missing channels also get native_rate = 0.0.
 *   Sleep stages follow the signal blocks, as before.
 *
 * ===== Timestamp channel (index 0) =====
 *
 *   For CHAOS .dat files, the timestamp channel carries the real monitor
 *   wall-clock converted to seconds-from-start. Upsampled block is at
 *   1 kHz. Raw block stores (row_time, row_time_seconds) pairs -- trivially
 *   monotonic, but included for symmetry with other raw blocks.
 *   For EDF files, the timestamp channel is synthetic: the upsampled block
 *   is k / 1000 for k = 0..N-1 at 1 kHz, and the raw block stores
 *   (k/ecg_rate, k/ecg_rate) pairs at the primary ECG rate. Native rate
 *   stored is the primary ECG rate.
 *
 * @author Mira Welner
 * @email MEW386@pitt.edu
 * @date   2026-04-24
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

// 128 uint32-sized fields = 512 bytes
//   4 scalar rate fields
// + 41 upsampled-block sizes (uint32)
// + 41 raw-block sizes        (uint32)
// + 41 native sampling rates  (float32)
// + 1 sleep-size
static const int NUM_CHANNELS = 41;
static const int NUM_HEADER_FIELDS = 4 + 3 * NUM_CHANNELS + 1;   // = 128
static const std::streamoff HEADER_SIZE = NUM_HEADER_FIELDS * 4;  // = 512 bytes
static const double SLEEP_STATE_LENGTH = 30.0;
static const std::string CONFIG_PATH = "config.csv";
static const double final_sampling_rate = 1000.0;
static const double BOOLEAN_RATE = 1.0;  // 1 Hz channels -do not upsample

// Channel indices (one source of truth, used throughout).
enum ChannelIdx {
    CH_TIMESTAMP = 0,
    CH_ECG1, CH_ECG2, CH_ECG3, CH_PPG,
    CH_ACCEL_X, CH_ACCEL_Y, CH_ACCEL_Z,
    CH_MARKER, CH_TEMP, CH_PACEMAKER,
    CH_EOG_L, CH_EOG_R, CH_EMG,
    CH_EEG1, CH_EEG2, CH_EEG3,
    CH_PRES, CH_FLOW, CH_THOR, CH_ABDO,
    CH_LEG, CH_THERM, CH_POS,
    CH_EKG_OFF, CH_EOG_L_OFF, CH_EOG_R_OFF, CH_EMG_OFF,
    CH_EEG1_OFF, CH_EEG2_OFF, CH_EEG3_OFF,
    CH_OXSTATUS, CH_SPO2, CH_HR, CH_DHR,
    CH_RESP, CH_ABP, CH_EEG4,
    CH_ART, CH_ART_PULM,
    CH_RESERVED_40
};

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
// ============================================================================
// One-pass column store for a CHAOS .dat file. Reading the same file 80+ times
// to extract 41 channels is the main bottleneck in make_binfile_dat. This
// struct is populated by `prescan_dat_columns()` -- it walks the file once,
// parses each line once, and caches every column's populated values + their
// row indices keyed by header-column index. Downstream readers do O(1) lookups.
// ============================================================================
struct PrescannedDat {
    std::vector<std::string> headers;     // header cell names (case kept)
    size_t totalRows = 0;                 // number of data rows in the file
    // For each header column index: the row indices where the cell was
    // populated, and the parsed double value at that row.
    std::vector<std::vector<size_t>> rowIdxPerCol;
    std::vector<std::vector<double>> valuePerCol;
};

// Open `path` and parse it ONCE into a PrescannedDat.
// All numeric parsing happens here; subsequent column lookups are zero-IO.
static PrescannedDat prescan_dat_columns(const std::filesystem::path& path) {
    PrescannedDat out;
    std::ifstream in(path);
    if (!in) return out;

    out.headers = find_real_header(in);
    if (out.headers.empty()) return out;

    const size_t nCols = out.headers.size();
    out.rowIdxPerCol.assign(nCols, {});
    out.valuePerCol.assign(nCols, {});
    for (size_t c = 0; c < nCols; ++c) {
        out.rowIdxPerCol[c].reserve(1 << 16);
        out.valuePerCol[c].reserve(1 << 16);
    }

    std::string line;
    size_t rowIdx = 0;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        const size_t lim = std::min(row.size(), nCols);
        for (size_t c = 0; c < lim; ++c) {
            if (row[c].empty()) continue;
            try {
                double v = std::stod(row[c]);
                out.rowIdxPerCol[c].push_back(rowIdx);
                out.valuePerCol[c].push_back(v);
            }
            catch (...) {
                // ignore malformed cell
            }
        }
        ++rowIdx;
        if ((rowIdx % 1000000) == 0) {
            std::cout << "    [prescan] " << rowIdx / 1000000 << "M rows...\n" << std::flush;
        }
    }
    std::cout << "    [prescan] final row count: " << rowIdx << "\n";
    out.totalRows = rowIdx;
    return out;
}

// Match a label against the pre-parsed header. Returns column index or -1.
static int find_column_index(const PrescannedDat& dat,
    const std::string& label, bool exactMatch)
{
    if (label.empty()) return -1;
    for (int i = 0; i < (int)dat.headers.size(); ++i) {
        bool match = exactMatch ? equals_ci(dat.headers[i], label)
            : contains(dat.headers[i], label);
        if (match) return i;
    }
    return -1;
}


// ============================================================================
// Reconstruct an "interpolated" full-length column from prescanned data.
// Behaves identically to the old read_dat_column_interpolated but does NO
// file IO -- everything comes from the in-memory PrescannedDat.
// ============================================================================
static std::vector<double> column_interpolated(const PrescannedDat& dat,
    const std::string& label, bool exactMatch)
{
    int colIdx = find_column_index(dat, label, exactMatch);
    if (colIdx < 0) return {};

    const auto& rowIdx = dat.rowIdxPerCol[colIdx];
    const auto& vals = dat.valuePerCol[colIdx];
    if (rowIdx.empty()) return {};

    const size_t n = dat.totalRows;
    std::vector<double> values(n, 0.0);
    std::vector<bool>   present(n, false);
    for (size_t k = 0; k < rowIdx.size(); ++k) {
        if (rowIdx[k] < n) {
            values[rowIdx[k]] = vals[k];
            present[rowIdx[k]] = true;
        }
    }

    size_t firstPresent = n, lastPresent = n;
    for (size_t i = 0; i < n; ++i) if (present[i]) { firstPresent = i; break; }
    if (firstPresent == n) return {};
    for (size_t i = n; i-- > 0; ) if (present[i]) { lastPresent = i; break; }

    for (size_t i = 0; i < firstPresent; ++i) {
        values[i] = values[firstPresent]; present[i] = true;
    }
    for (size_t i = lastPresent + 1; i < n; ++i) {
        values[i] = values[lastPresent]; present[i] = true;
    }

    size_t i = firstPresent;
    while (i < lastPresent) {
        if (present[i + 1]) { ++i; continue; }
        size_t j = i + 2;
        while (j <= lastPresent && !present[j]) ++j;
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

// Raw populated cells + their row indices, from prescanned data. No IO.
static void column_raw_with_indices(const PrescannedDat& dat,
    const std::string& label, bool exactMatch,
    std::vector<double>& outValues,
    std::vector<size_t>& outRowIndices)
{
    outValues.clear();
    outRowIndices.clear();
    int colIdx = find_column_index(dat, label, exactMatch);
    if (colIdx < 0) return;
    outValues = dat.valuePerCol[colIdx];
    outRowIndices = dat.rowIdxPerCol[colIdx];
}


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
// Read one column from a CHAOS .dat file as TRULY RAW: only the cells that
// were actually populated, paired with their row indices. No interpolation,
// no gap filling, no modification. Caller converts row index -> seconds
// using the file's row rate.
//
// - If `exactMatch` is true, header names must match `label` exactly
//   (case-insensitive).
// - If false, `contains()` matching is used.
//
// Writes into `outValues` / `outRowIndices` in parallel (index k of each
// vector refers to the same sample). Clears both outputs on entry.
// Leaves both empty if the column is absent or entirely empty.
// ============================================================================
static void read_dat_column_raw_with_indices(
    const std::filesystem::path& path,
    const std::string& label,
    bool exactMatch,
    std::vector<double>& outValues,
    std::vector<size_t>& outRowIndices)
{
    outValues.clear();
    outRowIndices.clear();

    std::ifstream in(path);
    if (!in || label.empty()) return;

    std::vector<std::string> hdrs = find_real_header(in);
    if (hdrs.empty()) return;

    int colIdx = -1;
    for (int i = 0; i < (int)hdrs.size(); ++i) {
        bool match = exactMatch ? equals_ci(hdrs[i], label)
            : contains(hdrs[i], label);
        if (match) { colIdx = i; break; }
    }
    if (colIdx < 0) return;

    outValues.reserve(1 << 18);
    outRowIndices.reserve(1 << 18);

    size_t rowIdx = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::vector<std::string> row = parse_csv_row(line);
        if (colIdx < (int)row.size() && !row[colIdx].empty()) {
            try {
                outValues.push_back(std::stod(row[colIdx]));
                outRowIndices.push_back(rowIdx);
            }
            catch (...) {
                // Skip malformed cells entirely -- we don't want them
                // polluting the raw block with fake values.
            }
        }
        ++rowIdx;
    }
}


// ============================================================================
// Write an EDF channel to binary output as a (upsampled, raw-pairs) pair.
//   - Upsampled block:  native samples resampled to TARGET_RATE (or a copy of
//                       raw if the channel is at/below BOOLEAN_RATE and
//                       `skip_resample` is true).
//   - Raw block:        (timestamp, value) pairs interleaved. EDF samples are
//                       uniformly spaced at `old_rate`, so t_k = k / old_rate.
//                       Storing timestamps keeps the on-disk format uniform
//                       with the CHAOS .dat path, where samples are irregular.
// sizeRawOut is the number of PAIRS (not doubles). Byte length = 2 * sizeRawOut * 8.
// nativeRateOut is the channel's native sampling rate in Hz, or 0.0 if absent.
// If the channel index is invalid, a single -1.0 is written for the upsampled
// slot and a single (-1.0, -1.0) sentinel pair for the raw slot.
// ============================================================================
void edf_to_bin(int handle, int idx, long long n, double old_rate,
    std::ofstream& out, uint32_t& sizeUpOut, uint32_t& sizeRawOut,
    float& nativeRateOut,
    bool skip_resample = false) {
    if (idx < 0 || n <= 0) {
        double dummy = -1.0;
        out.write((char*)&dummy, 8);                 // upsampled placeholder
        sizeUpOut = 1;
        double sentinel[2] = { -1.0, -1.0 };          // raw-pair sentinel
        out.write((char*)sentinel, 16);
        sizeRawOut = 1;                              // 1 pair
        nativeRateOut = 0.0f;
        return;
    }

    // Read native-rate samples.
    std::vector<double> raw(n);
    edfread_physical_samples(handle, idx, (int)n, raw.data());

    // Upsampled block first.
    if (!skip_resample) {
        std::vector<double> up = upsample(raw, old_rate);
        out.write((char*)up.data(), up.size() * 8);
        sizeUpOut = (uint32_t)up.size();
    }
    else {
        // Skip-resample channels: "upsampled" copy is identical to raw.
        out.write((char*)raw.data(), raw.size() * 8);
        sizeUpOut = (uint32_t)raw.size();
    }

    // Raw block: (t, v) pairs, uniformly spaced at old_rate.
    std::vector<double> pairs;
    pairs.reserve(raw.size() * 2);
    const double dt = (old_rate > 0.0) ? (1.0 / old_rate) : 0.0;
    for (size_t k = 0; k < raw.size(); ++k) {
        pairs.push_back(static_cast<double>(k) * dt);
        pairs.push_back(raw[k]);
    }
    out.write((char*)pairs.data(), pairs.size() * 8);
    sizeRawOut = (uint32_t)raw.size();   // pair count
    nativeRateOut = static_cast<float>(old_rate);
}


// ============================================================================
// Write a placeholder for a missing channel: a single -1.0 for the upsampled
// slot and a single (-1.0, -1.0) sentinel pair for the raw slot.
// Also sets native rate to 0.0 (= absent).
// ============================================================================
static void write_missing(std::ofstream& out,
    uint32_t& sizeUpOut, uint32_t& sizeRawOut,
    float& nativeRateOut) {
    double v = -1.0;
    out.write((char*)&v, 8);                         // upsampled placeholder
    sizeUpOut = 1;
    double sentinel[2] = { -1.0, -1.0 };              // raw-pair sentinel
    out.write((char*)sentinel, 16);
    sizeRawOut = 1;                                  // 1 pair
    nativeRateOut = 0.0f;
}


// ============================================================================
// Synthesize and write a timestamp channel for a file whose native samples
// are uniformly spaced (EDF). The upsampled block is t_k = k / TARGET_RATE
// at TARGET_RATE for the total file duration, and the raw block stores
// (k / native_rate, k / native_rate) pairs at the native rate -- trivially
// monotonic but uniform with the rest of the format.
//
// - durationSec: recording duration, in seconds (from the primary channel).
// - nativeRate:  the native rate whose samples anchor the raw timestamps.
// If nativeRate <= 0 or durationSec <= 0, writes a missing-channel placeholder.
// ============================================================================
static void write_synthetic_timestamp(std::ofstream& out,
    double durationSec, double nativeRate,
    uint32_t& sizeUpOut, uint32_t& sizeRawOut,
    float& nativeRateOut) {
    if (durationSec <= 0.0 || nativeRate <= 0.0) {
        write_missing(out, sizeUpOut, sizeRawOut, nativeRateOut);
        return;
    }

    // Upsampled block at TARGET_RATE (= final_sampling_rate)
    const size_t upLen = static_cast<size_t>(
        std::ceil(durationSec * final_sampling_rate));
    std::vector<double> up(upLen);
    const double dtUp = 1.0 / final_sampling_rate;
    for (size_t k = 0; k < upLen; ++k) up[k] = static_cast<double>(k) * dtUp;
    out.write((char*)up.data(), up.size() * 8);
    sizeUpOut = (uint32_t)up.size();

    // Raw block at native rate: (t, t) pairs.
    const size_t rawLen = static_cast<size_t>(
        std::floor(durationSec * nativeRate)) + 1;
    std::vector<double> pairs;
    pairs.reserve(rawLen * 2);
    const double dtRaw = 1.0 / nativeRate;
    for (size_t k = 0; k < rawLen; ++k) {
        const double t = static_cast<double>(k) * dtRaw;
        pairs.push_back(t);
        pairs.push_back(t);
    }
    out.write((char*)pairs.data(), pairs.size() * 8);
    sizeRawOut = (uint32_t)rawLen;
    nativeRateOut = static_cast<float>(nativeRate);
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
    // Bittium records pacemaker events at 8 Hz in a dedicated EDF channel.
    // Label varies across exports -- try a few common spellings before giving up.
    cm.pacemaker = find("Pacemaker");
    if (cm.pacemaker < 0) cm.pacemaker = find("Pace_Event");
    if (cm.pacemaker < 0) cm.pacemaker = find("Pace");

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

    uint32_t sizes_up[NUM_CHANNELS] = {};    // upsampled-block sizes
    uint32_t sizes_raw[NUM_CHANNELS] = {};   // raw-block sizes
    float    native_rates[NUM_CHANNELS] = {}; // per-channel native rate (Hz), 0 = absent

    // -------------------------------------------------------------------------
    // Channel 0: synthetic timestamp (seconds from start).
    // EDF uses uniform sampling, so there's no monitor clock to pull from --
    // we anchor the timestamp raw block to the primary ECG channel's rate so
    // it shares a grid with the strongest signal in the file.
    // -------------------------------------------------------------------------
    double tsNativeRate = (cfg.ecgRate > 0.0) ? cfg.ecgRate
        : edf_channel_rate(hdr.get(), cm.ecg1);
    double tsDurationSec = 0.0;
    if (cm.ecg1 >= 0 && tsNativeRate > 0.0) {
        tsDurationSec = (double)edf_samples(hdr.get(), cm.ecg1) / tsNativeRate;
    }
    write_synthetic_timestamp(out, tsDurationSec, tsNativeRate,
        sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
        native_rates[CH_TIMESTAMP]);

    // Helper: write one EDF channel as a (upsampled, raw) pair.
    // Determines skip_resample automatically: if the channel's native rate
    // is <= BOOLEAN_RATE (1 Hz), the "upsampled" copy is identical to the raw
    // samples (no resampling applied).
    auto writeChannel = [&](int chIdx, int sizeIdx, double rateOverride = 0.0) {
        double rate = (rateOverride > 0.0) ? rateOverride : edf_channel_rate(hdr.get(), chIdx);
        bool skip = (chIdx >= 0) && (rate <= BOOLEAN_RATE);
        edf_to_bin(hdr->handle, chIdx, edf_samples(hdr.get(), chIdx), rate, out,
            sizes_up[sizeIdx], sizes_raw[sizeIdx], native_rates[sizeIdx], skip);
        };

    // 1-3: ECG 1-3, 4: PPG  (upsample from ecgRate/ppgRate)
    writeChannel(cm.ecg1, CH_ECG1, cfg.ecgRate);
    writeChannel(cm.ecg2, CH_ECG2, cfg.ecgRate);
    writeChannel(cm.ecg3, CH_ECG3, cfg.ecgRate);
    writeChannel(cm.ppg, CH_PPG, cfg.ppgRate);

    // 5-7: Accelerometers (upsample from native rate, e.g. 25 Hz in CHAOS)
    writeChannel(cm.accel_x, CH_ACCEL_X);
    writeChannel(cm.accel_y, CH_ACCEL_Y);
    writeChannel(cm.accel_z, CH_ACCEL_Z);

    // 8-10: Marker (1 Hz, no upsample), Temperature (1 Hz, no upsample),
    //       Pacemaker (8 Hz in Bittium, no upsample -- see note below).
    writeChannel(cm.marker, CH_MARKER);
    writeChannel(cm.temp, CH_TEMP);
    // Pacemaker stays at 8 Hz (pacemaker_event_rate in the header). We pass
    // skip_resample via the <= BOOLEAN_RATE branch automatically by overriding
    // the rate to 8 Hz, but writeChannel's auto-skip only triggers at <= 1 Hz.
    // So pacemaker gets upsampled to 1 kHz in the "upsampled" block (fine --
    // it's a pulse train, upsample is essentially hold/linear), and the raw
    // block preserves the true 8 Hz sampling.
    writeChannel(cm.pacemaker, CH_PACEMAKER);

    // 11-12: EOG-L, EOG-R, 13: EMG (256 Hz in MESA -upsample)
    writeChannel(cm.eog_l, CH_EOG_L);
    writeChannel(cm.eog_r, CH_EOG_R);
    writeChannel(cm.emg, CH_EMG);

    // 14-16: EEG 1-3 (256 Hz in MESA -upsample)
    writeChannel(cm.eeg1, CH_EEG1);
    writeChannel(cm.eeg2, CH_EEG2);
    writeChannel(cm.eeg3, CH_EEG3);

    // 17-23: Pres, Flow, Thor, Abdo, Leg, Therm, Pos (32 Hz in MESA -upsample)
    writeChannel(cm.pres, CH_PRES);
    writeChannel(cm.flow, CH_FLOW);
    writeChannel(cm.thor, CH_THOR);
    writeChannel(cm.abdo, CH_ABDO);
    writeChannel(cm.leg, CH_LEG);
    writeChannel(cm.therm, CH_THERM);
    writeChannel(cm.pos, CH_POS);

    // 24-30: Offset channels (1 Hz -no upsample)
    writeChannel(cm.ekg_off, CH_EKG_OFF);
    writeChannel(cm.eog_l_off, CH_EOG_L_OFF);
    writeChannel(cm.eog_r_off, CH_EOG_R_OFF);
    writeChannel(cm.emg_off, CH_EMG_OFF);
    writeChannel(cm.eeg1_off, CH_EEG1_OFF);
    writeChannel(cm.eeg2_off, CH_EEG2_OFF);
    writeChannel(cm.eeg3_off, CH_EEG3_OFF);

    // 31-32: OxStatus (1 Hz -no upsample), SpO2 (1 Hz -no upsample)
    writeChannel(cm.oxstatus, CH_OXSTATUS);
    writeChannel(cm.spo2, CH_SPO2);

    // 33: HR (1 Hz -no upsample)
    writeChannel(cm.hr, CH_HR);

    // 34: DHR (256 Hz in MESA -upsample)
    writeChannel(cm.dhr, CH_DHR);

    // 35: Resp
    writeChannel(cm.resp, CH_RESP);

    // 36: ABP (Bittium only -not in MESA EDF)
    writeChannel(cm.abp, CH_ABP);

    // 37: EEG4 (Bittium only -not in MESA EDF)
    writeChannel(cm.eeg4, CH_EEG4);

    // 38: ART systemic (CHAOS only -not in MESA EDF)
    writeChannel(cm.art, CH_ART);

    // 39: ART_PULM pulmonary (CHAOS only -not in MESA EDF)
    writeChannel(cm.art_pulm, CH_ART_PULM);

    // 40: reserved slot -- always absent, keeps the channel count at 41.
    write_missing(out, sizes_up[CH_RESERVED_40], sizes_raw[CH_RESERVED_40],
        native_rates[CH_RESERVED_40]);

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

    // Upsampled-block sizes (41), raw-block sizes (41), native rates (41 floats).
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&sizes_up[i], 4);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&sizes_raw[i], 4);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&native_rates[i], 4);

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

    uint32_t sizes_up[NUM_CHANNELS] = {};     // upsampled-block sizes
    uint32_t sizes_raw[NUM_CHANNELS] = {};    // raw-block sizes
    float    native_rates[NUM_CHANNELS] = {}; // per-channel native rate (Hz)

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

    // ------------------------------------------------------------------
    // ONE-TIME prescan: parse the .dat file once, cache every column's
    // populated values + row indices in memory. All channel writes below
    // pull from this struct -- no further file IO, no further line parsing.
    // ------------------------------------------------------------------
    std::cout << "  [info] prescanning columns..." << std::flush;
    PrescannedDat prescan = prescan_dat_columns(path);
    std::cout << " done (" << prescan.totalRows << " rows, "
        << prescan.headers.size() << " cols)\n";

    // Helper: read a column two ways, then write it as a (upsampled, raw-pairs) pair.
    //   - Upsampled block:  built from `read_dat_column_interpolated()` (which
    //                       fills empty cells by linear interpolation so the
    //                       upsample has no holes), then resampled from row_rate
    //                       to TARGET_RATE.
    //   - Raw block:        built from `read_dat_column_raw_with_indices()`
    //                       (ONLY the cells that were actually populated, no
    //                       interpolation), emitted as interleaved (t, v) doubles.
    //                       t_k = row_index_k / row_rate, in seconds from start.
    // Writes a -1.0 in the upsampled slot and a (-1.0, -1.0) sentinel pair in
    // the raw slot if the column is absent or completely empty.
    // Helper: resample a channel to TARGET_RATE using ONLY its real populated
    // samples, placed at their true times (rawRowIdx[k] / row_rate seconds).
    // This avoids the "fake plateau" bug where a sparse channel (e.g. PPG
    // populated every 4th row in a 500 Hz file) is incorrectly treated as a
    // dense signal. Linear interpolation between consecutive real samples,
    // matching the style of upsample() in resample.hpp.
    auto resample_from_sparse = [&](const std::vector<double>& rawValues,
        const std::vector<size_t>& rawRowIdx,
        size_t totalRows) -> std::vector<double>
        {
            if (rawValues.size() < 2 || row_rate <= 0.0) return {};

            // Total recording duration, in seconds, shared across all channels
            // (this is what keeps ECG and PPG on the same time axis).
            const double totalDur = static_cast<double>(totalRows) / row_rate;
            const size_t outLen = static_cast<size_t>(
                std::ceil(totalDur * TARGET_RATE));
            if (outLen == 0) return {};

            std::vector<double> out(outLen);

            // Real sample times, in seconds, in ascending order.
            // t[k] = rawRowIdx[k] / row_rate
            const double inv_row_rate = 1.0 / row_rate;
            const size_t N = rawValues.size();

            size_t k = 0;   // index into rawValues/rawRowIdx; advances monotonically
            for (size_t m = 0; m < outLen; ++m) {
                const double t = static_cast<double>(m) / TARGET_RATE;

                // Advance k so that t[k] <= t < t[k+1] (or k == N-1 past the end).
                while (k + 1 < N &&
                    static_cast<double>(rawRowIdx[k + 1]) * inv_row_rate <= t) {
                    ++k;
                }

                if (k + 1 >= N) {
                    // At/past last real sample: hold.
                    out[m] = rawValues[N - 1];
                    continue;
                }

                const double t0 = static_cast<double>(rawRowIdx[k]) * inv_row_rate;
                const double t1 = static_cast<double>(rawRowIdx[k + 1]) * inv_row_rate;

                if (t <= t0) {
                    // Before first real sample: hold first value.
                    out[m] = rawValues[k];
                    continue;
                }

                const double span = t1 - t0;
                const double f = (span > 0.0) ? (t - t0) / span : 0.0;
                out[m] = rawValues[k] * (1.0 - f) + rawValues[k + 1] * f;
            }
            return out;
        };

    // Helper that infers a column's native rate from its populated-cell
    // density. For a dense channel (ECG, populated every row), this returns
    // ~row_rate. For a sparse channel (PPG, populated every 4th row), it
    // returns ~row_rate / 4 -- i.e. the channel's real native sampling rate,
    // not the file's row rate. This matches the per-channel "Original
    // Sampling Rate" that downstream consumers expect.
    auto infer_native_rate = [&](const std::vector<size_t>& rawRowIdx) -> double {
        if (rawRowIdx.size() < 2 || row_rate <= 0.0) return 0.0;
        const size_t span = rawRowIdx.back() - rawRowIdx.front();
        if (span == 0) return 0.0;
        // Average rows-between-samples over the whole column.
        const double rowsPerSample =
            static_cast<double>(span) / (rawRowIdx.size() - 1);
        if (rowsPerSample <= 0.0) return 0.0;
        return row_rate / rowsPerSample;
        };

    auto writeCol = [&](const std::string& label, bool exact,
        uint32_t& sizeUpOut, uint32_t& sizeRawOut, float& nativeRateOut) {
            // Truly raw: only populated cells + their row indices.
            std::vector<double> rawValues;
            std::vector<size_t> rawRowIdx;
            column_raw_with_indices(prescan, label, exact,
                rawValues, rawRowIdx);

            if (rawValues.empty()) {
                double v = -1.0;
                out.write((char*)&v, 8);                     // upsampled placeholder
                sizeUpOut = 1;
                double sentinel[2] = { -1.0, -1.0 };          // raw-pair sentinel
                out.write((char*)sentinel, 16);
                sizeRawOut = 1;
                nativeRateOut = 0.0f;
                return;
            }

            // Total row count comes from the prescan -- no second file pass.
            const size_t totalRows = prescan.totalRows;

            // Upsampled block: built from the REAL populated samples at their
            // REAL times -- no fake "one value per row" assumption. This keeps
            // ECG (dense) and PPG (sparse) on a shared wall-clock time axis.
            std::vector<double> up =
                resample_from_sparse(rawValues, rawRowIdx, totalRows);

            if (up.empty()) {
                // Fallback: not enough real samples to resample meaningfully.
                double v = -1.0;
                out.write((char*)&v, 8);
                sizeUpOut = 1;
            }
            else {
                out.write((char*)up.data(), up.size() * 8);
                sizeUpOut = (uint32_t)up.size();
            }

            // Raw block: interleaved (t, v) pairs at the channel's actual
            // populated-cell times. Irregular spacing, dropouts, and jitter
            // are all preserved exactly.
            std::vector<double> pairs;
            pairs.reserve(rawValues.size() * 2);
            const double dt = (row_rate > 0.0) ? (1.0 / row_rate) : 0.0;
            for (size_t k = 0; k < rawValues.size(); ++k) {
                pairs.push_back(static_cast<double>(rawRowIdx[k]) * dt);
                pairs.push_back(rawValues[k]);
            }
            out.write((char*)pairs.data(), pairs.size() * 8);
            sizeRawOut = (uint32_t)rawValues.size();   // pair count
            nativeRateOut = static_cast<float>(infer_native_rate(rawRowIdx));
        };

    // -------------------------------------------------------------------------
    // 0: Timestamp channel. The CHAOS Monitor TimeStamp is conceptually a
    // 500 Hz column populated every row, but we store seconds-from-start as
    // the value (not the raw wall-clock string). Row rate is authoritative.
    // -------------------------------------------------------------------------
    {
        // Use prescan totalRows -- no extra file IO.
        const size_t totalRows = prescan.totalRows;
        const double durationSec = (row_rate > 0.0)
            ? static_cast<double>(totalRows) / row_rate
            : 0.0;
        write_synthetic_timestamp(out, durationSec, row_rate,
            sizes_up[CH_TIMESTAMP], sizes_raw[CH_TIMESTAMP],
            native_rates[CH_TIMESTAMP]);
    }

    // 1-3: ECG 1-3, 4: PPG
    writeCol(cfg.ecg1Label, false, sizes_up[CH_ECG1], sizes_raw[CH_ECG1], native_rates[CH_ECG1]);
    writeCol(cfg.ecg2Label, false, sizes_up[CH_ECG2], sizes_raw[CH_ECG2], native_rates[CH_ECG2]);
    writeCol(cfg.ecg3Label, false, sizes_up[CH_ECG3], sizes_raw[CH_ECG3], native_rates[CH_ECG3]);
    writeCol(cfg.ppgLabel, false, sizes_up[CH_PPG], sizes_raw[CH_PPG], native_rates[CH_PPG]);

    // 5-7: Accelerometers (not present in CHAOS .dat)
    write_missing(out, sizes_up[CH_ACCEL_X], sizes_raw[CH_ACCEL_X], native_rates[CH_ACCEL_X]);
    write_missing(out, sizes_up[CH_ACCEL_Y], sizes_raw[CH_ACCEL_Y], native_rates[CH_ACCEL_Y]);
    write_missing(out, sizes_up[CH_ACCEL_Z], sizes_raw[CH_ACCEL_Z], native_rates[CH_ACCEL_Z]);

    // 8-10: Marker, Temp, Pacemaker (not present)
    write_missing(out, sizes_up[CH_MARKER], sizes_raw[CH_MARKER], native_rates[CH_MARKER]);
    write_missing(out, sizes_up[CH_TEMP], sizes_raw[CH_TEMP], native_rates[CH_TEMP]);
    write_missing(out, sizes_up[CH_PACEMAKER], sizes_raw[CH_PACEMAKER], native_rates[CH_PACEMAKER]);

    // 11-13: EOG-L, EOG-R, EMG (not present)
    write_missing(out, sizes_up[CH_EOG_L], sizes_raw[CH_EOG_L], native_rates[CH_EOG_L]);
    write_missing(out, sizes_up[CH_EOG_R], sizes_raw[CH_EOG_R], native_rates[CH_EOG_R]);
    write_missing(out, sizes_up[CH_EMG], sizes_raw[CH_EMG], native_rates[CH_EMG]);

    // 14-16: EEG 1-3 (contains-match; the real columns have _LBL suffix)
    writeCol("NLS_EEG_NAMES_EEG_CHAN1", false, sizes_up[CH_EEG1], sizes_raw[CH_EEG1], native_rates[CH_EEG1]);
    writeCol("NLS_EEG_NAMES_EEG_CHAN2", false, sizes_up[CH_EEG2], sizes_raw[CH_EEG2], native_rates[CH_EEG2]);
    writeCol("NLS_EEG_NAMES_EEG_CHAN3", false, sizes_up[CH_EEG3], sizes_raw[CH_EEG3], native_rates[CH_EEG3]);

    // 17: Pres / CVP
    writeCol("NLS_NOM_PRESS_BLD_VEN_CENT", false,
        sizes_up[CH_PRES], sizes_raw[CH_PRES], native_rates[CH_PRES]);

    // 18-23: Flow, Thor, Abdo, Leg, Therm, Pos (not present)
    for (int i = CH_FLOW; i <= CH_POS; ++i)
        write_missing(out, sizes_up[i], sizes_raw[i], native_rates[i]);

    // 24-30: MESA offset channels (not present)
    for (int i = CH_EKG_OFF; i <= CH_EEG3_OFF; ++i)
        write_missing(out, sizes_up[i], sizes_raw[i], native_rates[i]);

    // 31-34: OxStatus, SpO2, HR, DHR (not present)
    for (int i = CH_OXSTATUS; i <= CH_DHR; ++i)
        write_missing(out, sizes_up[i], sizes_raw[i], native_rates[i]);

    // 35: Resp
    writeCol("NLS_NOM_RESP", false,
        sizes_up[CH_RESP], sizes_raw[CH_RESP], native_rates[CH_RESP]);
    // 36: ABP (exact match -- distinguishes from ART / ART_PULM)
    writeCol("NLS_NOM_PRESS_BLD_ART_ABP", true,
        sizes_up[CH_ABP], sizes_raw[CH_ABP], native_rates[CH_ABP]);
    // 37: EEG4 (contains-match)
    writeCol("NLS_EEG_NAMES_EEG_CHAN4", false,
        sizes_up[CH_EEG4], sizes_raw[CH_EEG4], native_rates[CH_EEG4]);
    // 38: ART systemic arterial (exact match)
    writeCol("NLS_NOM_PRESS_BLD_ART", true,
        sizes_up[CH_ART], sizes_raw[CH_ART], native_rates[CH_ART]);
    // 39: ART_PULM pulmonary arterial (exact match)
    writeCol("NLS_NOM_PRESS_BLD_ART_PULM", true,
        sizes_up[CH_ART_PULM], sizes_raw[CH_ART_PULM], native_rates[CH_ART_PULM]);

    // 40: reserved slot
    write_missing(out, sizes_up[CH_RESERVED_40], sizes_raw[CH_RESERVED_40],
        native_rates[CH_RESERVED_40]);

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

    // Upsampled-block sizes (41), raw-block sizes (41), native rates (41 floats).
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&sizes_up[i], 4);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&sizes_raw[i], 4);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        out.write((char*)&native_rates[i], 4);

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