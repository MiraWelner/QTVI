/**
* @file   beat_log.cpp
* @brief  See beat_log.hpp.
*/
#include "beat_log.hpp"

#include <fstream>
#include <iomanip>

bool beat_log::writeCsv(const std::string& path) const {
    std::ofstream f(path);
    if (!f.is_open()) return false;

    f << "beat,ecg1_x,ecg1_y,ecg2_x,ecg2_y,ecg3_x,ecg3_y,ppg_x,ppg_y,abp_x,abp_y,"
        "blanking_ecg1,threshold_ecg1,blanking_ecg2,threshold_ecg2,"
        "blanking_ecg3,threshold_ecg3,blanking_ppg,threshold_ppg,blanking_abp,threshold_abp\n";
    f << std::fixed << std::setprecision(6);

    for (size_t i = 0; i < m_beats.size(); ++i) {
        f << i;
        for (const sample& s : m_beats[i].chan)      // 10 cols: x,y per channel
            f << ',' << s.x << ',' << s.y;
        for (const sample& s : m_beats[i].chan)      // 10 cols: blanking,threshold per channel
            f << ',' << s.blanking << ',' << s.threshold;
        f << '\n';
    }
    return true;
}