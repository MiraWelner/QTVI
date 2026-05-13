#pragma once

#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <vector>
#include <optional>

struct DatasetConfig {
    QString dataType;          // col 0   DATA_TYPE
    QString templatePath;      // col 12  template_path
    QString markingPath;       // col 13  markings_path
    double  ecgSamplingRate;   // col 3   ecg_rate
    double  ppgSamplingRate;   // col 4   ppg_rate     (0 if absent)
    double  upsampledRate;     // col 5   upsampled_rate
    int     binMinutes;        // col 6   bin_size_minutes
};

namespace config_detail {
    inline double toDouble(const QString& s) {
        bool ok = false;
        double v = s.trimmed().toDouble(&ok);
        return ok ? v : 0.0;
    }
    // Safe column access: returns "" if the row is short.
    inline QString col(const QStringList& cols, int i) {
        return (i < cols.size()) ? cols[i].trimmed() : QString();
    }
}

inline std::vector<DatasetConfig> readConfig(const QString& csvPath) {
    using namespace config_detail;

    std::vector<DatasetConfig> result;
    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return result;

    QTextStream in(&file);
    bool first = true;
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (first) { first = false; continue; }   // skip header
        if (line.isEmpty()) continue;

        QStringList cols = line.split(',');

        DatasetConfig cfg;
        cfg.dataType = col(cols, 0);
        cfg.ecgSamplingRate = toDouble(col(cols, 3));
        cfg.ppgSamplingRate = toDouble(col(cols, 4));
        cfg.upsampledRate = toDouble(col(cols, 5));
        cfg.binMinutes = static_cast<int>(toDouble(col(cols, 6)));
        cfg.templatePath = col(cols, 12);
        cfg.markingPath = col(cols, 13);
        result.push_back(cfg);
    }
    return result;
}

inline std::optional<DatasetConfig> findConfig(const std::vector<DatasetConfig>& configs,
    const QString& dataType) {
    for (const auto& c : configs) {
        if (c.dataType.compare(dataType, Qt::CaseInsensitive) == 0)
            return c;
    }
    return std::nullopt;
}