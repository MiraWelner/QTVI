#pragma once

#include <QString>
#include <QStringList>
#include <QFile>
#include <QTextStream>
#include <vector>
#include <optional>

struct DatasetConfig {
    QString dataType;         // col 0
    QString templatePath;     // col 8
    QString markingPath;      // col 9
    double ecgSamplingRate;   // col 14
    double ppgSamplingRate;   // col 15
};

inline std::vector<DatasetConfig> readConfig(const QString& csvPath) {
    std::vector<DatasetConfig> result;
    QFile file(csvPath);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) return result;

    QTextStream in(&file);
    bool first = true;
    while (!in.atEnd()) {
        QString line = in.readLine().trimmed();
        if (first) { first = false; continue; }
        if (line.isEmpty()) continue;

        QStringList cols = line.split(',');
        if (cols.size() < 16) continue;

        DatasetConfig cfg;
        cfg.dataType = cols[0].trimmed();
        cfg.templatePath = cols[8].trimmed();
        cfg.markingPath = cols[9].trimmed();
        cfg.ecgSamplingRate = cols[14].trimmed().toDouble();
        cfg.ppgSamplingRate = cols[15].trimmed().toDouble();
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