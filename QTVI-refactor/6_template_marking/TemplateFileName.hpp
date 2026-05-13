#pragma once

#include <QString>
#include <QStringList>
#include <optional>

struct TemplateFileName {
    QString subjectId;
    QString date;
    int     rateHz;
    int     binMinutes;
};

// Parse a stem like "3010023_20110817_1000_005_template_info"
// (suffix already stripped by caller — pass the QFileInfo::baseName-ish stem
// minus the trailing "_template_info").
//
// Layout is fixed: {subjectId}_{date}_{rate}_{binMinutes}
// subjectId itself may contain underscores — we take the LAST three
// underscore-separated tokens as date/rate/bin and everything before
// as the subject id.
inline std::optional<TemplateFileName> parseTemplateFileName(const QString& stem) {
    QStringList parts = stem.split('_');
    if (parts.size() < 4) return std::nullopt;

    bool ok1 = false, ok2 = false;
    int rate = parts[parts.size() - 2].toInt(&ok1);
    int binMin = parts[parts.size() - 1].toInt(&ok2);
    if (!ok1 || !ok2) return std::nullopt;

    QString date = parts[parts.size() - 3];

    QStringList subjectParts = parts.mid(0, parts.size() - 3);
    if (subjectParts.isEmpty()) return std::nullopt;

    TemplateFileName out;
    out.subjectId = subjectParts.join('_');
    out.date = date;
    out.rateHz = rate;
    out.binMinutes = binMin;
    return out;
}