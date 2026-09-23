#pragma once

#include <QString>
#include <QStringList>
#include <optional>

struct TemplateFileName {
    QString subjectId;
    QString date;       // empty if the source filename had no date field
    int     rateHz;
    int     binMinutes;
};

// Parse a template-file stem (caller has already stripped the "_templates"
// suffix). The trailing two underscore-separated tokens are always
// {rate}_{binMinutes} and must both parse as ints. What precedes them
// depends on the dataset:
//
//   MESA-style    (4+ tokens):  "{subjectId}_{date}_{rate}_{binMinutes}"
//                                e.g. "3010023_20110817_1000_005"
//                                subjectId itself may contain underscores;
//                                everything before {date} is rejoined.
//
//   Bittium-style (3 tokens):   "{subjectId}_{rate}_{binMinutes}"
//                                e.g. "001-09-40-14_1000_005"
//                                No date field; `date` stays empty.
//                                The subject id may contain dashes -- those
//                                don't split here because we split on '_' only.
//
// Earlier this function rejected anything with fewer than 4 tokens, which
// silently dropped every Bittium file as "unparseable" in template_marking.
inline std::optional<TemplateFileName> parseTemplateFileName(const QString& stem) {
    QStringList parts = stem.split('_');
    if (parts.size() < 3) return std::nullopt;

    bool ok1 = false, ok2 = false;
    int rate = parts[parts.size() - 2].toInt(&ok1);
    int binMin = parts[parts.size() - 1].toInt(&ok2);
    if (!ok1 || !ok2) return std::nullopt;

    TemplateFileName out;
    out.rateHz = rate;
    out.binMinutes = binMin;

    if (parts.size() == 3) {
        // Bittium-style: {subjectId}_{rate}_{binMinutes}
        out.subjectId = parts[0];
        out.date = QString();
    }
    else {
        // MESA-style: {subjectId...}_{date}_{rate}_{binMinutes}
        out.date = parts[parts.size() - 3];
        QStringList subjectParts = parts.mid(0, parts.size() - 3);
        if (subjectParts.isEmpty()) return std::nullopt;
        out.subjectId = subjectParts.join('_');
    }
    return out;
}