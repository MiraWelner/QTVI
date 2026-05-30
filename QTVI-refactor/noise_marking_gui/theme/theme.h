#pragma once

class QApplication;

#include <QFont>

namespace Theme {
    void apply(QApplication& app);

    QFont chartTitleFont();   // small, bold — chart titles
    QFont chartAxisFont();    // smaller, regular — tick labels
}

namespace Theme {
    void apply(QApplication& app);
}