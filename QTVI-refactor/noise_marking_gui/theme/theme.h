#pragma once

class QApplication;

#include <QFont>

namespace Theme {
    void apply(QApplication& app);

    QFont chartTitleFont();
    QFont chartAxisFont();
}

namespace Theme {
    void apply(QApplication& app);
}