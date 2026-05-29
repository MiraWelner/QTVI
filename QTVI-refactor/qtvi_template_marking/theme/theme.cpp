#include "theme.h"

#include <QApplication>
#include <QFile>
#include <QPalette>

namespace Theme {

    void apply(QApplication& app)
    {
        app.setStyle("Fusion");

        QPalette p;

        QColor bg(248, 248, 246);
        QColor surface(255, 255, 252);
        QColor text(24, 24, 26);
        QColor accent(90, 95, 255);

        p.setColor(QPalette::Window, bg);
        p.setColor(QPalette::Base, surface);

        p.setColor(QPalette::WindowText, text);
        p.setColor(QPalette::Text, text);

        p.setColor(QPalette::Highlight, accent);
        p.setColor(QPalette::HighlightedText, Qt::white);

        app.setPalette(p);

        QFile file("theme/style.qss");

        if (file.open(QFile::ReadOnly | QFile::Text)) {
            app.setStyleSheet(file.readAll());
        }
    }

}