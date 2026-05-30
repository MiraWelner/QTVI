#include "theme.h"

#include <QApplication>
#include <QFile>
#include <QPalette>

namespace Theme {

    void apply(QApplication& app)
    {
        app.setStyle("Fusion");

        QPalette p;

        const QColor bg(250, 250, 249);       // warm off-white canvas
        const QColor surface(255, 255, 255);  // cards / inputs
        const QColor text(28, 28, 34);         // near-black, easy on the eye
        const QColor mutedText(120, 120, 130); // secondary labels
        const QColor accent(99, 91, 255);      // confident indigo
        const QColor border(226, 226, 231);

        // Surfaces
        p.setColor(QPalette::Window, bg);
        p.setColor(QPalette::Base, surface);
        p.setColor(QPalette::AlternateBase, QColor(247, 247, 248));
        p.setColor(QPalette::Button, surface);

        // Text
        p.setColor(QPalette::WindowText, text);
        p.setColor(QPalette::Text, text);
        p.setColor(QPalette::ButtonText, text);
        p.setColor(QPalette::PlaceholderText, mutedText);
        p.setColor(QPalette::BrightText, Qt::white);

        // Accent / selection
        p.setColor(QPalette::Highlight, accent);
        p.setColor(QPalette::HighlightedText, Qt::white);
        p.setColor(QPalette::Link, accent);
        p.setColor(QPalette::LinkVisited, QColor(82, 72, 235));

        // Tooltips
        p.setColor(QPalette::ToolTipBase, text);
        p.setColor(QPalette::ToolTipText, Qt::white);

        // Disabled states — keep things looking intentional, not broken
        p.setColor(QPalette::Disabled, QPalette::Text, QColor(176, 176, 184));
        p.setColor(QPalette::Disabled, QPalette::ButtonText, QColor(176, 176, 184));
        p.setColor(QPalette::Disabled, QPalette::WindowText, QColor(176, 176, 184));
        p.setColor(QPalette::Disabled, QPalette::Highlight, QColor(220, 220, 226));

        Q_UNUSED(border);
        app.setPalette(p);

        QFile file(":/theme/style.qss");
     

        if (file.open(QFile::ReadOnly | QFile::Text)) {
            app.setStyleSheet(file.readAll());
        }
        app.setFont(QFont("Inter", 8));
    }

    QFont Theme::chartTitleFont() {
        return QFont("Inter", 8, QFont::Bold);
    }

    QFont Theme::chartAxisFont() {
        return QFont("Inter", 7);
    }

}