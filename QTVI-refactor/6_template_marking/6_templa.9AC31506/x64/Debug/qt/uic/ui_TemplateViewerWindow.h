/********************************************************************************
** Form generated from reading UI file 'TemplateViewerWindow.ui'
**
** Created by: Qt User Interface Compiler version 6.10.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_TEMPLATEVIEWERWINDOW_H
#define UI_TEMPLATEVIEWERWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_TemplateViewerWindow
{
public:
    QAction *actionMoveSubsequent;
    QWidget *centralwidget;
    QVBoxLayout *mainLayout;
    QHBoxLayout *topBar;
    QLabel *statusLabel;
    QSpacerItem *topSpacer;
    QLabel *subjectLabel;
    QSpacerItem *topSpacer2;
    QPushButton *finishButton;
    QScrollArea *scrollArea;
    QWidget *scrollContents;
    QVBoxLayout *plotLayout;
    QMenuBar *menubar;
    QMenu *markersMenu;

    void setupUi(QMainWindow *TemplateViewerWindow)
    {
        if (TemplateViewerWindow->objectName().isEmpty())
            TemplateViewerWindow->setObjectName("TemplateViewerWindow");
        TemplateViewerWindow->resize(1200, 900);
        actionMoveSubsequent = new QAction(TemplateViewerWindow);
        actionMoveSubsequent->setObjectName("actionMoveSubsequent");
        actionMoveSubsequent->setCheckable(true);
        centralwidget = new QWidget(TemplateViewerWindow);
        centralwidget->setObjectName("centralwidget");
        mainLayout = new QVBoxLayout(centralwidget);
        mainLayout->setObjectName("mainLayout");
        topBar = new QHBoxLayout();
        topBar->setObjectName("topBar");
        statusLabel = new QLabel(centralwidget);
        statusLabel->setObjectName("statusLabel");

        topBar->addWidget(statusLabel);

        topSpacer = new QSpacerItem(0, 0, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        topBar->addItem(topSpacer);

        subjectLabel = new QLabel(centralwidget);
        subjectLabel->setObjectName("subjectLabel");

        topBar->addWidget(subjectLabel);

        topSpacer2 = new QSpacerItem(0, 0, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        topBar->addItem(topSpacer2);

        finishButton = new QPushButton(centralwidget);
        finishButton->setObjectName("finishButton");

        topBar->addWidget(finishButton);


        mainLayout->addLayout(topBar);

        scrollArea = new QScrollArea(centralwidget);
        scrollArea->setObjectName("scrollArea");
        scrollArea->setWidgetResizable(true);
        scrollContents = new QWidget();
        scrollContents->setObjectName("scrollContents");
        plotLayout = new QVBoxLayout(scrollContents);
        plotLayout->setObjectName("plotLayout");
        scrollArea->setWidget(scrollContents);

        mainLayout->addWidget(scrollArea);

        TemplateViewerWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(TemplateViewerWindow);
        menubar->setObjectName("menubar");
        markersMenu = new QMenu(menubar);
        markersMenu->setObjectName("markersMenu");
        TemplateViewerWindow->setMenuBar(menubar);

        menubar->addAction(markersMenu->menuAction());
        markersMenu->addAction(actionMoveSubsequent);

        retranslateUi(TemplateViewerWindow);

        QMetaObject::connectSlotsByName(TemplateViewerWindow);
    } // setupUi

    void retranslateUi(QMainWindow *TemplateViewerWindow)
    {
        TemplateViewerWindow->setWindowTitle(QCoreApplication::translate("TemplateViewerWindow", "Template Marking Viewer", nullptr));
        actionMoveSubsequent->setText(QCoreApplication::translate("TemplateViewerWindow", "Move Subsequent Dicrotic", nullptr));
        statusLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "Mode: Move Individual", nullptr));
        subjectLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "No subject loaded", nullptr));
        finishButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Finish && Save", nullptr));
        markersMenu->setTitle(QCoreApplication::translate("TemplateViewerWindow", "Markers", nullptr));
    } // retranslateUi

};

namespace Ui {
    class TemplateViewerWindow: public Ui_TemplateViewerWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TEMPLATEVIEWERWINDOW_H
