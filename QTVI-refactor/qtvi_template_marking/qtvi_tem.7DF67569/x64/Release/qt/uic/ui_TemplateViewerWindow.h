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
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
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
    QLabel *instructionsLabel;
    QHBoxLayout *topBar;
    QLabel *subjectLabel;
    QCheckBox *show_ecg;
    QCheckBox *show_ppg;
    QSpacerItem *sL;
    QPushButton *prevButton;
    QLabel *pageLabel;
    QPushButton *nextButton;
    QSpacerItem *sR;
    QRadioButton *MoveIndividual;
    QRadioButton *MoveSubsequent;
    QPushButton *finishButton;
    QScrollArea *scrollArea;
    QWidget *scrollContents;
    QGridLayout *plotGrid;
    QButtonGroup *buttonGroup;

    void setupUi(QMainWindow *TemplateViewerWindow)
    {
        if (TemplateViewerWindow->objectName().isEmpty())
            TemplateViewerWindow->setObjectName("TemplateViewerWindow");
        TemplateViewerWindow->resize(1500, 1000);
        actionMoveSubsequent = new QAction(TemplateViewerWindow);
        actionMoveSubsequent->setObjectName("actionMoveSubsequent");
        actionMoveSubsequent->setCheckable(true);
        actionMoveSubsequent->setChecked(true);
        centralwidget = new QWidget(TemplateViewerWindow);
        centralwidget->setObjectName("centralwidget");
        mainLayout = new QVBoxLayout(centralwidget);
        mainLayout->setSpacing(4);
        mainLayout->setObjectName("mainLayout");
        mainLayout->setContentsMargins(4, 4, 4, 4);
        instructionsLabel = new QLabel(centralwidget);
        instructionsLabel->setObjectName("instructionsLabel");

        mainLayout->addWidget(instructionsLabel);

        topBar = new QHBoxLayout();
        topBar->setObjectName("topBar");
        subjectLabel = new QLabel(centralwidget);
        subjectLabel->setObjectName("subjectLabel");
        QFont font;
        font.setBold(true);
        subjectLabel->setFont(font);

        topBar->addWidget(subjectLabel);

        show_ecg = new QCheckBox(centralwidget);
        show_ecg->setObjectName("show_ecg");
        show_ecg->setChecked(true);

        topBar->addWidget(show_ecg);

        show_ppg = new QCheckBox(centralwidget);
        show_ppg->setObjectName("show_ppg");
        show_ppg->setChecked(true);

        topBar->addWidget(show_ppg);

        sL = new QSpacerItem(0, 0, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        topBar->addItem(sL);

        prevButton = new QPushButton(centralwidget);
        prevButton->setObjectName("prevButton");
        prevButton->setEnabled(false);
        prevButton->setMinimumSize(QSize(80, 0));

        topBar->addWidget(prevButton);

        pageLabel = new QLabel(centralwidget);
        pageLabel->setObjectName("pageLabel");
        pageLabel->setMinimumSize(QSize(80, 0));
        pageLabel->setAlignment(Qt::AlignmentFlag::AlignCenter);

        topBar->addWidget(pageLabel);

        nextButton = new QPushButton(centralwidget);
        nextButton->setObjectName("nextButton");
        nextButton->setEnabled(false);
        nextButton->setMinimumSize(QSize(80, 0));

        topBar->addWidget(nextButton);

        sR = new QSpacerItem(0, 0, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        topBar->addItem(sR);

        MoveIndividual = new QRadioButton(centralwidget);
        buttonGroup = new QButtonGroup(TemplateViewerWindow);
        buttonGroup->setObjectName("buttonGroup");
        buttonGroup->addButton(MoveIndividual);
        MoveIndividual->setObjectName("MoveIndividual");

        topBar->addWidget(MoveIndividual);

        MoveSubsequent = new QRadioButton(centralwidget);
        buttonGroup->addButton(MoveSubsequent);
        MoveSubsequent->setObjectName("MoveSubsequent");
        MoveSubsequent->setChecked(true);

        topBar->addWidget(MoveSubsequent);

        finishButton = new QPushButton(centralwidget);
        finishButton->setObjectName("finishButton");

        topBar->addWidget(finishButton);


        mainLayout->addLayout(topBar);

        scrollArea = new QScrollArea(centralwidget);
        scrollArea->setObjectName("scrollArea");
        scrollArea->setWidgetResizable(true);
        scrollContents = new QWidget();
        scrollContents->setObjectName("scrollContents");
        scrollContents->setGeometry(QRect(0, 0, 1490, 940));
        plotGrid = new QGridLayout(scrollContents);
        plotGrid->setSpacing(2);
        plotGrid->setObjectName("plotGrid");
        plotGrid->setContentsMargins(2, 2, 2, 2);
        scrollArea->setWidget(scrollContents);

        mainLayout->addWidget(scrollArea);

        TemplateViewerWindow->setCentralWidget(centralwidget);

        retranslateUi(TemplateViewerWindow);
        QObject::connect(finishButton, SIGNAL(clicked()), TemplateViewerWindow, SLOT(onFinish()));
        QObject::connect(prevButton, SIGNAL(clicked()), TemplateViewerWindow, SLOT(onPrevPage()));
        QObject::connect(nextButton, SIGNAL(clicked()), TemplateViewerWindow, SLOT(onNextPage()));

        finishButton->setDefault(true);


        QMetaObject::connectSlotsByName(TemplateViewerWindow);
    } // setupUi

    void retranslateUi(QMainWindow *TemplateViewerWindow)
    {
        TemplateViewerWindow->setWindowTitle(QCoreApplication::translate("TemplateViewerWindow", "Template Marking Viewer", nullptr));
        actionMoveSubsequent->setText(QCoreApplication::translate("TemplateViewerWindow", "Move Subsequent Dicrotic", nullptr));
        instructionsLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "Single Right Click: Marks bad R peak                Double Right Click: Marks bad PPG if PPG present (throw away that template)           Click and Drag Markers to Move", nullptr));
        subjectLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "No subject loaded", nullptr));
        show_ecg->setText(QCoreApplication::translate("TemplateViewerWindow", "Show ECG Markers", nullptr));
        show_ppg->setText(QCoreApplication::translate("TemplateViewerWindow", "Show PPG Markers", nullptr));
        prevButton->setText(QCoreApplication::translate("TemplateViewerWindow", "<  Prev", nullptr));
#if QT_CONFIG(shortcut)
        prevButton->setShortcut(QCoreApplication::translate("TemplateViewerWindow", "Left", nullptr));
#endif // QT_CONFIG(shortcut)
        pageLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "1 / 1", nullptr));
        nextButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Next  >", nullptr));
#if QT_CONFIG(shortcut)
        nextButton->setShortcut(QCoreApplication::translate("TemplateViewerWindow", "Right", nullptr));
#endif // QT_CONFIG(shortcut)
        MoveIndividual->setText(QCoreApplication::translate("TemplateViewerWindow", "Move Individual", nullptr));
        MoveSubsequent->setText(QCoreApplication::translate("TemplateViewerWindow", "Move Subsequent", nullptr));
        finishButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Finish && Save", nullptr));
#if QT_CONFIG(shortcut)
        finishButton->setShortcut(QCoreApplication::translate("TemplateViewerWindow", "Return", nullptr));
#endif // QT_CONFIG(shortcut)
    } // retranslateUi

};

namespace Ui {
    class TemplateViewerWindow: public Ui_TemplateViewerWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TEMPLATEVIEWERWINDOW_H
