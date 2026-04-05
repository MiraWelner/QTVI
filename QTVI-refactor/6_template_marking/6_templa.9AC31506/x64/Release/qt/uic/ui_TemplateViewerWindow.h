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
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
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
    QLabel *label;
    QHBoxLayout *topBar;
    QLabel *subjectLabel;
    QSpacerItem *sL;
    QPushButton *prevButton;
    QLabel *pageLabel;
    QPushButton *nextButton;
    QSpacerItem *sR;
    QPushButton *modeButton;
    QPushButton *finishButton;
    QScrollArea *scrollArea;
    QWidget *scrollContents;
    QGridLayout *plotGrid;

    void setupUi(QMainWindow *TemplateViewerWindow)
    {
        if (TemplateViewerWindow->objectName().isEmpty())
            TemplateViewerWindow->setObjectName("TemplateViewerWindow");
        TemplateViewerWindow->resize(1200, 900);
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
        label = new QLabel(centralwidget);
        label->setObjectName("label");

        mainLayout->addWidget(label);

        topBar = new QHBoxLayout();
        topBar->setObjectName("topBar");
        subjectLabel = new QLabel(centralwidget);
        subjectLabel->setObjectName("subjectLabel");
        QFont font;
        font.setBold(true);
        subjectLabel->setFont(font);

        topBar->addWidget(subjectLabel);

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

        modeButton = new QPushButton(centralwidget);
        modeButton->setObjectName("modeButton");
        modeButton->setMinimumSize(QSize(190, 0));
        modeButton->setCheckable(true);
        modeButton->setChecked(true);

        topBar->addWidget(modeButton);

        finishButton = new QPushButton(centralwidget);
        finishButton->setObjectName("finishButton");

        topBar->addWidget(finishButton);


        mainLayout->addLayout(topBar);

        scrollArea = new QScrollArea(centralwidget);
        scrollArea->setObjectName("scrollArea");
        scrollArea->setWidgetResizable(true);
        scrollContents = new QWidget();
        scrollContents->setObjectName("scrollContents");
        scrollContents->setGeometry(QRect(0, 0, 1190, 840));
        plotGrid = new QGridLayout(scrollContents);
        plotGrid->setSpacing(2);
        plotGrid->setObjectName("plotGrid");
        plotGrid->setContentsMargins(2, 2, 2, 2);
        scrollArea->setWidget(scrollContents);

        mainLayout->addWidget(scrollArea);

        TemplateViewerWindow->setCentralWidget(centralwidget);

        retranslateUi(TemplateViewerWindow);

        QMetaObject::connectSlotsByName(TemplateViewerWindow);
    } // setupUi

    void retranslateUi(QMainWindow *TemplateViewerWindow)
    {
        TemplateViewerWindow->setWindowTitle(QCoreApplication::translate("TemplateViewerWindow", "Template Marking Viewer", nullptr));
        actionMoveSubsequent->setText(QCoreApplication::translate("TemplateViewerWindow", "Move Subsequent Dicrotic", nullptr));
        label->setText(QCoreApplication::translate("TemplateViewerWindow", "Simple Right Click: Marks bad R peak                Double Right Click: Marks bad PPG if PPG present (throw away that template)           The red line marks the dicrotic notch, click and drag to move", nullptr));
        subjectLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "No subject loaded", nullptr));
        prevButton->setText(QCoreApplication::translate("TemplateViewerWindow", "<  Prev", nullptr));
        pageLabel->setText(QCoreApplication::translate("TemplateViewerWindow", "1 / 1", nullptr));
        nextButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Next  >", nullptr));
        modeButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Mode: Move Subsequent", nullptr));
        finishButton->setText(QCoreApplication::translate("TemplateViewerWindow", "Finish && Save", nullptr));
    } // retranslateUi

};

namespace Ui {
    class TemplateViewerWindow: public Ui_TemplateViewerWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_TEMPLATEVIEWERWINDOW_H
