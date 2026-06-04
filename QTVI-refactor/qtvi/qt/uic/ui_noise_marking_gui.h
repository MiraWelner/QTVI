/********************************************************************************
** Form generated from reading UI file 'noise_marking_gui.ui'
**
** Created by: Qt User Interface Compiler version 6.10.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_NOISE_MARKING_GUI_H
#define UI_NOISE_MARKING_GUI_H

#include <QtCharts/QChartView>
#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_noise_marking_gui
{
public:
    QHBoxLayout *rootLayout;
    QVBoxLayout *leftPanel;
    QGridLayout *window_length_settings;
    QLabel *window_length_label;
    QRadioButton *rb_10s;
    QRadioButton *rb_30s;
    QRadioButton *rb_1m;
    QRadioButton *rb_5m;
    QRadioButton *rb_10m;
    QHBoxLayout *eight_hour_nav;
    QPushButton *prev8hours;
    QPushButton *next8hours;
    QHBoxLayout *undo_clearall;
    QPushButton *undo_button;
    QPushButton *clearall_button;
    QHBoxLayout *skip_and_save;
    QPushButton *skip_button;
    QPushButton *finalize_button;
    QHBoxLayout *ecg_marks;
    QPushButton *startNoiseECG;
    QPushButton *stopNoiseECG;
    QHBoxLayout *ppg_marks_row;
    QPushButton *startNoisePPG;
    QPushButton *stopNoisePPG;
    QHBoxLayout *typeRow;
    QLabel *marking_type_label;
    QComboBox *marking_type;
    QHBoxLayout *scrollRow;
    QLabel *seconds_to_scroll_label;
    QLineEdit *skip_interval_box;
    QLabel *marking_type_label_2;
    QGridLayout *marking_type_legend;
    QHBoxLayout *hboxLayout;
    QLabel *signif_color;
    QLabel *signif_desc;
    QHBoxLayout *hboxLayout1;
    QLabel *benign_color;
    QLabel *benign_desc;
    QHBoxLayout *hboxLayout2;
    QLabel *pvc_color;
    QLabel *pvc_desc;
    QHBoxLayout *hboxLayout3;
    QLabel *noise_color;
    QLabel *noise_desc;
    QHBoxLayout *hboxLayout4;
    QLabel *delay_color;
    QLabel *delay_desc;
    QHBoxLayout *hboxLayout5;
    QLabel *vt_color;
    QLabel *vt_desc;
    QHBoxLayout *hboxLayout6;
    QLabel *pac_color;
    QLabel *pac_desc;
    QHBoxLayout *hboxLayout7;
    QLabel *svt_color;
    QLabel *svt_desc;
    QHBoxLayout *hboxLayout8;
    QLabel *af_color;
    QLabel *af_desc;
    QLabel *sleep_states_label;
    QGridLayout *sleep_state_legend;
    QHBoxLayout *hboxLayout9;
    QLabel *n3_c;
    QLabel *n3_l;
    QHBoxLayout *hboxLayout10;
    QLabel *n2_c;
    QLabel *n2_l;
    QHBoxLayout *hboxLayout11;
    QLabel *n1_c;
    QLabel *n1_l;
    QHBoxLayout *hboxLayout12;
    QLabel *rem_c;
    QLabel *rem_l;
    QHBoxLayout *hboxLayout13;
    QLabel *wake_c;
    QLabel *wake_l;
    QSpacerItem *vspacer;
    QVBoxLayout *plotArea;
    QLabel *topLabel;
    QVBoxLayout *ampogram_and_sleepstates;
    QChartView *amp_ecg_axis;
    QChartView *amp_ppg_axis;
    QChartView *sleep_state_axis;
    QLabel *label;
    QChartView *ecg_axis;
    QChartView *ppg_axis;

    void setupUi(QDialog *noise_marking_gui)
    {
        if (noise_marking_gui->objectName().isEmpty())
            noise_marking_gui->setObjectName("noise_marking_gui");
        noise_marking_gui->resize(1284, 1200);
        rootLayout = new QHBoxLayout(noise_marking_gui);
        rootLayout->setSpacing(15);
        rootLayout->setObjectName("rootLayout");
        rootLayout->setContentsMargins(10, 10, 10, 10);
        leftPanel = new QVBoxLayout();
        leftPanel->setSpacing(8);
        leftPanel->setObjectName("leftPanel");
        window_length_settings = new QGridLayout();
        window_length_settings->setObjectName("window_length_settings");
        window_length_label = new QLabel(noise_marking_gui);
        window_length_label->setObjectName("window_length_label");

        window_length_settings->addWidget(window_length_label, 0, 0, 1, 2);

        rb_10s = new QRadioButton(noise_marking_gui);
        rb_10s->setObjectName("rb_10s");

        window_length_settings->addWidget(rb_10s, 1, 0, 1, 1);

        rb_30s = new QRadioButton(noise_marking_gui);
        rb_30s->setObjectName("rb_30s");

        window_length_settings->addWidget(rb_30s, 1, 1, 1, 1);

        rb_1m = new QRadioButton(noise_marking_gui);
        rb_1m->setObjectName("rb_1m");

        window_length_settings->addWidget(rb_1m, 2, 0, 1, 1);

        rb_5m = new QRadioButton(noise_marking_gui);
        rb_5m->setObjectName("rb_5m");

        window_length_settings->addWidget(rb_5m, 2, 1, 1, 1);

        rb_10m = new QRadioButton(noise_marking_gui);
        rb_10m->setObjectName("rb_10m");

        window_length_settings->addWidget(rb_10m, 3, 1, 1, 1);


        leftPanel->addLayout(window_length_settings);

        eight_hour_nav = new QHBoxLayout();
        eight_hour_nav->setObjectName("eight_hour_nav");
        prev8hours = new QPushButton(noise_marking_gui);
        prev8hours->setObjectName("prev8hours");

        eight_hour_nav->addWidget(prev8hours);

        next8hours = new QPushButton(noise_marking_gui);
        next8hours->setObjectName("next8hours");

        eight_hour_nav->addWidget(next8hours);


        leftPanel->addLayout(eight_hour_nav);

        undo_clearall = new QHBoxLayout();
        undo_clearall->setObjectName("undo_clearall");
        undo_button = new QPushButton(noise_marking_gui);
        undo_button->setObjectName("undo_button");

        undo_clearall->addWidget(undo_button);

        clearall_button = new QPushButton(noise_marking_gui);
        clearall_button->setObjectName("clearall_button");

        undo_clearall->addWidget(clearall_button);


        leftPanel->addLayout(undo_clearall);

        skip_and_save = new QHBoxLayout();
        skip_and_save->setObjectName("skip_and_save");
        skip_button = new QPushButton(noise_marking_gui);
        skip_button->setObjectName("skip_button");

        skip_and_save->addWidget(skip_button);

        finalize_button = new QPushButton(noise_marking_gui);
        finalize_button->setObjectName("finalize_button");

        skip_and_save->addWidget(finalize_button);


        leftPanel->addLayout(skip_and_save);

        ecg_marks = new QHBoxLayout();
        ecg_marks->setObjectName("ecg_marks");
        startNoiseECG = new QPushButton(noise_marking_gui);
        startNoiseECG->setObjectName("startNoiseECG");

        ecg_marks->addWidget(startNoiseECG);

        stopNoiseECG = new QPushButton(noise_marking_gui);
        stopNoiseECG->setObjectName("stopNoiseECG");

        ecg_marks->addWidget(stopNoiseECG);


        leftPanel->addLayout(ecg_marks);

        ppg_marks_row = new QHBoxLayout();
        ppg_marks_row->setObjectName("ppg_marks_row");
        startNoisePPG = new QPushButton(noise_marking_gui);
        startNoisePPG->setObjectName("startNoisePPG");

        ppg_marks_row->addWidget(startNoisePPG);

        stopNoisePPG = new QPushButton(noise_marking_gui);
        stopNoisePPG->setObjectName("stopNoisePPG");

        ppg_marks_row->addWidget(stopNoisePPG);


        leftPanel->addLayout(ppg_marks_row);

        typeRow = new QHBoxLayout();
        typeRow->setObjectName("typeRow");
        marking_type_label = new QLabel(noise_marking_gui);
        marking_type_label->setObjectName("marking_type_label");

        typeRow->addWidget(marking_type_label);

        marking_type = new QComboBox(noise_marking_gui);
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->addItem(QString());
        marking_type->setObjectName("marking_type");

        typeRow->addWidget(marking_type);


        leftPanel->addLayout(typeRow);

        scrollRow = new QHBoxLayout();
        scrollRow->setObjectName("scrollRow");
        seconds_to_scroll_label = new QLabel(noise_marking_gui);
        seconds_to_scroll_label->setObjectName("seconds_to_scroll_label");

        scrollRow->addWidget(seconds_to_scroll_label);

        skip_interval_box = new QLineEdit(noise_marking_gui);
        skip_interval_box->setObjectName("skip_interval_box");
        skip_interval_box->setMaximumSize(QSize(50, 25));

        scrollRow->addWidget(skip_interval_box);


        leftPanel->addLayout(scrollRow);

        marking_type_label_2 = new QLabel(noise_marking_gui);
        marking_type_label_2->setObjectName("marking_type_label_2");
        QFont font;
        font.setBold(true);
        marking_type_label_2->setFont(font);

        leftPanel->addWidget(marking_type_label_2);

        marking_type_legend = new QGridLayout();
        marking_type_legend->setSpacing(4);
        marking_type_legend->setObjectName("marking_type_legend");
        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName("hboxLayout");
        signif_color = new QLabel(noise_marking_gui);
        signif_color->setObjectName("signif_color");
        signif_color->setMinimumSize(QSize(15, 15));
        signif_color->setMaximumSize(QSize(15, 15));

        hboxLayout->addWidget(signif_color);

        signif_desc = new QLabel(noise_marking_gui);
        signif_desc->setObjectName("signif_desc");

        hboxLayout->addWidget(signif_desc);


        marking_type_legend->addLayout(hboxLayout, 6, 0, 1, 2);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setObjectName("hboxLayout1");
        benign_color = new QLabel(noise_marking_gui);
        benign_color->setObjectName("benign_color");
        benign_color->setMinimumSize(QSize(15, 15));
        benign_color->setMaximumSize(QSize(15, 15));

        hboxLayout1->addWidget(benign_color);

        benign_desc = new QLabel(noise_marking_gui);
        benign_desc->setObjectName("benign_desc");

        hboxLayout1->addWidget(benign_desc);


        marking_type_legend->addLayout(hboxLayout1, 5, 0, 1, 2);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setObjectName("hboxLayout2");
        pvc_color = new QLabel(noise_marking_gui);
        pvc_color->setObjectName("pvc_color");
        pvc_color->setMinimumSize(QSize(15, 15));
        pvc_color->setMaximumSize(QSize(15, 15));

        hboxLayout2->addWidget(pvc_color);

        pvc_desc = new QLabel(noise_marking_gui);
        pvc_desc->setObjectName("pvc_desc");

        hboxLayout2->addWidget(pvc_desc);


        marking_type_legend->addLayout(hboxLayout2, 2, 1, 1, 1);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setObjectName("hboxLayout3");
        noise_color = new QLabel(noise_marking_gui);
        noise_color->setObjectName("noise_color");
        noise_color->setMinimumSize(QSize(15, 15));
        noise_color->setMaximumSize(QSize(15, 15));

        hboxLayout3->addWidget(noise_color);

        noise_desc = new QLabel(noise_marking_gui);
        noise_desc->setObjectName("noise_desc");

        hboxLayout3->addWidget(noise_desc);


        marking_type_legend->addLayout(hboxLayout3, 0, 0, 1, 1);

        hboxLayout4 = new QHBoxLayout();
        hboxLayout4->setObjectName("hboxLayout4");
        delay_color = new QLabel(noise_marking_gui);
        delay_color->setObjectName("delay_color");
        delay_color->setMinimumSize(QSize(15, 15));
        delay_color->setMaximumSize(QSize(15, 15));

        hboxLayout4->addWidget(delay_color);

        delay_desc = new QLabel(noise_marking_gui);
        delay_desc->setObjectName("delay_desc");

        hboxLayout4->addWidget(delay_desc);


        marking_type_legend->addLayout(hboxLayout4, 1, 0, 1, 1);

        hboxLayout5 = new QHBoxLayout();
        hboxLayout5->setObjectName("hboxLayout5");
        vt_color = new QLabel(noise_marking_gui);
        vt_color->setObjectName("vt_color");
        vt_color->setMinimumSize(QSize(15, 15));
        vt_color->setMaximumSize(QSize(15, 15));

        hboxLayout5->addWidget(vt_color);

        vt_desc = new QLabel(noise_marking_gui);
        vt_desc->setObjectName("vt_desc");

        hboxLayout5->addWidget(vt_desc);


        marking_type_legend->addLayout(hboxLayout5, 0, 1, 1, 1);

        hboxLayout6 = new QHBoxLayout();
        hboxLayout6->setObjectName("hboxLayout6");
        pac_color = new QLabel(noise_marking_gui);
        pac_color->setObjectName("pac_color");
        pac_color->setMinimumSize(QSize(15, 15));
        pac_color->setMaximumSize(QSize(15, 15));

        hboxLayout6->addWidget(pac_color);

        pac_desc = new QLabel(noise_marking_gui);
        pac_desc->setObjectName("pac_desc");

        hboxLayout6->addWidget(pac_desc);


        marking_type_legend->addLayout(hboxLayout6, 1, 1, 1, 1);

        hboxLayout7 = new QHBoxLayout();
        hboxLayout7->setObjectName("hboxLayout7");
        svt_color = new QLabel(noise_marking_gui);
        svt_color->setObjectName("svt_color");
        svt_color->setMinimumSize(QSize(15, 15));
        svt_color->setMaximumSize(QSize(15, 15));

        hboxLayout7->addWidget(svt_color);

        svt_desc = new QLabel(noise_marking_gui);
        svt_desc->setObjectName("svt_desc");

        hboxLayout7->addWidget(svt_desc);

        hboxLayout8 = new QHBoxLayout();
        hboxLayout8->setObjectName("hboxLayout8");
        af_color = new QLabel(noise_marking_gui);
        af_color->setObjectName("af_color");
        af_color->setMinimumSize(QSize(15, 15));
        af_color->setMaximumSize(QSize(15, 15));

        hboxLayout8->addWidget(af_color);

        af_desc = new QLabel(noise_marking_gui);
        af_desc->setObjectName("af_desc");

        hboxLayout8->addWidget(af_desc);


        hboxLayout7->addLayout(hboxLayout8);


        marking_type_legend->addLayout(hboxLayout7, 2, 0, 1, 1);


        leftPanel->addLayout(marking_type_legend);

        sleep_states_label = new QLabel(noise_marking_gui);
        sleep_states_label->setObjectName("sleep_states_label");
        sleep_states_label->setFont(font);

        leftPanel->addWidget(sleep_states_label);

        sleep_state_legend = new QGridLayout();
        sleep_state_legend->setSpacing(4);
        sleep_state_legend->setObjectName("sleep_state_legend");
        hboxLayout9 = new QHBoxLayout();
        hboxLayout9->setObjectName("hboxLayout9");
        n3_c = new QLabel(noise_marking_gui);
        n3_c->setObjectName("n3_c");
        n3_c->setMinimumSize(QSize(15, 15));
        n3_c->setMaximumSize(QSize(15, 15));

        hboxLayout9->addWidget(n3_c);

        n3_l = new QLabel(noise_marking_gui);
        n3_l->setObjectName("n3_l");

        hboxLayout9->addWidget(n3_l);


        sleep_state_legend->addLayout(hboxLayout9, 2, 0, 1, 1);

        hboxLayout10 = new QHBoxLayout();
        hboxLayout10->setObjectName("hboxLayout10");
        n2_c = new QLabel(noise_marking_gui);
        n2_c->setObjectName("n2_c");
        n2_c->setMinimumSize(QSize(15, 15));
        n2_c->setMaximumSize(QSize(15, 15));

        hboxLayout10->addWidget(n2_c);

        n2_l = new QLabel(noise_marking_gui);
        n2_l->setObjectName("n2_l");

        hboxLayout10->addWidget(n2_l);


        sleep_state_legend->addLayout(hboxLayout10, 1, 1, 1, 1);

        hboxLayout11 = new QHBoxLayout();
        hboxLayout11->setObjectName("hboxLayout11");
        n1_c = new QLabel(noise_marking_gui);
        n1_c->setObjectName("n1_c");
        n1_c->setMinimumSize(QSize(15, 15));
        n1_c->setMaximumSize(QSize(15, 15));

        hboxLayout11->addWidget(n1_c);

        n1_l = new QLabel(noise_marking_gui);
        n1_l->setObjectName("n1_l");

        hboxLayout11->addWidget(n1_l);


        sleep_state_legend->addLayout(hboxLayout11, 1, 0, 1, 1);

        hboxLayout12 = new QHBoxLayout();
        hboxLayout12->setObjectName("hboxLayout12");
        rem_c = new QLabel(noise_marking_gui);
        rem_c->setObjectName("rem_c");
        rem_c->setMinimumSize(QSize(15, 15));
        rem_c->setMaximumSize(QSize(15, 15));

        hboxLayout12->addWidget(rem_c);

        rem_l = new QLabel(noise_marking_gui);
        rem_l->setObjectName("rem_l");

        hboxLayout12->addWidget(rem_l);


        sleep_state_legend->addLayout(hboxLayout12, 0, 1, 1, 1);

        hboxLayout13 = new QHBoxLayout();
        hboxLayout13->setObjectName("hboxLayout13");
        wake_c = new QLabel(noise_marking_gui);
        wake_c->setObjectName("wake_c");
        wake_c->setMinimumSize(QSize(15, 15));
        wake_c->setMaximumSize(QSize(15, 15));

        hboxLayout13->addWidget(wake_c);

        wake_l = new QLabel(noise_marking_gui);
        wake_l->setObjectName("wake_l");

        hboxLayout13->addWidget(wake_l);


        sleep_state_legend->addLayout(hboxLayout13, 0, 0, 1, 1);


        leftPanel->addLayout(sleep_state_legend);

        vspacer = new QSpacerItem(0, 0, QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Expanding);

        leftPanel->addItem(vspacer);


        rootLayout->addLayout(leftPanel);

        plotArea = new QVBoxLayout();
        plotArea->setSpacing(0);
        plotArea->setObjectName("plotArea");
        topLabel = new QLabel(noise_marking_gui);
        topLabel->setObjectName("topLabel");

        plotArea->addWidget(topLabel);

        ampogram_and_sleepstates = new QVBoxLayout();
        ampogram_and_sleepstates->setSpacing(0);
        ampogram_and_sleepstates->setObjectName("ampogram_and_sleepstates");
        amp_ecg_axis = new QChartView(noise_marking_gui);
        amp_ecg_axis->setObjectName("amp_ecg_axis");

        ampogram_and_sleepstates->addWidget(amp_ecg_axis);

        amp_ppg_axis = new QChartView(noise_marking_gui);
        amp_ppg_axis->setObjectName("amp_ppg_axis");

        ampogram_and_sleepstates->addWidget(amp_ppg_axis);

        sleep_state_axis = new QChartView(noise_marking_gui);
        sleep_state_axis->setObjectName("sleep_state_axis");
        sleep_state_axis->setMinimumSize(QSize(0, 60));
        label = new QLabel(sleep_state_axis);
        label->setObjectName("label");
        label->setGeometry(QRect(10, 40, 49, 16));
        label->setStyleSheet(QString::fromUtf8("border:none;"));

        ampogram_and_sleepstates->addWidget(sleep_state_axis);


        plotArea->addLayout(ampogram_and_sleepstates);

        ecg_axis = new QChartView(noise_marking_gui);
        ecg_axis->setObjectName("ecg_axis");

        plotArea->addWidget(ecg_axis);

        ppg_axis = new QChartView(noise_marking_gui);
        ppg_axis->setObjectName("ppg_axis");

        plotArea->addWidget(ppg_axis);

        plotArea->setStretch(1, 1);
        plotArea->setStretch(2, 2);
        plotArea->setStretch(3, 2);

        rootLayout->addLayout(plotArea);

        rootLayout->setStretch(1, 1);

        retranslateUi(noise_marking_gui);

        QMetaObject::connectSlotsByName(noise_marking_gui);
    } // setupUi

    void retranslateUi(QDialog *noise_marking_gui)
    {
        noise_marking_gui->setWindowTitle(QCoreApplication::translate("noise_marking_gui", "PPG/ECG Noise Marker", nullptr));
        window_length_label->setText(QCoreApplication::translate("noise_marking_gui", "Window Length", nullptr));
        rb_10s->setText(QCoreApplication::translate("noise_marking_gui", "10s", nullptr));
        rb_30s->setText(QCoreApplication::translate("noise_marking_gui", "30s", nullptr));
        rb_1m->setText(QCoreApplication::translate("noise_marking_gui", "1m", nullptr));
        rb_5m->setText(QCoreApplication::translate("noise_marking_gui", "5m", nullptr));
        rb_10m->setText(QCoreApplication::translate("noise_marking_gui", "10m", nullptr));
        prev8hours->setText(QCoreApplication::translate("noise_marking_gui", "Prev 8h", nullptr));
        next8hours->setText(QCoreApplication::translate("noise_marking_gui", "Next 8h", nullptr));
        undo_button->setText(QCoreApplication::translate("noise_marking_gui", "Undo", nullptr));
        clearall_button->setText(QCoreApplication::translate("noise_marking_gui", "Clear All", nullptr));
        skip_button->setText(QCoreApplication::translate("noise_marking_gui", "Skip", nullptr));
        finalize_button->setText(QCoreApplication::translate("noise_marking_gui", "Save", nullptr));
        startNoiseECG->setText(QCoreApplication::translate("noise_marking_gui", "Mark ECG Start", nullptr));
        stopNoiseECG->setText(QCoreApplication::translate("noise_marking_gui", "Mark ECG End", nullptr));
        startNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG Start", nullptr));
        stopNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG End", nullptr));
        marking_type_label->setText(QCoreApplication::translate("noise_marking_gui", "Mark Type", nullptr));
        marking_type->setItemText(0, QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        marking_type->setItemText(1, QCoreApplication::translate("noise_marking_gui", "Conduction Delay", nullptr));
        marking_type->setItemText(2, QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        marking_type->setItemText(3, QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        marking_type->setItemText(4, QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        marking_type->setItemText(5, QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        marking_type->setItemText(6, QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        marking_type->setItemText(7, QCoreApplication::translate("noise_marking_gui", "Benign Arrhythmia", nullptr));
        marking_type->setItemText(8, QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));

        seconds_to_scroll_label->setText(QCoreApplication::translate("noise_marking_gui", "Scroll Dist", nullptr));
        skip_interval_box->setText(QCoreApplication::translate("noise_marking_gui", "5.0", nullptr));
        marking_type_label_2->setText(QCoreApplication::translate("noise_marking_gui", "Marking Type Legend", nullptr));
        signif_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(0, 255, 255); border-radius:3px;", nullptr));
        signif_desc->setText(QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));
        benign_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(255, 128, 255); border-radius:3px;", nullptr));
        benign_desc->setText(QCoreApplication::translate("noise_marking_gui", "Benign Arrhythmia", nullptr));
        pvc_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(128, 255, 0); border-radius:3px;", nullptr));
        pvc_desc->setText(QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        noise_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: yellow;", nullptr));
        noise_desc->setText(QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        delay_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(128,0,128); border-radius:3px;", nullptr));
        delay_desc->setText(QCoreApplication::translate("noise_marking_gui", "Conduction Delay", nullptr));
        vt_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: green; border-radius:3px;", nullptr));
        vt_desc->setText(QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        pac_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(255,128, 0); border-radius:3px;", nullptr));
        pac_desc->setText(QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        svt_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: blue; border-radius:3px;", nullptr));
        svt_desc->setText(QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        af_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: red; border-radius:3px;", nullptr));
        af_desc->setText(QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        sleep_states_label->setText(QCoreApplication::translate("noise_marking_gui", "Sleep States Legend", nullptr));
        n3_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: red; border-radius:3px;", nullptr));
        n3_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM3", nullptr));
        n2_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: cyan; border-radius:3px;", nullptr));
        n2_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM2", nullptr));
        n1_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: blue; border-radius:3px;", nullptr));
        n1_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM1", nullptr));
        rem_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: green; border-radius:3px;", nullptr));
        rem_l->setText(QCoreApplication::translate("noise_marking_gui", "REM", nullptr));
        wake_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: black; border-radius:3px;", nullptr));
        wake_l->setText(QCoreApplication::translate("noise_marking_gui", "Wake", nullptr));
        topLabel->setText(QString());
        amp_ecg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        amp_ppg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        sleep_state_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        label->setText(QCoreApplication::translate("noise_marking_gui", "Time (h)", nullptr));
        ecg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 2px solid black;", nullptr));
        ppg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid black;", nullptr));
    } // retranslateUi

};

namespace Ui {
    class noise_marking_gui: public Ui_noise_marking_gui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_NOISE_MARKING_GUI_H
