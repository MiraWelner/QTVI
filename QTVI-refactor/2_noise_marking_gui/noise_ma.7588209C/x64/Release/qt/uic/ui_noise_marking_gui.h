/********************************************************************************
** Form generated from reading UI file 'noise_marking_gui.ui'
**
** Created by: Qt User Interface Compiler version 6.10.1
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
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_noise_marking_gui
{
public:
    QWidget *layoutWidget;
    QVBoxLayout *verticalLayout;
    QLabel *topLabel;
    QHBoxLayout *top_row;
    QGroupBox *legend;
    QVBoxLayout *leftBoxLayout;
    QHBoxLayout *nrem1_label;
    QLabel *label_9;
    QLabel *label_5;
    QHBoxLayout *nrem2_label;
    QLabel *label_8;
    QLabel *label_4;
    QHBoxLayout *nrem3_label;
    QLabel *label_7;
    QLabel *label_3;
    QHBoxLayout *rem_label;
    QLabel *label_6;
    QLabel *label_2;
    QHBoxLayout *wake_label;
    QLabel *wake_color;
    QLabel *wake_text;
    QVBoxLayout *amp_and_sleep_axes;
    QChartView *amp_ecg_axis;
    QChartView *amp_ppg_axis;
    QChartView *sleep_state_axis;
    QLabel *hours_label;
    QGroupBox *selection_type_legend;
    QVBoxLayout *verticalLayout_2;
    QHBoxLayout *noise_label_box;
    QLabel *noise_color;
    QLabel *noise_desc;
    QHBoxLayout *vt_label_box_2;
    QLabel *label_19;
    QLabel *label_20;
    QHBoxLayout *vt_label_box;
    QLabel *label_17;
    QLabel *label_18;
    QHBoxLayout *af_label_box;
    QLabel *af_color;
    QLabel *af_desc;
    QHBoxLayout *svt_label_box;
    QLabel *label_15;
    QLabel *label_16;
    QHBoxLayout *pvc_label_box;
    QLabel *wake_color_2;
    QLabel *wake_text_2;
    QHBoxLayout *pac_label_box;
    QLabel *pac_color;
    QLabel *pac_desc;
    QHBoxLayout *benign_label_box;
    QLabel *wake_color_5;
    QLabel *wake_text_5;
    QHBoxLayout *significant_arrhythmia;
    QLabel *wake_color_6;
    QLabel *wake_text_6;
    QVBoxLayout *data_plots;
    QChartView *ecg_axis;
    QChartView *ppg_axis;
    QHBoxLayout *horizontalLayout;
    QGroupBox *window_length;
    QGridLayout *gridLayout_buttons;
    QRadioButton *rb_10s;
    QRadioButton *rb_30s;
    QRadioButton *rb_1m;
    QRadioButton *rb_5m;
    QRadioButton *rb_10m;
    QVBoxLayout *clear_undo;
    QPushButton *undo_button;
    QPushButton *clearall_button;
    QVBoxLayout *save_skip;
    QPushButton *skip_button;
    QPushButton *finalize_button;
    QVBoxLayout *ecg_button_layout;
    QPushButton *startNoiseECG;
    QPushButton *stopNoiseECG;
    QVBoxLayout *ppg_startstop;
    QPushButton *startNoisePPG;
    QPushButton *stopNoisePPG;
    QVBoxLayout *marking_type_box;
    QLabel *marking_type_label;
    QComboBox *marking_type;
    QVBoxLayout *scroll_seconds_layout;
    QLabel *seconds_to_scroll_label;
    QLineEdit *skip_interval_box;
    QVBoxLayout *ampogram_shift;
    QPushButton *prev8hours;
    QPushButton *next8hours;
    QSpacerItem *horizontalSpacer_2;

    void setupUi(QDialog *noise_marking_gui)
    {
        if (noise_marking_gui->objectName().isEmpty())
            noise_marking_gui->setObjectName("noise_marking_gui");
        noise_marking_gui->resize(1026, 800);
        QSizePolicy sizePolicy(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(noise_marking_gui->sizePolicy().hasHeightForWidth());
        noise_marking_gui->setSizePolicy(sizePolicy);
        layoutWidget = new QWidget(noise_marking_gui);
        layoutWidget->setObjectName("layoutWidget");
        layoutWidget->setGeometry(QRect(0, 0, 1028, 801));
        verticalLayout = new QVBoxLayout(layoutWidget);
        verticalLayout->setObjectName("verticalLayout");
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        topLabel = new QLabel(layoutWidget);
        topLabel->setObjectName("topLabel");

        verticalLayout->addWidget(topLabel);

        top_row = new QHBoxLayout();
        top_row->setSpacing(0);
        top_row->setObjectName("top_row");
        legend = new QGroupBox(layoutWidget);
        legend->setObjectName("legend");
        QSizePolicy sizePolicy1(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(0);
        sizePolicy1.setHeightForWidth(legend->sizePolicy().hasHeightForWidth());
        legend->setSizePolicy(sizePolicy1);
        legend->setMinimumSize(QSize(35, 0));
        leftBoxLayout = new QVBoxLayout(legend);
        leftBoxLayout->setSpacing(1);
        leftBoxLayout->setObjectName("leftBoxLayout");
        leftBoxLayout->setContentsMargins(5, 2, 5, 2);
        nrem1_label = new QHBoxLayout();
        nrem1_label->setObjectName("nrem1_label");
        label_9 = new QLabel(legend);
        label_9->setObjectName("label_9");
        QSizePolicy sizePolicy2(QSizePolicy::Policy::Fixed, QSizePolicy::Policy::Fixed);
        sizePolicy2.setHorizontalStretch(0);
        sizePolicy2.setVerticalStretch(0);
        sizePolicy2.setHeightForWidth(label_9->sizePolicy().hasHeightForWidth());
        label_9->setSizePolicy(sizePolicy2);
        label_9->setMinimumSize(QSize(15, 15));
        label_9->setMaximumSize(QSize(15, 15));
        label_9->setStyleSheet(QString::fromUtf8("background-color: black;\n"
"border-radius:3px;"));

        nrem1_label->addWidget(label_9);

        label_5 = new QLabel(legend);
        label_5->setObjectName("label_5");
        QSizePolicy sizePolicy3(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Ignored);
        sizePolicy3.setHorizontalStretch(0);
        sizePolicy3.setVerticalStretch(0);
        sizePolicy3.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy3);

        nrem1_label->addWidget(label_5);


        leftBoxLayout->addLayout(nrem1_label);

        nrem2_label = new QHBoxLayout();
        nrem2_label->setObjectName("nrem2_label");
        label_8 = new QLabel(legend);
        label_8->setObjectName("label_8");
        QSizePolicy sizePolicy4(QSizePolicy::Policy::Fixed, QSizePolicy::Policy::Ignored);
        sizePolicy4.setHorizontalStretch(0);
        sizePolicy4.setVerticalStretch(0);
        sizePolicy4.setHeightForWidth(label_8->sizePolicy().hasHeightForWidth());
        label_8->setSizePolicy(sizePolicy4);
        label_8->setMinimumSize(QSize(15, 15));
        label_8->setMaximumSize(QSize(15, 15));
        label_8->setStyleSheet(QString::fromUtf8("background-color: green;\n"
"border-radius: 3px;"));

        nrem2_label->addWidget(label_8);

        label_4 = new QLabel(legend);
        label_4->setObjectName("label_4");

        nrem2_label->addWidget(label_4);


        leftBoxLayout->addLayout(nrem2_label);

        nrem3_label = new QHBoxLayout();
        nrem3_label->setObjectName("nrem3_label");
        label_7 = new QLabel(legend);
        label_7->setObjectName("label_7");
        sizePolicy2.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy2);
        label_7->setMinimumSize(QSize(15, 15));
        label_7->setMaximumSize(QSize(15, 15));
        label_7->setStyleSheet(QString::fromUtf8("background-color: blue;\n"
"border-radius: 3px;"));

        nrem3_label->addWidget(label_7);

        label_3 = new QLabel(legend);
        label_3->setObjectName("label_3");

        nrem3_label->addWidget(label_3);


        leftBoxLayout->addLayout(nrem3_label);

        rem_label = new QHBoxLayout();
        rem_label->setObjectName("rem_label");
        label_6 = new QLabel(legend);
        label_6->setObjectName("label_6");
        sizePolicy2.setHeightForWidth(label_6->sizePolicy().hasHeightForWidth());
        label_6->setSizePolicy(sizePolicy2);
        label_6->setMinimumSize(QSize(15, 15));
        label_6->setMaximumSize(QSize(15, 15));
        label_6->setStyleSheet(QString::fromUtf8("background-color: cyan;\n"
"border-radius: 3px;"));

        rem_label->addWidget(label_6);

        label_2 = new QLabel(legend);
        label_2->setObjectName("label_2");

        rem_label->addWidget(label_2);


        leftBoxLayout->addLayout(rem_label);

        wake_label = new QHBoxLayout();
        wake_label->setObjectName("wake_label");
        wake_color = new QLabel(legend);
        wake_color->setObjectName("wake_color");
        sizePolicy2.setHeightForWidth(wake_color->sizePolicy().hasHeightForWidth());
        wake_color->setSizePolicy(sizePolicy2);
        wake_color->setMinimumSize(QSize(15, 15));
        wake_color->setMaximumSize(QSize(15, 15));
        wake_color->setStyleSheet(QString::fromUtf8("background-color: red;\n"
"border-radius: 3px;"));

        wake_label->addWidget(wake_color);

        wake_text = new QLabel(legend);
        wake_text->setObjectName("wake_text");

        wake_label->addWidget(wake_text);

        wake_label->setStretch(1, 1);

        leftBoxLayout->addLayout(wake_label);


        top_row->addWidget(legend);

        amp_and_sleep_axes = new QVBoxLayout();
        amp_and_sleep_axes->setSpacing(0);
        amp_and_sleep_axes->setObjectName("amp_and_sleep_axes");
        amp_ecg_axis = new QChartView(layoutWidget);
        amp_ecg_axis->setObjectName("amp_ecg_axis");
        QSizePolicy sizePolicy5(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Preferred);
        sizePolicy5.setHorizontalStretch(0);
        sizePolicy5.setVerticalStretch(0);
        sizePolicy5.setHeightForWidth(amp_ecg_axis->sizePolicy().hasHeightForWidth());
        amp_ecg_axis->setSizePolicy(sizePolicy5);
        amp_ecg_axis->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255); border: 1px solid black;"));

        amp_and_sleep_axes->addWidget(amp_ecg_axis);

        amp_ppg_axis = new QChartView(layoutWidget);
        amp_ppg_axis->setObjectName("amp_ppg_axis");
        sizePolicy5.setHeightForWidth(amp_ppg_axis->sizePolicy().hasHeightForWidth());
        amp_ppg_axis->setSizePolicy(sizePolicy5);
        amp_ppg_axis->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255); border: 1px solid black;"));

        amp_and_sleep_axes->addWidget(amp_ppg_axis);

        sleep_state_axis = new QChartView(layoutWidget);
        sleep_state_axis->setObjectName("sleep_state_axis");
        sizePolicy5.setHeightForWidth(sleep_state_axis->sizePolicy().hasHeightForWidth());
        sleep_state_axis->setSizePolicy(sizePolicy5);
        sleep_state_axis->setMinimumSize(QSize(0, 70));
        sleep_state_axis->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255); border: 1px solid black;"));
        hours_label = new QLabel(sleep_state_axis);
        hours_label->setObjectName("hours_label");
        hours_label->setGeometry(QRect(10, 50, 71, 20));
        hours_label->setStyleSheet(QString::fromUtf8("border:none;"));

        amp_and_sleep_axes->addWidget(sleep_state_axis);


        top_row->addLayout(amp_and_sleep_axes);

        selection_type_legend = new QGroupBox(layoutWidget);
        selection_type_legend->setObjectName("selection_type_legend");
        QSizePolicy sizePolicy6(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Expanding);
        sizePolicy6.setHorizontalStretch(0);
        sizePolicy6.setVerticalStretch(0);
        sizePolicy6.setHeightForWidth(selection_type_legend->sizePolicy().hasHeightForWidth());
        selection_type_legend->setSizePolicy(sizePolicy6);
        selection_type_legend->setMinimumSize(QSize(0, 0));
        verticalLayout_2 = new QVBoxLayout(selection_type_legend);
        verticalLayout_2->setSpacing(2);
        verticalLayout_2->setObjectName("verticalLayout_2");
        verticalLayout_2->setContentsMargins(5, 2, 5, 2);
        noise_label_box = new QHBoxLayout();
        noise_label_box->setSpacing(2);
        noise_label_box->setObjectName("noise_label_box");
        noise_color = new QLabel(selection_type_legend);
        noise_color->setObjectName("noise_color");
        sizePolicy2.setHeightForWidth(noise_color->sizePolicy().hasHeightForWidth());
        noise_color->setSizePolicy(sizePolicy2);
        noise_color->setMinimumSize(QSize(15, 15));
        noise_color->setMaximumSize(QSize(15, 15));
        noise_color->setStyleSheet(QString::fromUtf8("background-color: yellow"));

        noise_label_box->addWidget(noise_color);

        noise_desc = new QLabel(selection_type_legend);
        noise_desc->setObjectName("noise_desc");

        noise_label_box->addWidget(noise_desc);

        noise_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(noise_label_box);

        vt_label_box_2 = new QHBoxLayout();
        vt_label_box_2->setObjectName("vt_label_box_2");
        label_19 = new QLabel(selection_type_legend);
        label_19->setObjectName("label_19");
        sizePolicy2.setHeightForWidth(label_19->sizePolicy().hasHeightForWidth());
        label_19->setSizePolicy(sizePolicy2);
        label_19->setMinimumSize(QSize(15, 15));
        label_19->setMaximumSize(QSize(15, 15));
        label_19->setStyleSheet(QString::fromUtf8("background-color: rgb(128,0,128);\n"
"border-radius: 3px;"));

        vt_label_box_2->addWidget(label_19);

        label_20 = new QLabel(selection_type_legend);
        label_20->setObjectName("label_20");

        vt_label_box_2->addWidget(label_20);


        verticalLayout_2->addLayout(vt_label_box_2);

        vt_label_box = new QHBoxLayout();
        vt_label_box->setSpacing(2);
        vt_label_box->setObjectName("vt_label_box");
        label_17 = new QLabel(selection_type_legend);
        label_17->setObjectName("label_17");
        sizePolicy2.setHeightForWidth(label_17->sizePolicy().hasHeightForWidth());
        label_17->setSizePolicy(sizePolicy2);
        label_17->setMinimumSize(QSize(15, 15));
        label_17->setMaximumSize(QSize(15, 15));
        label_17->setStyleSheet(QString::fromUtf8("background-color: green;\n"
"border-radius: 3px;"));

        vt_label_box->addWidget(label_17);

        label_18 = new QLabel(selection_type_legend);
        label_18->setObjectName("label_18");

        vt_label_box->addWidget(label_18);

        vt_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(vt_label_box);

        af_label_box = new QHBoxLayout();
        af_label_box->setSpacing(2);
        af_label_box->setObjectName("af_label_box");
        af_color = new QLabel(selection_type_legend);
        af_color->setObjectName("af_color");
        sizePolicy4.setHeightForWidth(af_color->sizePolicy().hasHeightForWidth());
        af_color->setSizePolicy(sizePolicy4);
        af_color->setMinimumSize(QSize(15, 15));
        af_color->setMaximumSize(QSize(15, 15));
        af_color->setStyleSheet(QString::fromUtf8("background-color: red;\n"
"border-radius: 3px;"));

        af_label_box->addWidget(af_color);

        af_desc = new QLabel(selection_type_legend);
        af_desc->setObjectName("af_desc");

        af_label_box->addWidget(af_desc);

        af_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(af_label_box);

        svt_label_box = new QHBoxLayout();
        svt_label_box->setSpacing(2);
        svt_label_box->setObjectName("svt_label_box");
        label_15 = new QLabel(selection_type_legend);
        label_15->setObjectName("label_15");
        sizePolicy2.setHeightForWidth(label_15->sizePolicy().hasHeightForWidth());
        label_15->setSizePolicy(sizePolicy2);
        label_15->setMinimumSize(QSize(15, 15));
        label_15->setMaximumSize(QSize(15, 15));
        label_15->setStyleSheet(QString::fromUtf8("background-color: blue;\n"
"border-radius: 3px;"));

        svt_label_box->addWidget(label_15);

        label_16 = new QLabel(selection_type_legend);
        label_16->setObjectName("label_16");

        svt_label_box->addWidget(label_16);

        svt_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(svt_label_box);

        pvc_label_box = new QHBoxLayout();
        pvc_label_box->setSpacing(2);
        pvc_label_box->setObjectName("pvc_label_box");
        wake_color_2 = new QLabel(selection_type_legend);
        wake_color_2->setObjectName("wake_color_2");
        sizePolicy2.setHeightForWidth(wake_color_2->sizePolicy().hasHeightForWidth());
        wake_color_2->setSizePolicy(sizePolicy2);
        wake_color_2->setMinimumSize(QSize(15, 15));
        wake_color_2->setMaximumSize(QSize(15, 15));
        wake_color_2->setStyleSheet(QString::fromUtf8("background-color: rgb(128, 255, 0);\n"
"border-radius: 3px;"));

        pvc_label_box->addWidget(wake_color_2);

        wake_text_2 = new QLabel(selection_type_legend);
        wake_text_2->setObjectName("wake_text_2");

        pvc_label_box->addWidget(wake_text_2);

        pvc_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(pvc_label_box);

        pac_label_box = new QHBoxLayout();
        pac_label_box->setSpacing(2);
        pac_label_box->setObjectName("pac_label_box");
        pac_color = new QLabel(selection_type_legend);
        pac_color->setObjectName("pac_color");
        sizePolicy2.setHeightForWidth(pac_color->sizePolicy().hasHeightForWidth());
        pac_color->setSizePolicy(sizePolicy2);
        pac_color->setMinimumSize(QSize(15, 15));
        pac_color->setMaximumSize(QSize(15, 15));
        pac_color->setStyleSheet(QString::fromUtf8("background-color: rgb(255,128, 0);\n"
"border-radius: 3px;"));

        pac_label_box->addWidget(pac_color);

        pac_desc = new QLabel(selection_type_legend);
        pac_desc->setObjectName("pac_desc");

        pac_label_box->addWidget(pac_desc);

        pac_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(pac_label_box);

        benign_label_box = new QHBoxLayout();
        benign_label_box->setSpacing(2);
        benign_label_box->setObjectName("benign_label_box");
        wake_color_5 = new QLabel(selection_type_legend);
        wake_color_5->setObjectName("wake_color_5");
        sizePolicy2.setHeightForWidth(wake_color_5->sizePolicy().hasHeightForWidth());
        wake_color_5->setSizePolicy(sizePolicy2);
        wake_color_5->setMinimumSize(QSize(15, 15));
        wake_color_5->setMaximumSize(QSize(15, 15));
        wake_color_5->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 128, 255);\n"
"border-radius: 3px;"));

        benign_label_box->addWidget(wake_color_5);

        wake_text_5 = new QLabel(selection_type_legend);
        wake_text_5->setObjectName("wake_text_5");

        benign_label_box->addWidget(wake_text_5);

        benign_label_box->setStretch(1, 1);

        verticalLayout_2->addLayout(benign_label_box);

        significant_arrhythmia = new QHBoxLayout();
        significant_arrhythmia->setSpacing(2);
        significant_arrhythmia->setObjectName("significant_arrhythmia");
        wake_color_6 = new QLabel(selection_type_legend);
        wake_color_6->setObjectName("wake_color_6");
        sizePolicy2.setHeightForWidth(wake_color_6->sizePolicy().hasHeightForWidth());
        wake_color_6->setSizePolicy(sizePolicy2);
        wake_color_6->setMinimumSize(QSize(15, 15));
        wake_color_6->setMaximumSize(QSize(15, 15));
        wake_color_6->setStyleSheet(QString::fromUtf8("background-color: rgb(0, 255, 255);\n"
"border-radius: 3px;"));

        significant_arrhythmia->addWidget(wake_color_6);

        wake_text_6 = new QLabel(selection_type_legend);
        wake_text_6->setObjectName("wake_text_6");

        significant_arrhythmia->addWidget(wake_text_6);

        significant_arrhythmia->setStretch(1, 1);

        verticalLayout_2->addLayout(significant_arrhythmia);


        top_row->addWidget(selection_type_legend);

        top_row->setStretch(1, 1);

        verticalLayout->addLayout(top_row);

        data_plots = new QVBoxLayout();
        data_plots->setSpacing(0);
        data_plots->setObjectName("data_plots");
        data_plots->setContentsMargins(-1, -1, 0, -1);
        ecg_axis = new QChartView(layoutWidget);
        ecg_axis->setObjectName("ecg_axis");
        sizePolicy.setHeightForWidth(ecg_axis->sizePolicy().hasHeightForWidth());
        ecg_axis->setSizePolicy(sizePolicy);
        ecg_axis->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255); border: 1px solid black;"));

        data_plots->addWidget(ecg_axis);

        ppg_axis = new QChartView(layoutWidget);
        ppg_axis->setObjectName("ppg_axis");
        QSizePolicy sizePolicy7(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);
        sizePolicy7.setHorizontalStretch(0);
        sizePolicy7.setVerticalStretch(0);
        sizePolicy7.setHeightForWidth(ppg_axis->sizePolicy().hasHeightForWidth());
        ppg_axis->setSizePolicy(sizePolicy7);
        ppg_axis->setStyleSheet(QString::fromUtf8("background-color: rgb(255, 255, 255); border: 1px solid black;"));

        data_plots->addWidget(ppg_axis);

        data_plots->setStretch(0, 1);
        data_plots->setStretch(1, 1);

        verticalLayout->addLayout(data_plots);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName("horizontalLayout");
        window_length = new QGroupBox(layoutWidget);
        window_length->setObjectName("window_length");
        QSizePolicy sizePolicy8(QSizePolicy::Policy::Minimum, QSizePolicy::Policy::Minimum);
        sizePolicy8.setHorizontalStretch(0);
        sizePolicy8.setVerticalStretch(0);
        sizePolicy8.setHeightForWidth(window_length->sizePolicy().hasHeightForWidth());
        window_length->setSizePolicy(sizePolicy8);
        window_length->setMaximumSize(QSize(200, 70));
        gridLayout_buttons = new QGridLayout(window_length);
        gridLayout_buttons->setSpacing(1);
        gridLayout_buttons->setObjectName("gridLayout_buttons");
        gridLayout_buttons->setContentsMargins(5, 5, 5, 5);
        rb_10s = new QRadioButton(window_length);
        rb_10s->setObjectName("rb_10s");

        gridLayout_buttons->addWidget(rb_10s, 0, 0, 1, 1);

        rb_30s = new QRadioButton(window_length);
        rb_30s->setObjectName("rb_30s");

        gridLayout_buttons->addWidget(rb_30s, 0, 1, 1, 1);

        rb_1m = new QRadioButton(window_length);
        rb_1m->setObjectName("rb_1m");

        gridLayout_buttons->addWidget(rb_1m, 0, 2, 1, 1);

        rb_5m = new QRadioButton(window_length);
        rb_5m->setObjectName("rb_5m");

        gridLayout_buttons->addWidget(rb_5m, 1, 0, 1, 1);

        rb_10m = new QRadioButton(window_length);
        rb_10m->setObjectName("rb_10m");

        gridLayout_buttons->addWidget(rb_10m, 1, 1, 1, 1);


        horizontalLayout->addWidget(window_length);

        clear_undo = new QVBoxLayout();
        clear_undo->setSpacing(0);
        clear_undo->setObjectName("clear_undo");
        undo_button = new QPushButton(layoutWidget);
        undo_button->setObjectName("undo_button");
        undo_button->setMaximumSize(QSize(100, 16777215));

        clear_undo->addWidget(undo_button);

        clearall_button = new QPushButton(layoutWidget);
        clearall_button->setObjectName("clearall_button");
        clearall_button->setMaximumSize(QSize(100, 16777215));

        clear_undo->addWidget(clearall_button);


        horizontalLayout->addLayout(clear_undo);

        save_skip = new QVBoxLayout();
        save_skip->setObjectName("save_skip");
        skip_button = new QPushButton(layoutWidget);
        skip_button->setObjectName("skip_button");
        skip_button->setMaximumSize(QSize(100, 16777215));

        save_skip->addWidget(skip_button);

        finalize_button = new QPushButton(layoutWidget);
        finalize_button->setObjectName("finalize_button");
        finalize_button->setMaximumSize(QSize(100, 16777215));

        save_skip->addWidget(finalize_button);


        horizontalLayout->addLayout(save_skip);

        ecg_button_layout = new QVBoxLayout();
        ecg_button_layout->setObjectName("ecg_button_layout");
        startNoiseECG = new QPushButton(layoutWidget);
        startNoiseECG->setObjectName("startNoiseECG");
        startNoiseECG->setMaximumSize(QSize(130, 16777215));

        ecg_button_layout->addWidget(startNoiseECG);

        stopNoiseECG = new QPushButton(layoutWidget);
        stopNoiseECG->setObjectName("stopNoiseECG");
        stopNoiseECG->setMaximumSize(QSize(130, 16777215));

        ecg_button_layout->addWidget(stopNoiseECG);


        horizontalLayout->addLayout(ecg_button_layout);

        ppg_startstop = new QVBoxLayout();
        ppg_startstop->setObjectName("ppg_startstop");
        startNoisePPG = new QPushButton(layoutWidget);
        startNoisePPG->setObjectName("startNoisePPG");
        startNoisePPG->setMaximumSize(QSize(130, 16777215));

        ppg_startstop->addWidget(startNoisePPG);

        stopNoisePPG = new QPushButton(layoutWidget);
        stopNoisePPG->setObjectName("stopNoisePPG");
        stopNoisePPG->setMaximumSize(QSize(130, 16777215));

        ppg_startstop->addWidget(stopNoisePPG);


        horizontalLayout->addLayout(ppg_startstop);

        marking_type_box = new QVBoxLayout();
        marking_type_box->setObjectName("marking_type_box");
        marking_type_label = new QLabel(layoutWidget);
        marking_type_label->setObjectName("marking_type_label");
        QSizePolicy sizePolicy9(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Minimum);
        sizePolicy9.setHorizontalStretch(0);
        sizePolicy9.setVerticalStretch(0);
        sizePolicy9.setHeightForWidth(marking_type_label->sizePolicy().hasHeightForWidth());
        marking_type_label->setSizePolicy(sizePolicy9);
        marking_type_label->setMaximumSize(QSize(100, 25));
        marking_type_label->setMargin(2);

        marking_type_box->addWidget(marking_type_label);

        marking_type = new QComboBox(layoutWidget);
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
        marking_type->setEnabled(true);
        QSizePolicy sizePolicy10(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Maximum);
        sizePolicy10.setHorizontalStretch(0);
        sizePolicy10.setVerticalStretch(0);
        sizePolicy10.setHeightForWidth(marking_type->sizePolicy().hasHeightForWidth());
        marking_type->setSizePolicy(sizePolicy10);
        marking_type->setMaximumSize(QSize(125, 16777215));

        marking_type_box->addWidget(marking_type);


        horizontalLayout->addLayout(marking_type_box);

        scroll_seconds_layout = new QVBoxLayout();
        scroll_seconds_layout->setObjectName("scroll_seconds_layout");
        seconds_to_scroll_label = new QLabel(layoutWidget);
        seconds_to_scroll_label->setObjectName("seconds_to_scroll_label");
        seconds_to_scroll_label->setMaximumSize(QSize(16777215, 25));

        scroll_seconds_layout->addWidget(seconds_to_scroll_label);

        skip_interval_box = new QLineEdit(layoutWidget);
        skip_interval_box->setObjectName("skip_interval_box");
        QSizePolicy sizePolicy11(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Fixed);
        sizePolicy11.setHorizontalStretch(0);
        sizePolicy11.setVerticalStretch(0);
        sizePolicy11.setHeightForWidth(skip_interval_box->sizePolicy().hasHeightForWidth());
        skip_interval_box->setSizePolicy(sizePolicy11);
        skip_interval_box->setMaximumSize(QSize(50, 25));
        skip_interval_box->setAutoFillBackground(false);

        scroll_seconds_layout->addWidget(skip_interval_box);


        horizontalLayout->addLayout(scroll_seconds_layout);

        ampogram_shift = new QVBoxLayout();
        ampogram_shift->setObjectName("ampogram_shift");
        prev8hours = new QPushButton(layoutWidget);
        prev8hours->setObjectName("prev8hours");

        ampogram_shift->addWidget(prev8hours);

        next8hours = new QPushButton(layoutWidget);
        next8hours->setObjectName("next8hours");

        ampogram_shift->addWidget(next8hours);


        horizontalLayout->addLayout(ampogram_shift);

        horizontalSpacer_2 = new QSpacerItem(40, 20, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        horizontalLayout->addItem(horizontalSpacer_2);


        verticalLayout->addLayout(horizontalLayout);

        verticalLayout->setStretch(2, 1);

        retranslateUi(noise_marking_gui);

        QMetaObject::connectSlotsByName(noise_marking_gui);
    } // setupUi

    void retranslateUi(QDialog *noise_marking_gui)
    {
        noise_marking_gui->setWindowTitle(QCoreApplication::translate("noise_marking_gui", "PPG/ECG Noise Marker", nullptr));
        topLabel->setText(QString());
        legend->setTitle(QCoreApplication::translate("noise_marking_gui", "Sleep States", nullptr));
        label_9->setText(QString());
        label_5->setText(QCoreApplication::translate("noise_marking_gui", "Wake", nullptr));
        label_8->setText(QString());
        label_4->setText(QCoreApplication::translate("noise_marking_gui", "REM", nullptr));
        label_7->setText(QString());
        label_3->setText(QCoreApplication::translate("noise_marking_gui", "NREM1", nullptr));
        label_6->setText(QString());
        label_2->setText(QCoreApplication::translate("noise_marking_gui", "NREM2", nullptr));
        wake_color->setText(QString());
        wake_text->setText(QCoreApplication::translate("noise_marking_gui", "NREM3", nullptr));
        hours_label->setText(QCoreApplication::translate("noise_marking_gui", "Time (h)", nullptr));
        selection_type_legend->setTitle(QCoreApplication::translate("noise_marking_gui", "Marking Type", nullptr));
        noise_color->setText(QString());
        noise_desc->setText(QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        label_19->setText(QString());
        label_20->setText(QCoreApplication::translate("noise_marking_gui", "Conduction Delay", nullptr));
        label_17->setText(QString());
        label_18->setText(QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        af_color->setText(QString());
        af_desc->setText(QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        label_15->setText(QString());
        label_16->setText(QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        wake_color_2->setText(QString());
        wake_text_2->setText(QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        pac_color->setText(QString());
        pac_desc->setText(QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        wake_color_5->setText(QString());
        wake_text_5->setText(QCoreApplication::translate("noise_marking_gui", "Benign Arrhythmia", nullptr));
        wake_color_6->setText(QString());
        wake_text_6->setText(QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));
        window_length->setTitle(QCoreApplication::translate("noise_marking_gui", "Window Length", nullptr));
        rb_10s->setText(QCoreApplication::translate("noise_marking_gui", "10s", nullptr));
        rb_30s->setText(QCoreApplication::translate("noise_marking_gui", "30s", nullptr));
        rb_1m->setText(QCoreApplication::translate("noise_marking_gui", "1m", nullptr));
        rb_5m->setText(QCoreApplication::translate("noise_marking_gui", "5m", nullptr));
        rb_10m->setText(QCoreApplication::translate("noise_marking_gui", "10m", nullptr));
        undo_button->setText(QCoreApplication::translate("noise_marking_gui", "Undo", nullptr));
        clearall_button->setText(QCoreApplication::translate("noise_marking_gui", "Clear All", nullptr));
        skip_button->setText(QCoreApplication::translate("noise_marking_gui", "Skip", nullptr));
        finalize_button->setText(QCoreApplication::translate("noise_marking_gui", "Save", nullptr));
        startNoiseECG->setText(QCoreApplication::translate("noise_marking_gui", "ECG Marking Start", nullptr));
        stopNoiseECG->setText(QCoreApplication::translate("noise_marking_gui", "ECG Marking End", nullptr));
        startNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "PPG Marking Start", nullptr));
        stopNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "PPG Marking End", nullptr));
        marking_type_label->setText(QCoreApplication::translate("noise_marking_gui", "  Marking Type", nullptr));
        marking_type->setItemText(0, QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        marking_type->setItemText(1, QCoreApplication::translate("noise_marking_gui", "Conduction Delay", nullptr));
        marking_type->setItemText(2, QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        marking_type->setItemText(3, QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        marking_type->setItemText(4, QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        marking_type->setItemText(5, QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        marking_type->setItemText(6, QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        marking_type->setItemText(7, QCoreApplication::translate("noise_marking_gui", "Benign Arrhythmia", nullptr));
        marking_type->setItemText(8, QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));

        seconds_to_scroll_label->setText(QCoreApplication::translate("noise_marking_gui", "Seconds to Scroll", nullptr));
        skip_interval_box->setText(QCoreApplication::translate("noise_marking_gui", "5.0", nullptr));
        prev8hours->setText(QCoreApplication::translate("noise_marking_gui", "Prev 8 Hours", nullptr));
        next8hours->setText(QCoreApplication::translate("noise_marking_gui", "Next 8 Hours", nullptr));
    } // retranslateUi

};

namespace Ui {
    class noise_marking_gui: public Ui_noise_marking_gui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_NOISE_MARKING_GUI_H
