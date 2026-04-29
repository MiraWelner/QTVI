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
#include <QtWidgets/QSpacerItem>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_noise_marking_gui
{
public:
    QVBoxLayout *dialogTopLayout;
    QHBoxLayout *file_selector_layout;
    QLabel *topLabel;
    QSpacerItem *horizontalSpacer;
    QPushButton *browse_file_button;
    QHBoxLayout *main_layout;
    QVBoxLayout *button_section;
    QVBoxLayout *ampogram_and_sleepstates;
    QChartView *ecg_ampogram_axis;
    QChartView *amp_ppg_axis;
    QLabel *label;
    QChartView *resp_cvp_axis;
    QHBoxLayout *window_length_buttons;
    QGridLayout *gridLayout_3;
    QPushButton *prev8hours;
    QComboBox *window_length_selector;
    QPushButton *start_all_mark;
    QPushButton *startNoisePPG;
    QLabel *window_length_label;
    QPushButton *stopNoiseABP;
    QPushButton *finalize_button;
    QPushButton *stop_all_mark;
    QPushButton *clearall_button;
    QPushButton *undo_button;
    QPushButton *skip_button;
    QPushButton *stopNoisePPG;
    QPushButton *startNoiseABP;
    QPushButton *next8hours;
    QLabel *label_4;
    QComboBox *scatter_line;
    QGridLayout *gridLayout_4;
    QGridLayout *marking_type_legend;
    QHBoxLayout *hboxLayout;
    QLabel *noise_color;
    QLabel *noise_desc;
    QHBoxLayout *hboxLayout1;
    QLabel *vt_color;
    QLabel *vt_desc;
    QHBoxLayout *hboxLayout2;
    QLabel *delay_color;
    QLabel *delay_desc;
    QHBoxLayout *hboxLayout3;
    QLabel *pac_color;
    QLabel *pac_desc;
    QHBoxLayout *hboxLayout4;
    QLabel *svt_color;
    QLabel *svt_desc;
    QHBoxLayout *hboxLayout5;
    QLabel *pvc_color;
    QLabel *pvc_desc;
    QHBoxLayout *hboxLayout6;
    QLabel *af_color;
    QLabel *af_desc;
    QHBoxLayout *hboxLayout7;
    QLabel *benign_color;
    QLabel *benign_desc;
    QHBoxLayout *hboxLayout8;
    QLabel *signif_color;
    QLabel *signif_desc;
    QLabel *marking_type_label_2;
    QLabel *sleep_states_label;
    QGridLayout *sleep_state_legend;
    QHBoxLayout *hboxLayout9;
    QLabel *wake_c;
    QLabel *wake_l;
    QHBoxLayout *hboxLayout10;
    QLabel *n3_c;
    QLabel *n3_l;
    QHBoxLayout *hboxLayout11;
    QLabel *n2_c;
    QLabel *n2_l;
    QHBoxLayout *hboxLayout12;
    QLabel *n1_c;
    QLabel *n1_l;
    QHBoxLayout *hboxLayout13;
    QLabel *rem_c;
    QLabel *rem_l;
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *end_ecg_all;
    QPushButton *stop_ecg1_mark;
    QPushButton *stop_ecg2_mark;
    QPushButton *stop_ecg3_mark;
    QLineEdit *skip_interval_box;
    QLabel *seconds_to_scroll_label;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *start_ecg_all;
    QPushButton *start_ecg1_mark;
    QPushButton *start_ecg2_mark;
    QPushButton *start_ecg3_mark;
    QLabel *marking_type_label;
    QLabel *label_3;
    QComboBox *marking_type;
    QLabel *label_2;
    QVBoxLayout *main_plots;
    QChartView *ecg_axis_1;
    QChartView *ecg_axis_2;
    QChartView *ecg_axis_3;
    QChartView *ppg_axis;
    QChartView *accel_or_abg_axis;

    void setupUi(QDialog *noise_marking_gui)
    {
        if (noise_marking_gui->objectName().isEmpty())
            noise_marking_gui->setObjectName("noise_marking_gui");
        noise_marking_gui->resize(1445, 1200);
        dialogTopLayout = new QVBoxLayout(noise_marking_gui);
        dialogTopLayout->setSpacing(0);
        dialogTopLayout->setObjectName("dialogTopLayout");
        dialogTopLayout->setContentsMargins(10, 10, 10, 10);
        file_selector_layout = new QHBoxLayout();
        file_selector_layout->setObjectName("file_selector_layout");
        topLabel = new QLabel(noise_marking_gui);
        topLabel->setObjectName("topLabel");

        file_selector_layout->addWidget(topLabel);

        horizontalSpacer = new QSpacerItem(40, 20, QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Minimum);

        file_selector_layout->addItem(horizontalSpacer);

        browse_file_button = new QPushButton(noise_marking_gui);
        browse_file_button->setObjectName("browse_file_button");

        file_selector_layout->addWidget(browse_file_button);


        dialogTopLayout->addLayout(file_selector_layout);

        main_layout = new QHBoxLayout();
        main_layout->setSpacing(0);
        main_layout->setObjectName("main_layout");
        button_section = new QVBoxLayout();
        button_section->setSpacing(10);
        button_section->setObjectName("button_section");
        button_section->setContentsMargins(5, 5, 5, 5);
        ampogram_and_sleepstates = new QVBoxLayout();
        ampogram_and_sleepstates->setSpacing(0);
        ampogram_and_sleepstates->setObjectName("ampogram_and_sleepstates");
        ecg_ampogram_axis = new QChartView(noise_marking_gui);
        ecg_ampogram_axis->setObjectName("ecg_ampogram_axis");
        QSizePolicy sizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(ecg_ampogram_axis->sizePolicy().hasHeightForWidth());
        ecg_ampogram_axis->setSizePolicy(sizePolicy);

        ampogram_and_sleepstates->addWidget(ecg_ampogram_axis);

        amp_ppg_axis = new QChartView(noise_marking_gui);
        amp_ppg_axis->setObjectName("amp_ppg_axis");
        sizePolicy.setHeightForWidth(amp_ppg_axis->sizePolicy().hasHeightForWidth());
        amp_ppg_axis->setSizePolicy(sizePolicy);

        ampogram_and_sleepstates->addWidget(amp_ppg_axis);

        label = new QLabel(noise_marking_gui);
        label->setObjectName("label");
        label->setStyleSheet(QString::fromUtf8("border:none;"));

        ampogram_and_sleepstates->addWidget(label);

        resp_cvp_axis = new QChartView(noise_marking_gui);
        resp_cvp_axis->setObjectName("resp_cvp_axis");
        sizePolicy.setHeightForWidth(resp_cvp_axis->sizePolicy().hasHeightForWidth());
        resp_cvp_axis->setSizePolicy(sizePolicy);

        ampogram_and_sleepstates->addWidget(resp_cvp_axis);

        ampogram_and_sleepstates->setStretch(0, 1);
        ampogram_and_sleepstates->setStretch(1, 1);

        button_section->addLayout(ampogram_and_sleepstates);

        window_length_buttons = new QHBoxLayout();
        window_length_buttons->setObjectName("window_length_buttons");
        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName("gridLayout_3");
        prev8hours = new QPushButton(noise_marking_gui);
        prev8hours->setObjectName("prev8hours");
        prev8hours->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(prev8hours, 0, 2, 1, 1);

        window_length_selector = new QComboBox(noise_marking_gui);
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->setObjectName("window_length_selector");

        gridLayout_3->addWidget(window_length_selector, 1, 0, 1, 1);

        start_all_mark = new QPushButton(noise_marking_gui);
        start_all_mark->setObjectName("start_all_mark");
        start_all_mark->setMinimumSize(QSize(0, 0));
        start_all_mark->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_3->addWidget(start_all_mark, 0, 5, 1, 1);

        startNoisePPG = new QPushButton(noise_marking_gui);
        startNoisePPG->setObjectName("startNoisePPG");

        gridLayout_3->addWidget(startNoisePPG, 0, 6, 1, 1);

        window_length_label = new QLabel(noise_marking_gui);
        window_length_label->setObjectName("window_length_label");
        QFont font;
        font.setBold(true);
        window_length_label->setFont(font);

        gridLayout_3->addWidget(window_length_label, 0, 0, 1, 1);

        stopNoiseABP = new QPushButton(noise_marking_gui);
        stopNoiseABP->setObjectName("stopNoiseABP");

        gridLayout_3->addWidget(stopNoiseABP, 1, 7, 1, 1);

        finalize_button = new QPushButton(noise_marking_gui);
        finalize_button->setObjectName("finalize_button");
        finalize_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(finalize_button, 1, 4, 1, 1);

        stop_all_mark = new QPushButton(noise_marking_gui);
        stop_all_mark->setObjectName("stop_all_mark");
        stop_all_mark->setMinimumSize(QSize(10, 0));
        stop_all_mark->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_3->addWidget(stop_all_mark, 1, 5, 1, 1);

        clearall_button = new QPushButton(noise_marking_gui);
        clearall_button->setObjectName("clearall_button");
        clearall_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(clearall_button, 1, 3, 1, 1);

        undo_button = new QPushButton(noise_marking_gui);
        undo_button->setObjectName("undo_button");
        undo_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(undo_button, 0, 3, 1, 1);

        skip_button = new QPushButton(noise_marking_gui);
        skip_button->setObjectName("skip_button");
        skip_button->setMinimumSize(QSize(30, 0));
        skip_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(skip_button, 0, 4, 1, 1);

        stopNoisePPG = new QPushButton(noise_marking_gui);
        stopNoisePPG->setObjectName("stopNoisePPG");

        gridLayout_3->addWidget(stopNoisePPG, 1, 6, 1, 1);

        startNoiseABP = new QPushButton(noise_marking_gui);
        startNoiseABP->setObjectName("startNoiseABP");

        gridLayout_3->addWidget(startNoiseABP, 0, 7, 1, 1);

        next8hours = new QPushButton(noise_marking_gui);
        next8hours->setObjectName("next8hours");
        next8hours->setMaximumSize(QSize(47, 16777215));

        gridLayout_3->addWidget(next8hours, 1, 2, 1, 1);

        label_4 = new QLabel(noise_marking_gui);
        label_4->setObjectName("label_4");
        label_4->setFont(font);

        gridLayout_3->addWidget(label_4, 0, 1, 1, 1);

        scatter_line = new QComboBox(noise_marking_gui);
        scatter_line->addItem(QString());
        scatter_line->addItem(QString());
        scatter_line->setObjectName("scatter_line");

        gridLayout_3->addWidget(scatter_line, 1, 1, 1, 1);


        window_length_buttons->addLayout(gridLayout_3);


        button_section->addLayout(window_length_buttons);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName("gridLayout_4");
        marking_type_legend = new QGridLayout();
        marking_type_legend->setSpacing(4);
        marking_type_legend->setObjectName("marking_type_legend");
        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName("hboxLayout");
        noise_color = new QLabel(noise_marking_gui);
        noise_color->setObjectName("noise_color");
        noise_color->setMinimumSize(QSize(15, 15));
        noise_color->setMaximumSize(QSize(15, 15));

        hboxLayout->addWidget(noise_color);

        noise_desc = new QLabel(noise_marking_gui);
        noise_desc->setObjectName("noise_desc");

        hboxLayout->addWidget(noise_desc);


        marking_type_legend->addLayout(hboxLayout, 0, 0, 1, 1);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setObjectName("hboxLayout1");
        vt_color = new QLabel(noise_marking_gui);
        vt_color->setObjectName("vt_color");
        vt_color->setMinimumSize(QSize(15, 15));
        vt_color->setMaximumSize(QSize(15, 15));

        hboxLayout1->addWidget(vt_color);

        vt_desc = new QLabel(noise_marking_gui);
        vt_desc->setObjectName("vt_desc");

        hboxLayout1->addWidget(vt_desc);


        marking_type_legend->addLayout(hboxLayout1, 0, 1, 1, 1);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setObjectName("hboxLayout2");
        delay_color = new QLabel(noise_marking_gui);
        delay_color->setObjectName("delay_color");
        delay_color->setMinimumSize(QSize(15, 15));
        delay_color->setMaximumSize(QSize(15, 15));

        hboxLayout2->addWidget(delay_color);

        delay_desc = new QLabel(noise_marking_gui);
        delay_desc->setObjectName("delay_desc");

        hboxLayout2->addWidget(delay_desc);


        marking_type_legend->addLayout(hboxLayout2, 1, 0, 1, 1);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setObjectName("hboxLayout3");
        pac_color = new QLabel(noise_marking_gui);
        pac_color->setObjectName("pac_color");
        pac_color->setMinimumSize(QSize(15, 15));
        pac_color->setMaximumSize(QSize(15, 15));

        hboxLayout3->addWidget(pac_color);

        pac_desc = new QLabel(noise_marking_gui);
        pac_desc->setObjectName("pac_desc");

        hboxLayout3->addWidget(pac_desc);


        marking_type_legend->addLayout(hboxLayout3, 1, 1, 1, 1);

        hboxLayout4 = new QHBoxLayout();
        hboxLayout4->setObjectName("hboxLayout4");
        svt_color = new QLabel(noise_marking_gui);
        svt_color->setObjectName("svt_color");
        svt_color->setMinimumSize(QSize(15, 15));
        svt_color->setMaximumSize(QSize(15, 15));

        hboxLayout4->addWidget(svt_color);

        svt_desc = new QLabel(noise_marking_gui);
        svt_desc->setObjectName("svt_desc");

        hboxLayout4->addWidget(svt_desc);


        marking_type_legend->addLayout(hboxLayout4, 2, 0, 1, 1);

        hboxLayout5 = new QHBoxLayout();
        hboxLayout5->setObjectName("hboxLayout5");
        pvc_color = new QLabel(noise_marking_gui);
        pvc_color->setObjectName("pvc_color");
        pvc_color->setMinimumSize(QSize(15, 15));
        pvc_color->setMaximumSize(QSize(15, 15));

        hboxLayout5->addWidget(pvc_color);

        pvc_desc = new QLabel(noise_marking_gui);
        pvc_desc->setObjectName("pvc_desc");

        hboxLayout5->addWidget(pvc_desc);


        marking_type_legend->addLayout(hboxLayout5, 2, 1, 1, 1);

        hboxLayout6 = new QHBoxLayout();
        hboxLayout6->setObjectName("hboxLayout6");
        af_color = new QLabel(noise_marking_gui);
        af_color->setObjectName("af_color");
        af_color->setMinimumSize(QSize(15, 15));
        af_color->setMaximumSize(QSize(15, 15));

        hboxLayout6->addWidget(af_color);

        af_desc = new QLabel(noise_marking_gui);
        af_desc->setObjectName("af_desc");

        hboxLayout6->addWidget(af_desc);


        marking_type_legend->addLayout(hboxLayout6, 3, 0, 1, 1);

        hboxLayout7 = new QHBoxLayout();
        hboxLayout7->setObjectName("hboxLayout7");
        benign_color = new QLabel(noise_marking_gui);
        benign_color->setObjectName("benign_color");
        benign_color->setMinimumSize(QSize(15, 15));
        benign_color->setMaximumSize(QSize(15, 15));

        hboxLayout7->addWidget(benign_color);

        benign_desc = new QLabel(noise_marking_gui);
        benign_desc->setObjectName("benign_desc");

        hboxLayout7->addWidget(benign_desc);


        marking_type_legend->addLayout(hboxLayout7, 3, 1, 1, 1);

        hboxLayout8 = new QHBoxLayout();
        hboxLayout8->setObjectName("hboxLayout8");
        signif_color = new QLabel(noise_marking_gui);
        signif_color->setObjectName("signif_color");
        signif_color->setMinimumSize(QSize(15, 15));
        signif_color->setMaximumSize(QSize(15, 15));

        hboxLayout8->addWidget(signif_color);

        signif_desc = new QLabel(noise_marking_gui);
        signif_desc->setObjectName("signif_desc");

        hboxLayout8->addWidget(signif_desc);


        marking_type_legend->addLayout(hboxLayout8, 4, 0, 1, 2);


        gridLayout_4->addLayout(marking_type_legend, 1, 0, 1, 1);

        marking_type_label_2 = new QLabel(noise_marking_gui);
        marking_type_label_2->setObjectName("marking_type_label_2");
        marking_type_label_2->setFont(font);

        gridLayout_4->addWidget(marking_type_label_2, 0, 0, 1, 1);

        sleep_states_label = new QLabel(noise_marking_gui);
        sleep_states_label->setObjectName("sleep_states_label");
        sleep_states_label->setFont(font);

        gridLayout_4->addWidget(sleep_states_label, 0, 2, 1, 1);

        sleep_state_legend = new QGridLayout();
        sleep_state_legend->setSpacing(4);
        sleep_state_legend->setObjectName("sleep_state_legend");
        hboxLayout9 = new QHBoxLayout();
        hboxLayout9->setObjectName("hboxLayout9");
        wake_c = new QLabel(noise_marking_gui);
        wake_c->setObjectName("wake_c");
        wake_c->setMinimumSize(QSize(15, 15));
        wake_c->setMaximumSize(QSize(15, 15));

        hboxLayout9->addWidget(wake_c);

        wake_l = new QLabel(noise_marking_gui);
        wake_l->setObjectName("wake_l");

        hboxLayout9->addWidget(wake_l);


        sleep_state_legend->addLayout(hboxLayout9, 0, 0, 1, 1);

        hboxLayout10 = new QHBoxLayout();
        hboxLayout10->setObjectName("hboxLayout10");
        n3_c = new QLabel(noise_marking_gui);
        n3_c->setObjectName("n3_c");
        n3_c->setMinimumSize(QSize(15, 15));
        n3_c->setMaximumSize(QSize(15, 15));

        hboxLayout10->addWidget(n3_c);

        n3_l = new QLabel(noise_marking_gui);
        n3_l->setObjectName("n3_l");

        hboxLayout10->addWidget(n3_l);


        sleep_state_legend->addLayout(hboxLayout10, 2, 0, 1, 1);

        hboxLayout11 = new QHBoxLayout();
        hboxLayout11->setObjectName("hboxLayout11");
        n2_c = new QLabel(noise_marking_gui);
        n2_c->setObjectName("n2_c");
        n2_c->setMinimumSize(QSize(15, 15));
        n2_c->setMaximumSize(QSize(15, 15));

        hboxLayout11->addWidget(n2_c);

        n2_l = new QLabel(noise_marking_gui);
        n2_l->setObjectName("n2_l");

        hboxLayout11->addWidget(n2_l);


        sleep_state_legend->addLayout(hboxLayout11, 1, 1, 1, 1);

        hboxLayout12 = new QHBoxLayout();
        hboxLayout12->setObjectName("hboxLayout12");
        n1_c = new QLabel(noise_marking_gui);
        n1_c->setObjectName("n1_c");
        n1_c->setMinimumSize(QSize(15, 15));
        n1_c->setMaximumSize(QSize(15, 15));

        hboxLayout12->addWidget(n1_c);

        n1_l = new QLabel(noise_marking_gui);
        n1_l->setObjectName("n1_l");

        hboxLayout12->addWidget(n1_l);


        sleep_state_legend->addLayout(hboxLayout12, 1, 0, 1, 1);

        hboxLayout13 = new QHBoxLayout();
        hboxLayout13->setObjectName("hboxLayout13");
        rem_c = new QLabel(noise_marking_gui);
        rem_c->setObjectName("rem_c");
        rem_c->setMinimumSize(QSize(15, 15));
        rem_c->setMaximumSize(QSize(15, 15));

        hboxLayout13->addWidget(rem_c);

        rem_l = new QLabel(noise_marking_gui);
        rem_l->setObjectName("rem_l");

        hboxLayout13->addWidget(rem_l);


        sleep_state_legend->addLayout(hboxLayout13, 0, 1, 1, 1);


        gridLayout_4->addLayout(sleep_state_legend, 1, 2, 1, 1);

        gridLayout = new QGridLayout();
        gridLayout->setObjectName("gridLayout");
        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName("horizontalLayout_4");
        end_ecg_all = new QPushButton(noise_marking_gui);
        end_ecg_all->setObjectName("end_ecg_all");
        end_ecg_all->setMinimumSize(QSize(0, 0));
        end_ecg_all->setMaximumSize(QSize(50, 16777215));

        horizontalLayout_4->addWidget(end_ecg_all);

        stop_ecg1_mark = new QPushButton(noise_marking_gui);
        stop_ecg1_mark->setObjectName("stop_ecg1_mark");
        stop_ecg1_mark->setMinimumSize(QSize(24, 0));
        stop_ecg1_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_4->addWidget(stop_ecg1_mark);

        stop_ecg2_mark = new QPushButton(noise_marking_gui);
        stop_ecg2_mark->setObjectName("stop_ecg2_mark");
        stop_ecg2_mark->setMinimumSize(QSize(24, 0));
        stop_ecg2_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_4->addWidget(stop_ecg2_mark);

        stop_ecg3_mark = new QPushButton(noise_marking_gui);
        stop_ecg3_mark->setObjectName("stop_ecg3_mark");
        stop_ecg3_mark->setMinimumSize(QSize(24, 0));
        stop_ecg3_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_4->addWidget(stop_ecg3_mark);


        gridLayout->addLayout(horizontalLayout_4, 3, 1, 1, 1);

        skip_interval_box = new QLineEdit(noise_marking_gui);
        skip_interval_box->setObjectName("skip_interval_box");
        skip_interval_box->setMaximumSize(QSize(50, 25));

        gridLayout->addWidget(skip_interval_box, 1, 0, 1, 1);

        seconds_to_scroll_label = new QLabel(noise_marking_gui);
        seconds_to_scroll_label->setObjectName("seconds_to_scroll_label");
        seconds_to_scroll_label->setFont(font);

        gridLayout->addWidget(seconds_to_scroll_label, 0, 0, 1, 1);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName("horizontalLayout_2");
        start_ecg_all = new QPushButton(noise_marking_gui);
        start_ecg_all->setObjectName("start_ecg_all");
        start_ecg_all->setMaximumSize(QSize(50, 16777215));

        horizontalLayout_2->addWidget(start_ecg_all);

        start_ecg1_mark = new QPushButton(noise_marking_gui);
        start_ecg1_mark->setObjectName("start_ecg1_mark");
        start_ecg1_mark->setMinimumSize(QSize(24, 0));
        start_ecg1_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_2->addWidget(start_ecg1_mark);

        start_ecg2_mark = new QPushButton(noise_marking_gui);
        start_ecg2_mark->setObjectName("start_ecg2_mark");
        start_ecg2_mark->setMinimumSize(QSize(24, 0));
        start_ecg2_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_2->addWidget(start_ecg2_mark);

        start_ecg3_mark = new QPushButton(noise_marking_gui);
        start_ecg3_mark->setObjectName("start_ecg3_mark");
        start_ecg3_mark->setMinimumSize(QSize(24, 0));
        start_ecg3_mark->setMaximumSize(QSize(24, 16777215));

        horizontalLayout_2->addWidget(start_ecg3_mark);


        gridLayout->addLayout(horizontalLayout_2, 3, 0, 1, 1);

        marking_type_label = new QLabel(noise_marking_gui);
        marking_type_label->setObjectName("marking_type_label");
        marking_type_label->setFont(font);

        gridLayout->addWidget(marking_type_label, 0, 1, 1, 1);

        label_3 = new QLabel(noise_marking_gui);
        label_3->setObjectName("label_3");
        label_3->setFont(font);

        gridLayout->addWidget(label_3, 2, 1, 1, 1);

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
        marking_type->setMaximumSize(QSize(100, 16777215));

        gridLayout->addWidget(marking_type, 1, 1, 1, 1);

        label_2 = new QLabel(noise_marking_gui);
        label_2->setObjectName("label_2");
        label_2->setFont(font);

        gridLayout->addWidget(label_2, 2, 0, 1, 1);


        gridLayout_4->addLayout(gridLayout, 1, 1, 1, 1);


        button_section->addLayout(gridLayout_4);


        main_layout->addLayout(button_section);

        main_plots = new QVBoxLayout();
        main_plots->setSpacing(0);
        main_plots->setObjectName("main_plots");
        ecg_axis_1 = new QChartView(noise_marking_gui);
        ecg_axis_1->setObjectName("ecg_axis_1");
        QSizePolicy sizePolicy1(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Expanding);
        sizePolicy1.setHorizontalStretch(0);
        sizePolicy1.setVerticalStretch(1);
        sizePolicy1.setHeightForWidth(ecg_axis_1->sizePolicy().hasHeightForWidth());
        ecg_axis_1->setSizePolicy(sizePolicy1);
        ecg_axis_1->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));

        main_plots->addWidget(ecg_axis_1);

        ecg_axis_2 = new QChartView(noise_marking_gui);
        ecg_axis_2->setObjectName("ecg_axis_2");
        sizePolicy1.setHeightForWidth(ecg_axis_2->sizePolicy().hasHeightForWidth());
        ecg_axis_2->setSizePolicy(sizePolicy1);
        ecg_axis_2->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));

        main_plots->addWidget(ecg_axis_2);

        ecg_axis_3 = new QChartView(noise_marking_gui);
        ecg_axis_3->setObjectName("ecg_axis_3");
        sizePolicy1.setHeightForWidth(ecg_axis_3->sizePolicy().hasHeightForWidth());
        ecg_axis_3->setSizePolicy(sizePolicy1);
        ecg_axis_3->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));

        main_plots->addWidget(ecg_axis_3);

        ppg_axis = new QChartView(noise_marking_gui);
        ppg_axis->setObjectName("ppg_axis");
        sizePolicy1.setHeightForWidth(ppg_axis->sizePolicy().hasHeightForWidth());
        ppg_axis->setSizePolicy(sizePolicy1);
        ppg_axis->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));

        main_plots->addWidget(ppg_axis);

        accel_or_abg_axis = new QChartView(noise_marking_gui);
        accel_or_abg_axis->setObjectName("accel_or_abg_axis");
        sizePolicy1.setHeightForWidth(accel_or_abg_axis->sizePolicy().hasHeightForWidth());
        accel_or_abg_axis->setSizePolicy(sizePolicy1);

        main_plots->addWidget(accel_or_abg_axis);


        main_layout->addLayout(main_plots);

        main_layout->setStretch(0, 1);
        main_layout->setStretch(1, 2);

        dialogTopLayout->addLayout(main_layout);


        retranslateUi(noise_marking_gui);

        QMetaObject::connectSlotsByName(noise_marking_gui);
    } // setupUi

    void retranslateUi(QDialog *noise_marking_gui)
    {
        noise_marking_gui->setWindowTitle(QCoreApplication::translate("noise_marking_gui", "PPG/ECG Noise Marker", nullptr));
        topLabel->setText(QString());
        browse_file_button->setText(QCoreApplication::translate("noise_marking_gui", "Browse Different File", nullptr));
        ecg_ampogram_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        amp_ppg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        label->setText(QCoreApplication::translate("noise_marking_gui", "Time (h) - 8 hours", nullptr));
        resp_cvp_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        prev8hours->setText(QCoreApplication::translate("noise_marking_gui", "Prev 8h", nullptr));
        window_length_selector->setItemText(0, QCoreApplication::translate("noise_marking_gui", "1 Second", nullptr));
        window_length_selector->setItemText(1, QCoreApplication::translate("noise_marking_gui", "3 Second", nullptr));
        window_length_selector->setItemText(2, QCoreApplication::translate("noise_marking_gui", "10 Seconds", nullptr));
        window_length_selector->setItemText(3, QCoreApplication::translate("noise_marking_gui", "30 Seconds", nullptr));
        window_length_selector->setItemText(4, QCoreApplication::translate("noise_marking_gui", "1 Minute", nullptr));
        window_length_selector->setItemText(5, QCoreApplication::translate("noise_marking_gui", "10 Minutes", nullptr));

        start_all_mark->setText(QCoreApplication::translate("noise_marking_gui", "Mark All Start", nullptr));
        startNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG Start", nullptr));
        window_length_label->setText(QCoreApplication::translate("noise_marking_gui", "Window Length:", nullptr));
        stopNoiseABP->setText(QCoreApplication::translate("noise_marking_gui", "Mark aBP End", nullptr));
        finalize_button->setText(QCoreApplication::translate("noise_marking_gui", "Save", nullptr));
        stop_all_mark->setText(QCoreApplication::translate("noise_marking_gui", "Mark All End", nullptr));
        clearall_button->setText(QCoreApplication::translate("noise_marking_gui", "Clear", nullptr));
        undo_button->setText(QCoreApplication::translate("noise_marking_gui", "Undo", nullptr));
        skip_button->setText(QCoreApplication::translate("noise_marking_gui", "Skip", nullptr));
        stopNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG End", nullptr));
        startNoiseABP->setText(QCoreApplication::translate("noise_marking_gui", "Mark aBP Start", nullptr));
        next8hours->setText(QCoreApplication::translate("noise_marking_gui", "Next 8h", nullptr));
        label_4->setText(QCoreApplication::translate("noise_marking_gui", "Plotting Style", nullptr));
        scatter_line->setItemText(0, QCoreApplication::translate("noise_marking_gui", "Line Plot", nullptr));
        scatter_line->setItemText(1, QCoreApplication::translate("noise_marking_gui", "Scatter Plot", nullptr));

        noise_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: yellow;", nullptr));
        noise_desc->setText(QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        vt_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: green; border-radius:3px;", nullptr));
        vt_desc->setText(QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        delay_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(128,0,128); border-radius:3px;", nullptr));
        delay_desc->setText(QCoreApplication::translate("noise_marking_gui", "Cond. Delay", nullptr));
        pac_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(255,128, 0); border-radius:3px;", nullptr));
        pac_desc->setText(QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        svt_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: blue; border-radius:3px;", nullptr));
        svt_desc->setText(QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        pvc_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(128, 255, 0); border-radius:3px;", nullptr));
        pvc_desc->setText(QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        af_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: red; border-radius:3px;", nullptr));
        af_desc->setText(QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        benign_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(255, 128, 255); border-radius:3px;", nullptr));
        benign_desc->setText(QCoreApplication::translate("noise_marking_gui", "Benign", nullptr));
        signif_color->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: rgb(0, 255, 255); border-radius:3px;", nullptr));
        signif_desc->setText(QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));
        marking_type_label_2->setText(QCoreApplication::translate("noise_marking_gui", "Marking Type Legend", nullptr));
        sleep_states_label->setText(QCoreApplication::translate("noise_marking_gui", "Sleep States Legend", nullptr));
        wake_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: black; border-radius:3px;", nullptr));
        wake_l->setText(QCoreApplication::translate("noise_marking_gui", "Wake", nullptr));
        n3_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: red; border-radius:3px;", nullptr));
        n3_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM3", nullptr));
        n2_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: cyan; border-radius:3px;", nullptr));
        n2_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM2", nullptr));
        n1_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: blue; border-radius:3px;", nullptr));
        n1_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM1", nullptr));
        rem_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: green; border-radius:3px;", nullptr));
        rem_l->setText(QCoreApplication::translate("noise_marking_gui", "REM", nullptr));
        end_ecg_all->setText(QCoreApplication::translate("noise_marking_gui", "All ECG", nullptr));
        stop_ecg1_mark->setText(QCoreApplication::translate("noise_marking_gui", "1", nullptr));
        stop_ecg2_mark->setText(QCoreApplication::translate("noise_marking_gui", "2", nullptr));
        stop_ecg3_mark->setText(QCoreApplication::translate("noise_marking_gui", "3", nullptr));
        skip_interval_box->setText(QCoreApplication::translate("noise_marking_gui", "5.0", nullptr));
        seconds_to_scroll_label->setText(QCoreApplication::translate("noise_marking_gui", "Scroll Dist", nullptr));
        start_ecg_all->setText(QCoreApplication::translate("noise_marking_gui", "All ECG", nullptr));
        start_ecg1_mark->setText(QCoreApplication::translate("noise_marking_gui", "1", nullptr));
        start_ecg2_mark->setText(QCoreApplication::translate("noise_marking_gui", "2", nullptr));
        start_ecg3_mark->setText(QCoreApplication::translate("noise_marking_gui", "3", nullptr));
        marking_type_label->setText(QCoreApplication::translate("noise_marking_gui", "Mark Type", nullptr));
        label_3->setText(QCoreApplication::translate("noise_marking_gui", "Mark ECG Channel Stop", nullptr));
        marking_type->setItemText(0, QCoreApplication::translate("noise_marking_gui", "Noise/Artifact", nullptr));
        marking_type->setItemText(1, QCoreApplication::translate("noise_marking_gui", "Conduction Delay", nullptr));
        marking_type->setItemText(2, QCoreApplication::translate("noise_marking_gui", "AF", nullptr));
        marking_type->setItemText(3, QCoreApplication::translate("noise_marking_gui", "SVT", nullptr));
        marking_type->setItemText(4, QCoreApplication::translate("noise_marking_gui", "VT", nullptr));
        marking_type->setItemText(5, QCoreApplication::translate("noise_marking_gui", "PVC", nullptr));
        marking_type->setItemText(6, QCoreApplication::translate("noise_marking_gui", "PAC", nullptr));
        marking_type->setItemText(7, QCoreApplication::translate("noise_marking_gui", "Benign Arrhythmia", nullptr));
        marking_type->setItemText(8, QCoreApplication::translate("noise_marking_gui", "Significant Arrhythmia", nullptr));

        label_2->setText(QCoreApplication::translate("noise_marking_gui", "Mark ECG Channel Start", nullptr));
        accel_or_abg_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid black;", nullptr));
    } // retranslateUi

};

namespace Ui {
    class noise_marking_gui: public Ui_noise_marking_gui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_NOISE_MARKING_GUI_H
