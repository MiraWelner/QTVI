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
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_noise_marking_gui
{
public:
    QVBoxLayout *verticalLayout;
    QHBoxLayout *buttons_and_ampograms;
    QGridLayout *buttons;
    QGridLayout *gridLayout_5;
    QPushButton *startNoiseABP;
    QPushButton *start_all_mark;
    QPushButton *stop_all_mark;
    QPushButton *stopNoiseABP;
    QPushButton *stopNoisePPG;
    QPushButton *startNoisePPG;
    QGridLayout *gridLayout_2;
    QPushButton *next8hours;
    QPushButton *prev8hours;
    QPushButton *undo_button;
    QPushButton *clearall_button;
    QPushButton *finalize_button;
    QPushButton *skip_button;
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
    QLabel *marking_type_label_2;
    QGridLayout *sleep_state_legend;
    QHBoxLayout *hboxLayout;
    QLabel *rem_c;
    QLabel *rem_l;
    QHBoxLayout *hboxLayout1;
    QLabel *n1_c;
    QLabel *n1_l;
    QHBoxLayout *hboxLayout2;
    QLabel *n3_c;
    QLabel *n3_l;
    QHBoxLayout *hboxLayout3;
    QLabel *n2_c;
    QLabel *n2_l;
    QHBoxLayout *hboxLayout4;
    QLabel *wake_c;
    QLabel *wake_l;
    QLabel *sleep_states_label;
    QGridLayout *gridLayout_3;
    QLabel *window_length_label;
    QComboBox *window_length_selector;
    QComboBox *scatter_line;
    QLabel *label_4;
    QGridLayout *marking_type_legend;
    QHBoxLayout *hboxLayout5;
    QLabel *noise_color;
    QLabel *noise_desc;
    QHBoxLayout *hboxLayout6;
    QLabel *vt_color;
    QLabel *vt_desc;
    QHBoxLayout *hboxLayout7;
    QLabel *delay_color;
    QLabel *delay_desc;
    QHBoxLayout *hboxLayout8;
    QLabel *pac_color;
    QLabel *pac_desc;
    QHBoxLayout *hboxLayout9;
    QLabel *svt_color;
    QLabel *svt_desc;
    QHBoxLayout *hboxLayout10;
    QLabel *pvc_color;
    QLabel *pvc_desc;
    QHBoxLayout *hboxLayout11;
    QLabel *af_color;
    QLabel *af_desc;
    QHBoxLayout *hboxLayout12;
    QLabel *benign_color;
    QLabel *benign_desc;
    QHBoxLayout *hboxLayout13;
    QLabel *signif_color;
    QLabel *signif_desc;
    QGridLayout *gridLayout_4;
    QPushButton *process_button;
    QPushButton *browse_file_button;
    QPushButton *save_current_plot;
    QCheckBox *show_peaks_check;
    QVBoxLayout *ampogram_and_sleepstates;
    QChartView *ecg_ampogram_axis;
    QChartView *ppg_ampogram_axis;
    QChartView *hyp_accel_resp_cvp_axis;
    QVBoxLayout *main_plots;
    QChartView *ecg_axis_1;
    QCheckBox *ecg_1_check;
    QDoubleSpinBox *ecg_1_gain;
    QChartView *ecg_axis_2;
    QDoubleSpinBox *ecg_2_gain;
    QCheckBox *ecg_2_check;
    QChartView *ecg_axis_3;
    QCheckBox *ecg_3_check;
    QDoubleSpinBox *ecg_3_gain;
    QChartView *ppg_axis;
    QDoubleSpinBox *ppg_gain;
    QCheckBox *ppg_check;
    QChartView *accel_or_abp_axis;
    QDoubleSpinBox *abp_gain;
    QCheckBox *abp_check;

    void setupUi(QDialog *noise_marking_gui)
    {
        if (noise_marking_gui->objectName().isEmpty())
            noise_marking_gui->setObjectName("noise_marking_gui");
        noise_marking_gui->resize(1445, 1200);
        verticalLayout = new QVBoxLayout(noise_marking_gui);
        verticalLayout->setSpacing(3);
        verticalLayout->setObjectName("verticalLayout");
        verticalLayout->setContentsMargins(3, 3, 3, 3);
        buttons_and_ampograms = new QHBoxLayout();
        buttons_and_ampograms->setObjectName("buttons_and_ampograms");
        buttons = new QGridLayout();
        buttons->setObjectName("buttons");
        gridLayout_5 = new QGridLayout();
        gridLayout_5->setObjectName("gridLayout_5");
        startNoiseABP = new QPushButton(noise_marking_gui);
        startNoiseABP->setObjectName("startNoiseABP");

        gridLayout_5->addWidget(startNoiseABP, 1, 2, 1, 1);

        start_all_mark = new QPushButton(noise_marking_gui);
        start_all_mark->setObjectName("start_all_mark");
        start_all_mark->setMinimumSize(QSize(0, 0));
        start_all_mark->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_5->addWidget(start_all_mark, 1, 0, 1, 1);

        stop_all_mark = new QPushButton(noise_marking_gui);
        stop_all_mark->setObjectName("stop_all_mark");
        stop_all_mark->setMinimumSize(QSize(10, 0));
        stop_all_mark->setMaximumSize(QSize(16777215, 16777215));

        gridLayout_5->addWidget(stop_all_mark, 2, 0, 1, 1);

        stopNoiseABP = new QPushButton(noise_marking_gui);
        stopNoiseABP->setObjectName("stopNoiseABP");

        gridLayout_5->addWidget(stopNoiseABP, 2, 2, 1, 1);

        stopNoisePPG = new QPushButton(noise_marking_gui);
        stopNoisePPG->setObjectName("stopNoisePPG");

        gridLayout_5->addWidget(stopNoisePPG, 2, 1, 1, 1);

        startNoisePPG = new QPushButton(noise_marking_gui);
        startNoisePPG->setObjectName("startNoisePPG");

        gridLayout_5->addWidget(startNoisePPG, 1, 1, 1, 1);


        buttons->addLayout(gridLayout_5, 0, 1, 1, 1);

        gridLayout_2 = new QGridLayout();
        gridLayout_2->setObjectName("gridLayout_2");
        next8hours = new QPushButton(noise_marking_gui);
        next8hours->setObjectName("next8hours");
        next8hours->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(next8hours, 1, 0, 1, 1);

        prev8hours = new QPushButton(noise_marking_gui);
        prev8hours->setObjectName("prev8hours");
        prev8hours->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(prev8hours, 0, 0, 1, 1);

        undo_button = new QPushButton(noise_marking_gui);
        undo_button->setObjectName("undo_button");
        undo_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(undo_button, 0, 1, 1, 1);

        clearall_button = new QPushButton(noise_marking_gui);
        clearall_button->setObjectName("clearall_button");
        clearall_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(clearall_button, 1, 1, 1, 1);

        finalize_button = new QPushButton(noise_marking_gui);
        finalize_button->setObjectName("finalize_button");
        finalize_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(finalize_button, 1, 2, 1, 1);

        skip_button = new QPushButton(noise_marking_gui);
        skip_button->setObjectName("skip_button");
        skip_button->setMinimumSize(QSize(30, 0));
        skip_button->setMaximumSize(QSize(47, 16777215));

        gridLayout_2->addWidget(skip_button, 0, 2, 1, 1);


        buttons->addLayout(gridLayout_2, 0, 2, 1, 1);

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
        QFont font;
        font.setBold(true);
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


        buttons->addLayout(gridLayout, 3, 1, 1, 1);

        marking_type_label_2 = new QLabel(noise_marking_gui);
        marking_type_label_2->setObjectName("marking_type_label_2");
        marking_type_label_2->setFont(font);

        buttons->addWidget(marking_type_label_2, 2, 0, 1, 1);

        sleep_state_legend = new QGridLayout();
        sleep_state_legend->setSpacing(4);
        sleep_state_legend->setObjectName("sleep_state_legend");
        hboxLayout = new QHBoxLayout();
        hboxLayout->setObjectName("hboxLayout");
        rem_c = new QLabel(noise_marking_gui);
        rem_c->setObjectName("rem_c");
        rem_c->setMinimumSize(QSize(15, 15));
        rem_c->setMaximumSize(QSize(15, 15));

        hboxLayout->addWidget(rem_c);

        rem_l = new QLabel(noise_marking_gui);
        rem_l->setObjectName("rem_l");

        hboxLayout->addWidget(rem_l);


        sleep_state_legend->addLayout(hboxLayout, 0, 1, 1, 1);

        hboxLayout1 = new QHBoxLayout();
        hboxLayout1->setObjectName("hboxLayout1");
        n1_c = new QLabel(noise_marking_gui);
        n1_c->setObjectName("n1_c");
        n1_c->setMinimumSize(QSize(15, 15));
        n1_c->setMaximumSize(QSize(15, 15));

        hboxLayout1->addWidget(n1_c);

        n1_l = new QLabel(noise_marking_gui);
        n1_l->setObjectName("n1_l");

        hboxLayout1->addWidget(n1_l);


        sleep_state_legend->addLayout(hboxLayout1, 1, 0, 1, 1);

        hboxLayout2 = new QHBoxLayout();
        hboxLayout2->setObjectName("hboxLayout2");
        n3_c = new QLabel(noise_marking_gui);
        n3_c->setObjectName("n3_c");
        n3_c->setMinimumSize(QSize(15, 15));
        n3_c->setMaximumSize(QSize(15, 15));

        hboxLayout2->addWidget(n3_c);

        n3_l = new QLabel(noise_marking_gui);
        n3_l->setObjectName("n3_l");

        hboxLayout2->addWidget(n3_l);


        sleep_state_legend->addLayout(hboxLayout2, 2, 0, 1, 1);

        hboxLayout3 = new QHBoxLayout();
        hboxLayout3->setObjectName("hboxLayout3");
        n2_c = new QLabel(noise_marking_gui);
        n2_c->setObjectName("n2_c");
        n2_c->setMinimumSize(QSize(15, 15));
        n2_c->setMaximumSize(QSize(15, 15));

        hboxLayout3->addWidget(n2_c);

        n2_l = new QLabel(noise_marking_gui);
        n2_l->setObjectName("n2_l");

        hboxLayout3->addWidget(n2_l);


        sleep_state_legend->addLayout(hboxLayout3, 1, 1, 1, 1);

        hboxLayout4 = new QHBoxLayout();
        hboxLayout4->setObjectName("hboxLayout4");
        wake_c = new QLabel(noise_marking_gui);
        wake_c->setObjectName("wake_c");
        wake_c->setMinimumSize(QSize(15, 15));
        wake_c->setMaximumSize(QSize(15, 15));

        hboxLayout4->addWidget(wake_c);

        wake_l = new QLabel(noise_marking_gui);
        wake_l->setObjectName("wake_l");

        hboxLayout4->addWidget(wake_l);


        sleep_state_legend->addLayout(hboxLayout4, 0, 0, 1, 1);


        buttons->addLayout(sleep_state_legend, 3, 2, 1, 1);

        sleep_states_label = new QLabel(noise_marking_gui);
        sleep_states_label->setObjectName("sleep_states_label");
        sleep_states_label->setFont(font);

        buttons->addWidget(sleep_states_label, 2, 2, 1, 1);

        gridLayout_3 = new QGridLayout();
        gridLayout_3->setObjectName("gridLayout_3");
        window_length_label = new QLabel(noise_marking_gui);
        window_length_label->setObjectName("window_length_label");
        window_length_label->setFont(font);

        gridLayout_3->addWidget(window_length_label, 0, 0, 1, 1);

        window_length_selector = new QComboBox(noise_marking_gui);
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->addItem(QString());
        window_length_selector->setObjectName("window_length_selector");

        gridLayout_3->addWidget(window_length_selector, 1, 0, 1, 1);

        scatter_line = new QComboBox(noise_marking_gui);
        scatter_line->addItem(QString());
        scatter_line->addItem(QString());
        scatter_line->setObjectName("scatter_line");

        gridLayout_3->addWidget(scatter_line, 1, 1, 1, 1);

        label_4 = new QLabel(noise_marking_gui);
        label_4->setObjectName("label_4");
        label_4->setFont(font);

        gridLayout_3->addWidget(label_4, 0, 1, 1, 1);


        buttons->addLayout(gridLayout_3, 0, 0, 1, 1);

        marking_type_legend = new QGridLayout();
        marking_type_legend->setSpacing(4);
        marking_type_legend->setObjectName("marking_type_legend");
        hboxLayout5 = new QHBoxLayout();
        hboxLayout5->setObjectName("hboxLayout5");
        noise_color = new QLabel(noise_marking_gui);
        noise_color->setObjectName("noise_color");
        noise_color->setMinimumSize(QSize(15, 15));
        noise_color->setMaximumSize(QSize(15, 15));

        hboxLayout5->addWidget(noise_color);

        noise_desc = new QLabel(noise_marking_gui);
        noise_desc->setObjectName("noise_desc");

        hboxLayout5->addWidget(noise_desc);


        marking_type_legend->addLayout(hboxLayout5, 0, 0, 1, 1);

        hboxLayout6 = new QHBoxLayout();
        hboxLayout6->setObjectName("hboxLayout6");
        vt_color = new QLabel(noise_marking_gui);
        vt_color->setObjectName("vt_color");
        vt_color->setMinimumSize(QSize(15, 15));
        vt_color->setMaximumSize(QSize(15, 15));

        hboxLayout6->addWidget(vt_color);

        vt_desc = new QLabel(noise_marking_gui);
        vt_desc->setObjectName("vt_desc");

        hboxLayout6->addWidget(vt_desc);


        marking_type_legend->addLayout(hboxLayout6, 0, 1, 1, 1);

        hboxLayout7 = new QHBoxLayout();
        hboxLayout7->setObjectName("hboxLayout7");
        delay_color = new QLabel(noise_marking_gui);
        delay_color->setObjectName("delay_color");
        delay_color->setMinimumSize(QSize(15, 15));
        delay_color->setMaximumSize(QSize(15, 15));

        hboxLayout7->addWidget(delay_color);

        delay_desc = new QLabel(noise_marking_gui);
        delay_desc->setObjectName("delay_desc");

        hboxLayout7->addWidget(delay_desc);


        marking_type_legend->addLayout(hboxLayout7, 1, 0, 1, 1);

        hboxLayout8 = new QHBoxLayout();
        hboxLayout8->setObjectName("hboxLayout8");
        pac_color = new QLabel(noise_marking_gui);
        pac_color->setObjectName("pac_color");
        pac_color->setMinimumSize(QSize(15, 15));
        pac_color->setMaximumSize(QSize(15, 15));

        hboxLayout8->addWidget(pac_color);

        pac_desc = new QLabel(noise_marking_gui);
        pac_desc->setObjectName("pac_desc");

        hboxLayout8->addWidget(pac_desc);


        marking_type_legend->addLayout(hboxLayout8, 1, 1, 1, 1);

        hboxLayout9 = new QHBoxLayout();
        hboxLayout9->setObjectName("hboxLayout9");
        svt_color = new QLabel(noise_marking_gui);
        svt_color->setObjectName("svt_color");
        svt_color->setMinimumSize(QSize(15, 15));
        svt_color->setMaximumSize(QSize(15, 15));

        hboxLayout9->addWidget(svt_color);

        svt_desc = new QLabel(noise_marking_gui);
        svt_desc->setObjectName("svt_desc");

        hboxLayout9->addWidget(svt_desc);


        marking_type_legend->addLayout(hboxLayout9, 2, 0, 1, 1);

        hboxLayout10 = new QHBoxLayout();
        hboxLayout10->setObjectName("hboxLayout10");
        pvc_color = new QLabel(noise_marking_gui);
        pvc_color->setObjectName("pvc_color");
        pvc_color->setMinimumSize(QSize(15, 15));
        pvc_color->setMaximumSize(QSize(15, 15));

        hboxLayout10->addWidget(pvc_color);

        pvc_desc = new QLabel(noise_marking_gui);
        pvc_desc->setObjectName("pvc_desc");

        hboxLayout10->addWidget(pvc_desc);


        marking_type_legend->addLayout(hboxLayout10, 2, 1, 1, 1);

        hboxLayout11 = new QHBoxLayout();
        hboxLayout11->setObjectName("hboxLayout11");
        af_color = new QLabel(noise_marking_gui);
        af_color->setObjectName("af_color");
        af_color->setMinimumSize(QSize(15, 15));
        af_color->setMaximumSize(QSize(15, 15));

        hboxLayout11->addWidget(af_color);

        af_desc = new QLabel(noise_marking_gui);
        af_desc->setObjectName("af_desc");

        hboxLayout11->addWidget(af_desc);


        marking_type_legend->addLayout(hboxLayout11, 3, 0, 1, 1);

        hboxLayout12 = new QHBoxLayout();
        hboxLayout12->setObjectName("hboxLayout12");
        benign_color = new QLabel(noise_marking_gui);
        benign_color->setObjectName("benign_color");
        benign_color->setMinimumSize(QSize(15, 15));
        benign_color->setMaximumSize(QSize(15, 15));

        hboxLayout12->addWidget(benign_color);

        benign_desc = new QLabel(noise_marking_gui);
        benign_desc->setObjectName("benign_desc");

        hboxLayout12->addWidget(benign_desc);


        marking_type_legend->addLayout(hboxLayout12, 3, 1, 1, 1);

        hboxLayout13 = new QHBoxLayout();
        hboxLayout13->setObjectName("hboxLayout13");
        signif_color = new QLabel(noise_marking_gui);
        signif_color->setObjectName("signif_color");
        signif_color->setMinimumSize(QSize(15, 15));
        signif_color->setMaximumSize(QSize(15, 15));

        hboxLayout13->addWidget(signif_color);

        signif_desc = new QLabel(noise_marking_gui);
        signif_desc->setObjectName("signif_desc");

        hboxLayout13->addWidget(signif_desc);


        marking_type_legend->addLayout(hboxLayout13, 4, 0, 1, 2);


        buttons->addLayout(marking_type_legend, 3, 0, 1, 1);

        gridLayout_4 = new QGridLayout();
        gridLayout_4->setObjectName("gridLayout_4");
        process_button = new QPushButton(noise_marking_gui);
        process_button->setObjectName("process_button");

        gridLayout_4->addWidget(process_button, 0, 0, 1, 1);

        browse_file_button = new QPushButton(noise_marking_gui);
        browse_file_button->setObjectName("browse_file_button");

        gridLayout_4->addWidget(browse_file_button, 1, 0, 1, 1);

        save_current_plot = new QPushButton(noise_marking_gui);
        save_current_plot->setObjectName("save_current_plot");

        gridLayout_4->addWidget(save_current_plot, 0, 1, 1, 1);

        show_peaks_check = new QCheckBox(noise_marking_gui);
        show_peaks_check->setObjectName("show_peaks_check");

        gridLayout_4->addWidget(show_peaks_check, 1, 1, 1, 1);


        buttons->addLayout(gridLayout_4, 2, 1, 1, 1);


        buttons_and_ampograms->addLayout(buttons);

        ampogram_and_sleepstates = new QVBoxLayout();
        ampogram_and_sleepstates->setSpacing(0);
        ampogram_and_sleepstates->setObjectName("ampogram_and_sleepstates");
        ampogram_and_sleepstates->setContentsMargins(-1, -1, -1, 0);
        ecg_ampogram_axis = new QChartView(noise_marking_gui);
        ecg_ampogram_axis->setObjectName("ecg_ampogram_axis");
        QSizePolicy sizePolicy(QSizePolicy::Policy::Preferred, QSizePolicy::Policy::Expanding);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(1);
        sizePolicy.setHeightForWidth(ecg_ampogram_axis->sizePolicy().hasHeightForWidth());
        ecg_ampogram_axis->setSizePolicy(sizePolicy);

        ampogram_and_sleepstates->addWidget(ecg_ampogram_axis);

        ppg_ampogram_axis = new QChartView(noise_marking_gui);
        ppg_ampogram_axis->setObjectName("ppg_ampogram_axis");
        sizePolicy.setHeightForWidth(ppg_ampogram_axis->sizePolicy().hasHeightForWidth());
        ppg_ampogram_axis->setSizePolicy(sizePolicy);

        ampogram_and_sleepstates->addWidget(ppg_ampogram_axis);

        hyp_accel_resp_cvp_axis = new QChartView(noise_marking_gui);
        hyp_accel_resp_cvp_axis->setObjectName("hyp_accel_resp_cvp_axis");
        sizePolicy.setHeightForWidth(hyp_accel_resp_cvp_axis->sizePolicy().hasHeightForWidth());
        hyp_accel_resp_cvp_axis->setSizePolicy(sizePolicy);
        hyp_accel_resp_cvp_axis->setMinimumSize(QSize(0, 85));

        ampogram_and_sleepstates->addWidget(hyp_accel_resp_cvp_axis);

        ampogram_and_sleepstates->setStretch(2, 1);

        buttons_and_ampograms->addLayout(ampogram_and_sleepstates);

        buttons_and_ampograms->setStretch(1, 1);

        verticalLayout->addLayout(buttons_and_ampograms);

        main_plots = new QVBoxLayout();
        main_plots->setSpacing(0);
        main_plots->setObjectName("main_plots");
        ecg_axis_1 = new QChartView(noise_marking_gui);
        ecg_axis_1->setObjectName("ecg_axis_1");
        QSizePolicy sizePolicy1(QSizePolicy::Policy::Expanding, QSizePolicy::Policy::Expanding);
        sizePolicy1.setHorizontalStretch(1);
        sizePolicy1.setVerticalStretch(1);
        sizePolicy1.setHeightForWidth(ecg_axis_1->sizePolicy().hasHeightForWidth());
        ecg_axis_1->setSizePolicy(sizePolicy1);
        ecg_axis_1->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));
        ecg_1_check = new QCheckBox(ecg_axis_1);
        ecg_1_check->setObjectName("ecg_1_check");
        ecg_1_check->setGeometry(QRect(2, 7, 66, 18));
        ecg_1_check->setMinimumSize(QSize(0, 10));
        ecg_1_gain = new QDoubleSpinBox(ecg_axis_1);
        ecg_1_gain->setObjectName("ecg_1_gain");
        ecg_1_gain->setGeometry(QRect(71, 6, 66, 20));
        ecg_1_gain->setMinimumSize(QSize(0, 10));

        main_plots->addWidget(ecg_axis_1);

        ecg_axis_2 = new QChartView(noise_marking_gui);
        ecg_axis_2->setObjectName("ecg_axis_2");
        sizePolicy1.setHeightForWidth(ecg_axis_2->sizePolicy().hasHeightForWidth());
        ecg_axis_2->setSizePolicy(sizePolicy1);
        ecg_axis_2->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));
        ecg_2_gain = new QDoubleSpinBox(ecg_axis_2);
        ecg_2_gain->setObjectName("ecg_2_gain");
        ecg_2_gain->setGeometry(QRect(80, 5, 66, 20));
        ecg_2_gain->setMinimumSize(QSize(0, 10));
        ecg_2_check = new QCheckBox(ecg_axis_2);
        ecg_2_check->setObjectName("ecg_2_check");
        ecg_2_check->setGeometry(QRect(1, 6, 66, 18));
        ecg_2_check->setMinimumSize(QSize(0, 10));

        main_plots->addWidget(ecg_axis_2);

        ecg_axis_3 = new QChartView(noise_marking_gui);
        ecg_axis_3->setObjectName("ecg_axis_3");
        sizePolicy1.setHeightForWidth(ecg_axis_3->sizePolicy().hasHeightForWidth());
        ecg_axis_3->setSizePolicy(sizePolicy1);
        ecg_axis_3->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));
        ecg_3_check = new QCheckBox(ecg_axis_3);
        ecg_3_check->setObjectName("ecg_3_check");
        ecg_3_check->setGeometry(QRect(1, 6, 66, 18));
        ecg_3_check->setMinimumSize(QSize(0, 10));
        ecg_3_gain = new QDoubleSpinBox(ecg_axis_3);
        ecg_3_gain->setObjectName("ecg_3_gain");
        ecg_3_gain->setGeometry(QRect(80, 5, 66, 20));
        ecg_3_gain->setMinimumSize(QSize(0, 10));

        main_plots->addWidget(ecg_axis_3);

        ppg_axis = new QChartView(noise_marking_gui);
        ppg_axis->setObjectName("ppg_axis");
        sizePolicy1.setHeightForWidth(ppg_axis->sizePolicy().hasHeightForWidth());
        ppg_axis->setSizePolicy(sizePolicy1);
        ppg_axis->setStyleSheet(QString::fromUtf8("background-color: white; border: 1px solid black;"));
        ppg_gain = new QDoubleSpinBox(ppg_axis);
        ppg_gain->setObjectName("ppg_gain");
        ppg_gain->setGeometry(QRect(80, 5, 66, 20));
        ppg_gain->setMinimumSize(QSize(0, 10));
        ppg_check = new QCheckBox(ppg_axis);
        ppg_check->setObjectName("ppg_check");
        ppg_check->setGeometry(QRect(1, 6, 66, 18));
        ppg_check->setMinimumSize(QSize(0, 10));

        main_plots->addWidget(ppg_axis);

        accel_or_abp_axis = new QChartView(noise_marking_gui);
        accel_or_abp_axis->setObjectName("accel_or_abp_axis");
        sizePolicy1.setHeightForWidth(accel_or_abp_axis->sizePolicy().hasHeightForWidth());
        accel_or_abp_axis->setSizePolicy(sizePolicy1);
        abp_gain = new QDoubleSpinBox(accel_or_abp_axis);
        abp_gain->setObjectName("abp_gain");
        abp_gain->setGeometry(QRect(80, 1, 66, 20));
        abp_check = new QCheckBox(accel_or_abp_axis);
        abp_check->setObjectName("abp_check");
        abp_check->setGeometry(QRect(1, 2, 66, 18));

        main_plots->addWidget(accel_or_abp_axis);

        main_plots->setStretch(0, 1);
        main_plots->setStretch(2, 1);
        main_plots->setStretch(3, 1);
        main_plots->setStretch(4, 1);

        verticalLayout->addLayout(main_plots);

        verticalLayout->setStretch(0, 1);
        verticalLayout->setStretch(1, 4);

        retranslateUi(noise_marking_gui);

        QMetaObject::connectSlotsByName(noise_marking_gui);
    } // setupUi

    void retranslateUi(QDialog *noise_marking_gui)
    {
        noise_marking_gui->setWindowTitle(QCoreApplication::translate("noise_marking_gui", "PPG/ECG Noise Marker", nullptr));
        startNoiseABP->setText(QCoreApplication::translate("noise_marking_gui", "Mark aBP Start", nullptr));
        start_all_mark->setText(QCoreApplication::translate("noise_marking_gui", "Mark All Start", nullptr));
        stop_all_mark->setText(QCoreApplication::translate("noise_marking_gui", "Mark All End", nullptr));
        stopNoiseABP->setText(QCoreApplication::translate("noise_marking_gui", "Mark aBP End", nullptr));
        stopNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG End", nullptr));
        startNoisePPG->setText(QCoreApplication::translate("noise_marking_gui", "Mark PPG Start", nullptr));
        next8hours->setText(QCoreApplication::translate("noise_marking_gui", "Next 8h", nullptr));
        prev8hours->setText(QCoreApplication::translate("noise_marking_gui", "Prev 8h", nullptr));
        undo_button->setText(QCoreApplication::translate("noise_marking_gui", "Undo", nullptr));
        clearall_button->setText(QCoreApplication::translate("noise_marking_gui", "Clear", nullptr));
        finalize_button->setText(QCoreApplication::translate("noise_marking_gui", "Save", nullptr));
        skip_button->setText(QCoreApplication::translate("noise_marking_gui", "Skip", nullptr));
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
        marking_type->setItemText(0, QCoreApplication::translate("noise_marking_gui", "1) Noise/Artifact", nullptr));
        marking_type->setItemText(1, QCoreApplication::translate("noise_marking_gui", "2) Cond. Delay", nullptr));
        marking_type->setItemText(2, QCoreApplication::translate("noise_marking_gui", "3) AF", nullptr));
        marking_type->setItemText(3, QCoreApplication::translate("noise_marking_gui", "4) SVT", nullptr));
        marking_type->setItemText(4, QCoreApplication::translate("noise_marking_gui", "5) VT", nullptr));
        marking_type->setItemText(5, QCoreApplication::translate("noise_marking_gui", "6) PVC", nullptr));
        marking_type->setItemText(6, QCoreApplication::translate("noise_marking_gui", "7) PAC", nullptr));
        marking_type->setItemText(7, QCoreApplication::translate("noise_marking_gui", "8) Benign Arr.", nullptr));
        marking_type->setItemText(8, QCoreApplication::translate("noise_marking_gui", "9) Significant Arr.", nullptr));

        label_2->setText(QCoreApplication::translate("noise_marking_gui", "Mark ECG Channel Start", nullptr));
        marking_type_label_2->setText(QCoreApplication::translate("noise_marking_gui", "Marking Type Legend", nullptr));
        rem_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: green; border-radius:3px;", nullptr));
        rem_l->setText(QCoreApplication::translate("noise_marking_gui", "REM", nullptr));
        n1_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: blue; border-radius:3px;", nullptr));
        n1_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM1", nullptr));
        n3_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: red; border-radius:3px;", nullptr));
        n3_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM3", nullptr));
        n2_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: cyan; border-radius:3px;", nullptr));
        n2_l->setText(QCoreApplication::translate("noise_marking_gui", "NREM2", nullptr));
        wake_c->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: black; border-radius:3px;", nullptr));
        wake_l->setText(QCoreApplication::translate("noise_marking_gui", "Wake", nullptr));
        sleep_states_label->setText(QCoreApplication::translate("noise_marking_gui", "Sleep States Legend", nullptr));
        window_length_label->setText(QCoreApplication::translate("noise_marking_gui", "Window Length:", nullptr));
        window_length_selector->setItemText(0, QCoreApplication::translate("noise_marking_gui", "1 Second", nullptr));
        window_length_selector->setItemText(1, QCoreApplication::translate("noise_marking_gui", "3 Second", nullptr));
        window_length_selector->setItemText(2, QCoreApplication::translate("noise_marking_gui", "10 Seconds", nullptr));
        window_length_selector->setItemText(3, QCoreApplication::translate("noise_marking_gui", "30 Seconds", nullptr));
        window_length_selector->setItemText(4, QCoreApplication::translate("noise_marking_gui", "1 Minute", nullptr));
        window_length_selector->setItemText(5, QCoreApplication::translate("noise_marking_gui", "2 Minutes", nullptr));
        window_length_selector->setItemText(6, QCoreApplication::translate("noise_marking_gui", "5 Minutes", nullptr));

        scatter_line->setItemText(0, QCoreApplication::translate("noise_marking_gui", "Line Plot", nullptr));
        scatter_line->setItemText(1, QCoreApplication::translate("noise_marking_gui", "Scatter Plot", nullptr));

        label_4->setText(QCoreApplication::translate("noise_marking_gui", "Plotting Style", nullptr));
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
        process_button->setText(QCoreApplication::translate("noise_marking_gui", "Process Output", nullptr));
        browse_file_button->setText(QCoreApplication::translate("noise_marking_gui", "Browse Different File", nullptr));
        save_current_plot->setText(QCoreApplication::translate("noise_marking_gui", "Save Current Plot", nullptr));
        show_peaks_check->setText(QCoreApplication::translate("noise_marking_gui", "Show Peaks?", nullptr));
        ecg_ampogram_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        ppg_ampogram_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        hyp_accel_resp_cvp_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid gray;", nullptr));
        ecg_1_check->setText(QCoreApplication::translate("noise_marking_gui", "Fix Scale", nullptr));
        ecg_2_check->setText(QCoreApplication::translate("noise_marking_gui", "Fix Scale", nullptr));
        ecg_3_check->setText(QCoreApplication::translate("noise_marking_gui", "Fix Scale", nullptr));
        ppg_check->setText(QCoreApplication::translate("noise_marking_gui", "Fix Scale", nullptr));
        accel_or_abp_axis->setStyleSheet(QCoreApplication::translate("noise_marking_gui", "background-color: white; border: 1px solid black;", nullptr));
        abp_check->setText(QCoreApplication::translate("noise_marking_gui", "Fix Scale", nullptr));
    } // retranslateUi

};

namespace Ui {
    class noise_marking_gui: public Ui_noise_marking_gui {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_NOISE_MARKING_GUI_H
