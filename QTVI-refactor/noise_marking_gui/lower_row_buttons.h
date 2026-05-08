/**
 * @file   lower_row_buttons.h
 * @brief  Handles the bottom toolbar: undo, clear, skip, save, marking
 *         start/stop for each ECG channel and PPG, mark-all, and
 *         window-size radio buttons.
 */
#pragma once

#include <QObject>

class noise_marking_gui;

class lower_row_buttons : public QObject {
    Q_OBJECT

public:
    explicit lower_row_buttons(noise_marking_gui* parent);
    ~lower_row_buttons() = default;
    void setupConnections();

private:
    noise_marking_gui* m_gui;

    void handle_undo_button();
    void handle_clearall_button();
    void handle_finalize_button();
    void handle_skip_button();

    void handle_allmarkingstart_button();
    void handle_allmarkingstop_button();

    void handle_ppgmarkingstart_button();
    void handle_ppgmarkingstop_button();

    void handle_window_toggle(bool checked, double duration);
};