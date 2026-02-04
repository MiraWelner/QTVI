// lower_row_buttons.h
#pragma once

#include <QObject>

// Forward declarations
class noise_marking_gui;

class lower_row_buttons : public QObject {
    Q_OBJECT

public:
    explicit lower_row_buttons(noise_marking_gui* parent);
    void setupConnections();
    ~lower_row_buttons() = default;

    // Button handlers
    void handle_undo_button();
    void handle_clearall_button();
    void handle_finalize_button();
    void handle_skip_button();
    void handle_ecgmarkingstart_button();
    void handle_ecgmarkingstop_button();
    void handle_ppgmarkingstart_button();
    void handle_ppgmarkingstop_button();

    void handle_10s_window_toggle(bool checked);
    void handle_30s_window_toggle(bool checked);
    void handle_1m_window_toggle(bool checked);
    void handle_5m_window_toggle(bool checked);
    void handle_10m_window_toggle(bool checked);

private:
    noise_marking_gui* m_gui;
};