/**
 * @file   user_control_handler.h
 * @brief  Event handlers for the noise-marking GUI,
 *         clear, finalize, skip, save, per-channel marking start/stop,
 *         mark-all, and the window-length selector.
 *
 *         Owns no data of its own; every action forwards to the parent
 *         noise_marking_gui via m_gui.
 */
#pragma once

#include <QObject>

class noise_marking_gui;

class user_control_handler : public QObject {
    Q_OBJECT

public:
    explicit user_control_handler(noise_marking_gui* parent);
    ~user_control_handler() = default;

    /// Wire up every button on the bottom toolbar to its handler.
    void setupConnections();

private:
    // --- Button Presses
    void handle_clearall_button();
    void handle_finalize_button();
    void handle_skip_button();
    void save_current_plot();
    void save_current_csv();

    // --- Window-length selector ---------------------------------------------
    void handle_window_toggle(bool checked, double duration);

    noise_marking_gui* m_gui;
};