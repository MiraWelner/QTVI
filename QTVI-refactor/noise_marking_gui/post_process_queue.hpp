/**
 * @file   post_process_queue.hpp
 * @brief  Single-worker FIFO queue that runs post-processing in the
 *         background, so the GUI can move on to marking the next file
 *         while the previous one is being annealed / peak-found /
 *         templated. Files are locked while queued or in flight;
 *         isLocked() is the way to check from the GUI.
 *
 * @author Mira Welner
 * @email  MEW386@pitt.edu
 */
#pragma once

#include <condition_variable>
#include <deque>
#include <filesystem>
#include <mutex>
#include <set>
#include <string>
#include <thread>

#include "config_entry.hpp"

class PostProcessQueue {
public:
    PostProcessQueue() = default;
    ~PostProcessQueue();

    PostProcessQueue(const PostProcessQueue&) = delete;
    PostProcessQueue& operator=(const PostProcessQueue&) = delete;

    /// Queue a file for background post-processing. The file is locked
    /// immediately and stays locked until processOneFile() completes.
    void enqueue(const std::filesystem::path& binPath, const config_entry& cfg);

    /// True if `binPath` is currently queued or being processed.
    bool isLocked(const std::filesystem::path& binPath) const;

    /// Jobs currently queued or in flight. For status display only.
    int pendingCount() const;

    /// Block until all queued jobs finish and the worker thread exits.
    /// Called automatically by the destructor; call manually at program
    /// shutdown if you want explicit control over the wait.
    void drain();

private:
    struct Job {
        std::filesystem::path binPath;
        config_entry          cfg;
    };

    void run();

    mutable std::mutex      m_mu;
    std::condition_variable m_cv;
    std::deque<Job>         m_queue;
    std::set<std::string>   m_locked;
    std::thread             m_worker;
    bool                    m_stop = false;
};