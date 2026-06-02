/**
 * @file   post_process_queue.cpp
 * @brief  Implementation of PostProcessQueue. See header for usage.
 */
#include "post_process_queue.hpp"

#include <iostream>
#include "post_process.hpp"

PostProcessQueue::~PostProcessQueue() { drain(); }

void PostProcessQueue::enqueue(const std::filesystem::path& binPath,
    const config_entry& cfg) {
    const std::string key = binPath.string();
    {
        std::lock_guard<std::mutex> lk(m_mu);
        m_locked.insert(key);
        m_queue.push_back({ binPath, cfg });
    }
    m_cv.notify_one();

    // Lazy worker start. Once started, the thread sleeps on the cv
    // when idle; we don't churn threads per job.
    if (!m_worker.joinable()) {
        m_worker = std::thread([this] { run(); });
    }
}

bool PostProcessQueue::isLocked(const std::filesystem::path& binPath) const {
    std::lock_guard<std::mutex> lk(m_mu);
    return m_locked.count(binPath.string()) > 0;
}

int PostProcessQueue::pendingCount() const {
    std::lock_guard<std::mutex> lk(m_mu);
    return m_locked.size();
}

void PostProcessQueue::drain() {
    {
        std::lock_guard<std::mutex> lk(m_mu);
        m_stop = true;
    }
    m_cv.notify_one();
    if (m_worker.joinable()) m_worker.join();
}

void PostProcessQueue::run() {
    for (;;) {
        Job job;
        {
            std::unique_lock<std::mutex> lk(m_mu);
            m_cv.wait(lk, [this] { return m_stop || !m_queue.empty(); });
            if (m_queue.empty()) {
                if (m_stop) return;
                continue;
            }
            job = std::move(m_queue.front());
            m_queue.pop_front();
            // Note: keep the path in m_locked until the job FINISHES,
            // not just dequeues. Removing on dequeue would unlock the
            // file while processing was still mid-flight.
        }
        try {
            post_process_detail::processOneFile(job.cfg, job.binPath);
        }
        catch (const std::exception& e) {
            std::cerr << "  background error: " << e.what() << "\n";
        }
        catch (...) {
            std::cerr << "  background error (unknown).\n";
        }

        {
            std::lock_guard<std::mutex> lk(m_mu);
            m_locked.erase(job.binPath.string());
        }
    }
}