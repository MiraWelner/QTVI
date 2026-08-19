/**
 * @file   beat_classifier.hpp
 * @brief  Automated pre-marking classifier (4.7.2) and the Stage 2 envelope
 *         update on user action (4.7.3).
 *
 *         PHASE SCOPE. The ONNX population model, the per-session adapter and
 *         the daily retraining are later-phase work. In this phase the
 *         cold-start behaviour is the correct behaviour: preMarkAll returns
 *         cls UNKNOWN with confidence 0 for every beat. What is live now is
 *         the envelope: buildEnvelope, scoreBeat, and the rebuild-and-rescore
 *         cycle when a beat is reclassified.
 *
 *         Four corrections against the literal transcription:
 *
 *          1. IT LINKS. removeBeatFromPool, addBeatToPool, getBeat and
 *             extractFeatures were declared here and defined nowhere in the
 *             tree, so any translation unit calling onUserMark failed at link.
 *             The pool is now explicit state (EnvelopeSession) instead of free
 *             functions over a hidden matrix, and extractFeatures --
 *             later-phase, feature vector for the model -- is not referenced.
 *
 *          2. INDEX ALIASING. preMarkAll sorted marks by confidence while
 *             onUserMark did marks[beatIdx] and paired marks[i] with beat i.
 *             After the sort slot i was not beat i. The sort is gone;
 *             reviewOrder() returns the least-confident-first ordering as a
 *             separate index vector, which is what a review queue needs
 *             anyway. beatIndex is now populated, so the mapping is
 *             recoverable either way.
 *
 *          3. POOL MEMBERSHIP BY BEAT INDEX. removeBeatFromPool(pool, beatIdx)
 *             indexed a compacted pool with a global beat index; the two
 *             diverge as soon as one beat leaves, and a removed beat could
 *             never be re-added because its row was gone. Membership is now a
 *             flag over the full beat list.
 *
 *          4. THE RE-SCORE FEEDS SOMETHING. Step 4 computed scoreBeat and
 *             discarded it via (void)bm, taking confidence from the ONNX
 *             posterior alone -- so in cold start the whole loop was a no-op.
 *             The band scores are now retained, and reviewOrder falls back to
 *             ranking on them while confidence is uniformly 0.
 *
 * @date   2026-08-19
 */
#pragma once

#include "morphology_envelope.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef DEX_HAVE_ONNXRUNTIME
#include <onnxruntime_cxx_api.h>
#endif

 // ---------------------------------------------------------------------------
 // Automated pre-marking classifier (Section 4.7.2)
 // ---------------------------------------------------------------------------

struct PreMark {
    int beatIndex = -1;
    enum Class { SINUS, PVC, PAC, ARTIFACT, CONDUCTION, SVT, VT, AF, UNKNOWN } cls = UNKNOWN;
    double confidence = 0.0;   // 0 to 1; stays 0 until the model is wired
    bool reviewed = false;     // an operator verdict is final
};

inline const char* preMarkClassName(PreMark::Class c) {
    switch (c) {
    case PreMark::SINUS:      return "sinus";
    case PreMark::PVC:        return "pvc";
    case PreMark::PAC:        return "pac";
    case PreMark::ARTIFACT:   return "artifact";
    case PreMark::CONDUCTION: return "conduction";
    case PreMark::SVT:        return "svt";
    case PreMark::VT:         return "vt";
    case PreMark::AF:         return "af";
    default:                  return "unknown";
    }
}

class BeatClassifier {
#ifdef DEX_HAVE_ONNXRUNTIME
    Ort::Session populationModel;   // daily-retrained ONNX model (later phase)
#endif
public:
#ifdef DEX_HAVE_ONNXRUNTIME
    explicit BeatClassifier(const std::string& onnxPath)
        : populationModel(Ort::Env(), onnxPath.c_str(), Ort::SessionOptions()) {
    }
#else
    explicit BeatClassifier(const std::string& onnxPath = {}) { (void)onnxPath; }
#endif

    // Cold start: one mark per beat, in beat order, cls UNKNOWN, confidence 0,
    // beatIndex populated. The envelope arguments are accepted so the
    // signature does not change when the model lands; they are unused now.
    std::vector<PreMark> preMarkAll(const std::vector<std::vector<double>>& beats,
        const MorphologyEnvelope& indivEnv,
        const MorphologyEnvelope& subgroupEnv,
        const MorphologyEnvelope& populationEnv)
    {
        (void)indivEnv; (void)subgroupEnv; (void)populationEnv;
        std::vector<PreMark> marks(beats.size());
        for (size_t i = 0; i < marks.size(); ++i)
            marks[i].beatIndex = static_cast<int>(i);
        return marks;
    }

    // Raw population-model posterior. Returns NaN without a model, so callers
    // can test for it rather than consuming a fabricated 0.
    double getPopulationPosterior(int beatIdx) const {
        (void)beatIdx;
        return std::numeric_limits<double>::quiet_NaN();
    }
};

// ---------------------------------------------------------------------------
// Per-session adaptation layer (Section 4.7.3) -- LATER PHASE
// ---------------------------------------------------------------------------
// Compiles, and is referenced by nothing above. It cannot contribute before a
// population posterior exists to overlay: with all-zero weights
// sessionPosterior is exactly 0.5, so adjustedPosterior would pull every
// confidence toward 0.5 on no evidence.

class SessionAdapter {
    std::vector<double> weights;      // logistic regression weights (~100 features)
    double learningRate = 0.01;
    int corrections = 0;
public:
    explicit SessionAdapter(size_t nFeatures = 100) : weights(nFeatures, 0.0) {}

    void update(const std::vector<double>& features, int trueClass, int predictedClass) {
        if (trueClass == predictedClass) return;   // no correction needed
        ++corrections;
        const size_t n = std::min(weights.size(), features.size());
        for (size_t i = 0; i < n; ++i)
            weights[i] += learningRate * (trueClass - predictedClass) * features[i];
    }

    double adjustedPosterior(double populationPosterior,
        const std::vector<double>& features,
        int nCorrections) const
    {
        double sessionLogit = 0;
        const size_t n = std::min(weights.size(), features.size());
        for (size_t i = 0; i < n; ++i) sessionLogit += weights[i] * features[i];
        const double sessionPosterior = 1.0 / (1.0 + std::exp(-sessionLogit));
        const double alpha = std::min(0.5, nCorrections / 100.0);
        return (1.0 - alpha) * populationPosterior + alpha * sessionPosterior;
    }

    int correctionCount() const { return corrections; }
};

// ---------------------------------------------------------------------------
// Envelope update on user action (Section 4.7.3 Stage 2)
// ---------------------------------------------------------------------------

// The state Stage 2 operates on: one bin's beats, which of them currently
// contribute to the corridor, the corridor, and the per-beat band scores
// against it. `beats` is borrowed -- the caller owns the matrix (typically
// out.beats.per_channel_beats[ch][bin]).
struct EnvelopeSession {
    const std::vector<std::vector<double>>* beats = nullptr;
    std::vector<char>            inPool;    // 1 = contributes to the envelope
    std::vector<PreMark>         marks;     // beat order, parallel to *beats
    std::vector<BandMatchResult> scores;    // beat order, against env
    MorphologyEnvelope           env;
    int pEnd = -1, qrsStart = -1, qrsEnd = -1, tBegin = -1, tEnd = -1;
    int W = 0;
    int minPoolBeats = 20;                  // refuse to shrink below this

    int nBeats() const { return beats ? static_cast<int>(beats->size()) : 0; }
    int poolSize() const {
        return static_cast<int>(std::count(inPool.begin(), inPool.end(), static_cast<char>(1)));
    }
    // Beats not yet ruled on -- the convergence readout.
    int unreviewedCount() const {
        int n = 0;
        for (const auto& m : marks) if (!m.reviewed) ++n;
        return n;
    }
};

// Copy out the pool rows. buildEnvelope takes a vector of vectors; at ~400
// doubles per beat and a few thousand beats per bin this is a few MB, inside
// the 10 ms the spec budgets for a rebuild.
inline std::vector<std::vector<double>> poolBeats(const EnvelopeSession& s) {
    std::vector<std::vector<double>> pool;
    if (!s.beats) return pool;
    pool.reserve(static_cast<size_t>(s.poolSize()));
    for (size_t i = 0; i < s.beats->size() && i < s.inPool.size(); ++i)
        if (s.inPool[i]) pool.push_back((*s.beats)[i]);
    return pool;
}

inline void rebuildEnvelope(EnvelopeSession& s) {
    const auto pool = poolBeats(s);
    if (pool.empty()) return;                    // keep the last good corridor
    s.env = buildEnvelope(pool, s.W);
}

// Re-score every beat not yet ruled on. Reviewed beats keep the scores they
// held when reviewed, so a displayed number cannot shift under a verdict
// already given.
inline void rescoreUnreviewed(EnvelopeSession& s) {
    if (!s.beats) return;
    const int n = s.nBeats();
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < n; ++i) {
        if (s.marks[i].reviewed) continue;
        s.scores[i] = scoreBeat((*s.beats)[i], s.env, s.pEnd, s.qrsStart,
            s.qrsEnd, s.tBegin, s.tEnd);
        // confidence is deliberately NOT touched: there is no model and no
        // adapter in this phase, so it stays 0 rather than being synthesised
        // from the band percentages.
    }
}

// Initial state. Every beat starts in the pool -- with no classifier output
// there is no basis for excluding one; operator marks are what take beats out.
inline EnvelopeSession makeSession(const std::vector<std::vector<double>>& beats,
    int pEnd, int qrsStart, int qrsEnd, int tBegin, int tEnd,
    BeatClassifier& classifier, int minPoolBeats = 20)
{
    EnvelopeSession s;
    if (beats.empty() || beats[0].empty()) return s;
    s.beats = &beats;
    s.W = static_cast<int>(beats[0].size());
    s.pEnd = pEnd; s.qrsStart = qrsStart; s.qrsEnd = qrsEnd;
    s.tBegin = tBegin; s.tEnd = tEnd;
    s.minPoolBeats = minPoolBeats;
    s.inPool.assign(beats.size(), 1);
    s.scores.assign(beats.size(), BandMatchResult{});
    rebuildEnvelope(s);
    s.marks = classifier.preMarkAll(beats, s.env, MorphologyEnvelope{},
        MorphologyEnvelope{});
    rescoreUnreviewed(s);
    return s;
}

// Steps 1, 2 and 4 of Stage 2. Step 3 (session adapter) is later-phase; step 5
// (GUI highlight) belongs to the caller.
//
// Returns true if the corridor was rebuilt. Idempotent: re-marking a beat with
// the class it already holds records the review but moves no membership, so
// the envelope does not drift under repeated clicks.
inline bool onUserMark(EnvelopeSession& s, int beatIdx, PreMark::Class newClass)
{
    if (!s.beats || beatIdx < 0 || beatIdx >= s.nBeats()) return false;

    // 1. Membership follows the verdict: SINUS contributes, nothing else does.
    const char want = (newClass == PreMark::SINUS) ? 1 : 0;
    const bool membershipChanged = (s.inPool[beatIdx] != want);

    // The operator's verdict: certain, and final.
    s.marks[beatIdx].cls = newClass;
    s.marks[beatIdx].confidence = 1.0;
    s.marks[beatIdx].reviewed = true;

    if (!membershipChanged) return false;

    // Refuse to shrink the pool past the point where a percentile corridor
    // means anything. Without this the operator can reach a state where the
    // corridor rejects everything and there is no way back.
    if (want == 0 && s.poolSize() - 1 < s.minPoolBeats) return false;
    s.inPool[beatIdx] = want;

    // 2. Rebuild, then 4. re-score the beats still awaiting review.
    rebuildEnvelope(s);
    rescoreUnreviewed(s);
    return true;
}

// Review order, least confident first. Returned as indices INTO marks, so
// marks stays in beat order. While every confidence is 0 (cold start) the
// ordering falls back to the weakest QRS band match, which is the informative
// quantity in that regime -- see the dilution note in premark_beats.hpp.
inline std::vector<int> reviewOrder(const EnvelopeSession& s) {
    std::vector<int> order(s.marks.size());
    std::iota(order.begin(), order.end(), 0);

    bool anyConfidence = false;
    for (const auto& m : s.marks)
        if (m.confidence > 0.0 && !m.reviewed) { anyConfidence = true; break; }

    if (anyConfidence) {
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            return s.marks[a].confidence < s.marks[b].confidence;
            });
    }
    else {
        // NaN sorts last: an unscorable band is not evidence of a bad beat.
        std::stable_sort(order.begin(), order.end(), [&](int a, int b) {
            const double qa = s.scores[a].pct_QRS, qb = s.scores[b].pct_QRS;
            if (std::isnan(qa)) return false;
            if (std::isnan(qb)) return true;
            return qa < qb;
            });
    }
    return order;
}

// Mean per-sample sd over defined columns -- the tightening readout ("SD
// decreases after PVC removal"). Call before and after an onUserMark that
// removed an ectopic beat.
inline double envelopeMeanSd(const MorphologyEnvelope& env) {
    double s = 0.0; int n = 0;
    for (size_t c = 0; c < env.sd.size(); ++c)
        if (!std::isnan(env.sd[c]) && env.sd[c] > 0.0) { s += env.sd[c]; ++n; }
    return n ? s / n : std::numeric_limits<double>::quiet_NaN();
}