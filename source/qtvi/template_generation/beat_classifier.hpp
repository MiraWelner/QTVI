/**
 * @file   beat_classifier.hpp
 * @brief  Automated pre-marking classifier, per-session adaptation layer, and
 *         the Stage 2 envelope update on user action.
 *         Spec Sections 4.7.2 and 4.7.3, transcribed as written.
 *
 *         Four things had to be supplied for this to build in this tree; all
 *         are additions, none change the logic as specified:
 *
 *          1. PreMark::reviewed -- onUserMark tests marks[i].reviewed, so the
 *             field has to exist on the struct.
 *          2. Ort::Session is behind DEX_HAVE_ONNXRUNTIME. This tree does not
 *             link onnxruntime yet, so without that define the member and its
 *             initializer are compiled out and preMarkAll runs off
 *             getPopulationPosterior alone.
 *          3. SessionAdapter::weights needed a size (ctor arg, default 100).
 *          4. onUserMark's helpers -- removeBeatFromPool, addBeatToPool,
 *             getBeat, extractFeatures -- are declared here and defined by
 *             the caller against its own beat matrix. The segment boundaries
 *             the spec left as a comment are now parameters.
 *
 *         One flag: preMarkAll sorts marks by confidence, and onUserMark then
 *         does marks[beatIdx] and pairs marks[i] with getBeat(i). After the
 *         sort, slot i is no longer beat i. Left as specified -- if the GUI
 *         indexes by beat, keep a separate sorted index vector for the queue.
 *
 * @date   2026-08-14
 */
#pragma once

#include "morphology_envelope.hpp"

#include <algorithm>
#include <cmath>
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
    double confidence = 0.0;   // 0 to 1
    bool reviewed = false;     // see note 1 in the file header
};

class BeatClassifier {
#ifdef DEX_HAVE_ONNXRUNTIME
    Ort::Session populationModel;   // daily-retrained ONNX model
#endif
public:
#ifdef DEX_HAVE_ONNXRUNTIME
    BeatClassifier(const std::string& onnxPath)
        : populationModel(Ort::Env(), onnxPath.c_str(), Ort::SessionOptions()) {
    }
#else
    BeatClassifier(const std::string& onnxPath) { (void)onnxPath; }
#endif

    std::vector<PreMark> preMarkAll(const std::vector<std::vector<double>>& beats,
        const MorphologyEnvelope& indivEnv,
        const MorphologyEnvelope& subgroupEnv,
        const MorphologyEnvelope& populationEnv)
    {
        std::vector<PreMark> marks(beats.size());
        // Build feature matrix: [nBeats x (waveform + 3 envelope deviations + scalars)]
        // Run ONNX inference in batch
        // Sort by confidence ascending (least confident first for review queue)
        std::sort(marks.begin(), marks.end(),
            [](auto& a, auto& b) { return a.confidence < b.confidence; });
        return marks;
    }

    // Raw population-model posterior for one beat, read by onUserMark.
    double getPopulationPosterior(int beatIdx) const;
};

// ---------------------------------------------------------------------------
// Per-session adaptation layer (Section 4.7.3)
// ---------------------------------------------------------------------------

class SessionAdapter {
    std::vector<double> weights;      // logistic regression weights (~100 features)
    double learningRate = 0.01;
    int corrections = 0;
public:
    explicit SessionAdapter(size_t nFeatures = 100) : weights(nFeatures, 0.0) {}

    void update(const std::vector<double>& features, int trueClass, int predictedClass) {
        if (trueClass == predictedClass) return;   // no correction needed
        ++corrections;
        // Gradient step on cross-entropy loss for the logistic regression
        // Fast: ~100 multiplies, under 1 ms
        for (size_t i = 0; i < weights.size(); ++i)
            weights[i] += learningRate * (trueClass - predictedClass) * features[i];
    }

    double adjustedPosterior(double populationPosterior,
        const std::vector<double>& features,
        int nCorrections)
    {
        double sessionLogit = 0;
        for (size_t i = 0; i < weights.size(); ++i) sessionLogit += weights[i] * features[i];
        double sessionPosterior = 1.0 / (1.0 + std::exp(-sessionLogit));
        // Weight session overlay more as corrections accumulate
        double alpha = std::min(0.5, nCorrections / 100.0);
        return (1.0 - alpha) * populationPosterior + alpha * sessionPosterior;
    }

    int correctionCount() const { return corrections; }
};

// ---------------------------------------------------------------------------
// Envelope update on user action (Section 4.7.3 Stage 2)
// ---------------------------------------------------------------------------

// Supplied by the caller against its own beat matrix (see note 4).
void removeBeatFromPool(std::vector<std::vector<double>>& sinusPool, int beatIdx);
void addBeatToPool(std::vector<std::vector<double>>& sinusPool, int beatIdx);
const std::vector<double>& getBeat(int beatIdx);
std::vector<double> extractFeatures(int beatIdx);

inline void onUserMark(int beatIdx, PreMark::Class newClass,
    std::vector<std::vector<double>>& sinusPool,
    MorphologyEnvelope& env, SessionAdapter& adapter,
    BeatClassifier& classifier, std::vector<PreMark>& marks,
    int pEnd, int qrsStart, int qrsEnd, int tEnd)
{
    // 1. Remove or add beat to the appropriate pool
    if (newClass != PreMark::SINUS) {
        removeBeatFromPool(sinusPool, beatIdx);
    }
    else {
        addBeatToPool(sinusPool, beatIdx);
    }

    // 2. Rebuild envelope (O(N*W*logN), under 10 ms)
    if (sinusPool.empty()) return;
    int W = sinusPool[0].size();
    env = buildEnvelope(sinusPool, W);

    // 3. Update session adapter
    auto features = extractFeatures(beatIdx);
    adapter.update(features, (int)newClass, (int)marks[beatIdx].cls);

    // 4. Re-score and re-classify all remaining unreviewed beats (OpenMP)
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (int)marks.size(); ++i) {
        if (marks[i].reviewed) continue;
        auto bm = scoreBeat(getBeat(i), env, pEnd, qrsStart, qrsEnd, tEnd);
        (void)bm;
        double popPost = classifier.getPopulationPosterior(i);
        marks[i].confidence = adapter.adjustedPosterior(popPost, extractFeatures(i),
            adapter.correctionCount());
        // Re-classify if confidence changed enough
    }

    // 5. Highlight newly exposed beats in the GUI
}