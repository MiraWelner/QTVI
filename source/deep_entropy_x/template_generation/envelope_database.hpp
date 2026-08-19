/**
 * @file   envelope_database.hpp
 * @brief  Three-level envelope database: individual / subgroup / population.
 *         Spec Section 4.7.4, transcribed as written.
 *
 *         The PostgreSQL / TimescaleDB / FAISS backends named in the spec are
 *         not linked in this tree, so the methods are declared here and left
 *         for the storage layer to define. Nothing above this header changes
 *         when they land.
 *
 * @date   2026-08-14
 */
#pragma once

#include "beat_classifier.hpp"      // PreMark::Class
#include "morphology_envelope.hpp"

#include <cstdint>
#include <map>
#include <string>
#include <vector>

struct EnvelopeRecord {
    std::string patientId;
    int64_t recordingTimestamp;
    MorphologyEnvelope sinusEnv;
    std::map<PreMark::Class, MorphologyEnvelope> arrhythmiaEnvs;
    int totalBeats, reviewedBeats, autoClassifiedBeats;
};

class EnvelopeDatabase {
    // PostgreSQL: envelope metadata
    // TimescaleDB: per-sample bounds as time-series
    // FAISS: envelope shape vectors for subgroup clustering
public:
    void storeIndividual(const EnvelopeRecord& rec);
    MorphologyEnvelope getIndividual(const std::string& patientId);
    MorphologyEnvelope getSubgroup(const std::vector<float>& profileVector, int k = 10);
    MorphologyEnvelope getPopulation();
    std::vector<EnvelopeRecord> getHistory(const std::string& patientId);
};