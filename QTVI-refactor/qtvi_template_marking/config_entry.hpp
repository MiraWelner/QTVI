#pragma once
/**
 * @file   config_entry.hpp
 * @brief  Fields the template-marking GUI needs at runtime.
 *
 *         The CSV on disk (config.csv) has more columns than this --
 *         most belong to sibling programs in the pipeline. The loader
 *         reads only the rate / bin-length cells and ignores the rest.
 *
 *         Paths:
 *           - output_path     : root chosen by the user at startup
 *                               (the noise-marking GUI's output folder).
 *           - template_path   : where *_templates.bin files are read from.
 *           - qtvi_marker_path: where *_template_markings.bin files are
 *                               written. Created inside output_path as
 *                               "QTVI_markings/".
 */

#include <string>

struct config_entry {
    // CSV-driven filters used to skip template files whose filename-encoded
    // rate or bin length don't match the dataset's expected values.
    double finalSamplingRate = 0.0;
    double bin_length_minutes = 0.0;

    // Folder chosen by the user; the two paths below are derived from it.
    std::string output_path;
    std::string template_path;
    std::string qtvi_marker_path;
};