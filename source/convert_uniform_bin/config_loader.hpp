#pragma once
/**
 * @file   bin_handler.hpp
 * @brief  Source-file discovery and config loading for the QTVi pipeline.
 *         The `ConfigEntry` struct lives in config_entry.hpp (Qt-free)
 *         so file_to_bin can include it without pulling in Qt.
 */

#include <QString>
#include <QStringList>
#include <filesystem>

#include "config_entry.hpp"

 /// Load a single dataset row from config.csv. dataType: 1=MESA, 2=Bittium,
 /// 3=CHAOS. Returns false if the file can't be opened or the dataset row
 /// isn't found.
bool load_config(int dataType, config_entry& out);

/// Recursively list every source file (matching cfg.dataset_type) under
/// cfg.originalFilePath. Result is sorted alphabetically.
QStringList discoverSourceFiles(const config_entry& cfg);

/**
* @brief if the paths are not in the config file, prompt for them
*/
bool promptForMissingPaths(config_entry& cfg);