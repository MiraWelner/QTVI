#pragma once
/**
 * @file   config_loader.cpp
 * @brief  Loads the config.csv file, parses it based on the datset type selected by the user, and fills up a config_entry struct with the relevant paths, rates, and channel labels.
 *         The channel labels (eg. "ECG_1" vs "EKG") are dataset-specific but not in the config file, so they are assigned in apply_dataset_specific_channel_labels() based on the dataset type.
 *         The output paths are either found in the config file, or prompted for manually if the config file cells are blank. The output_path is used to create the subfolders where the
 *         specific types of out put are found.
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