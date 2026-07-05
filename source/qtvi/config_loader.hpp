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

bool load_config(int dataType, config_entry& out);
