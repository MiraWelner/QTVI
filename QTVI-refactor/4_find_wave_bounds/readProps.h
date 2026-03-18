// ============================================================================
// File: readProps.h
// Read key=value properties from a configuration file
// ============================================================================
#pragma once

#include "SignalProcessingTypes.h"

map<string, string> readProps(const string& inputFile);