// ============================================================================
// File: readProps.h
// Read properties from configuration file
// ============================================================================
#ifndef READPROPS_H
#define READPROPS_H

#include "SignalProcessingTypes.h"

// readProps - read properties from input file and return a map
map<string, string> readProps(const string& inputFile);

#endif // READPROPS_H
