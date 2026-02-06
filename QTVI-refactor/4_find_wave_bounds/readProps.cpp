// ============================================================================
// File: readProps.cpp
// ============================================================================
#include "readProps.h"

map<string, string> readProps(const string& inputFile) {
    map<string, string> props;
    std::ifstream fid(inputFile);

    if (!fid.is_open()) {
        throw std::runtime_error("Cannot open file: " + inputFile);
    }

    string line;
    while (std::getline(fid, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Skip empty lines or comments
        if (line.empty() || line[0] == '#') {
            continue;
        }

        // Split by '='
        size_t pos = line.find('=');
        if (pos != string::npos) {
            string key = line.substr(0, pos);
            string value = line.substr(pos + 1);

            // Trim key and value
            key.erase(0, key.find_first_not_of(" \t\r\n"));
            key.erase(key.find_last_not_of(" \t\r\n") + 1);
            value.erase(0, value.find_first_not_of(" \t\r\n"));
            value.erase(value.find_last_not_of(" \t\r\n") + 1);

            props[key] = value;
        }
    }

    fid.close();
    return props;
}