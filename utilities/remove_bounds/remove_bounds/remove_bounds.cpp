#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>

        namespace fs = std::filesystem;

        // Helper to format EDF header strings (8 characters, space padded)
        std::string formatEDF(std::string val) {
            while (val.length() < 8) val += " ";
            return val.substr(0, 8);
        }

        int main() {
            // ========================================================================
            // HARDCODED PATH
            // ========================================================================
            std::string input_path = R"(D:\USERS\MiraWelner\QTVI\QTVI-data-files\0_original_files\mesa_files\3010112_20110725_EDF\3010112_20110725.edf)";
            // ========================================================================

            std::ifstream in(input_path, std::ios::binary);
            if (!in.is_open()) {
                std::cerr << "Error: Could not open file." << std::endl;
                return 1;
            }

            fs::path p(input_path);
            std::string output_path = (p.parent_path() / (p.stem().string() + "_unbounded.edf")).string();
            std::ofstream out(output_path, std::ios::binary);

            // 1. Copy Fixed Header (256 bytes)
            char fixed_header[256];
            in.read(fixed_header, 256);
            int ns = std::stoi(std::string(fixed_header + 252, 4));
            out.write(fixed_header, 256);

            // 2. Read the entire variable header into memory (256 bytes * ns)
            std::vector<char> v_header(256 * ns);
            in.read(v_header.data(), 256 * ns);

            // 3. Modify limits for every signal in the header
            for (int i = 0; i < ns; ++i) {
                int offset = i * 256;

                // Field Offsets within the 256-byte signal block:
                // Label: 0, Transducer: 16, Phys Dim: 96
                // Phys Min: 104, Phys Max: 112, Dig Min: 120, Dig Max: 128

                std::string p_min = formatEDF("-32768");
                std::string p_max = formatEDF("32767");
                std::string d_min = formatEDF("-32768");
                std::string d_max = formatEDF("32767");

                // Overwrite Physical Limits
                memcpy(&v_header[offset + 104], p_min.c_str(), 8);
                memcpy(&v_header[offset + 112], p_max.c_str(), 8);

                // Overwrite Digital Limits (Crucial for MESA data)
                memcpy(&v_header[offset + 120], d_min.c_str(), 8);
                memcpy(&v_header[offset + 128], d_max.c_str(), 8);
            }

            // 4. Write modified variable header
            out.write(v_header.data(), 256 * ns);

            // 5. Copy the actual signal data
            out << in.rdbuf();

            std::cout << "Successfully 'unbounded' MESA EDF." << std::endl;
            std::cout << "Physical and Digital limits synced to [-32768, 32767]." << std::endl;
            return 0;
        }
