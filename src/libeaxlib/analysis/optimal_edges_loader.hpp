#pragma once
#include <vector>
#include <string>
#include <utility>
#include <set>

namespace eax {
    namespace analysis {
        const std::set<std::pair<size_t, size_t>>& load_optimal_edges_from_directory(const std::string& dir);
        std::vector<std::string> list_files_in_directory(const std::string& dir);
        std::set<std::pair<size_t, size_t>> load_optimal_edges_from_file(const std::string& file_name);
    }
}
