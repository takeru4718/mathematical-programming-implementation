#pragma once
#include <vector>
#include <string>
#include <utility>

namespace eax {
    namespace analysis {
        const std::vector<std::pair<size_t, size_t>>& load_optimal_edges(const std::string& file_name);
    }
}
