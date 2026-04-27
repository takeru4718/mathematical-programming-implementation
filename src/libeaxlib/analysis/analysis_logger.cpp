#include "analysis_logger.hpp"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace eax::analysis {

void AnalysisLogger::append_abcycle_optimal_edges(std::size_t abcycle_size,
                                    std::size_t num_optimal_edges_decreased,
                                    std::size_t num_optimal_edges_increased ) const {
    if (!config_.enable_abcycle_optimal_edges_metrics) return;
    
    const bool need_header = 
        !std::filesystem::exists(config_.abcycle_metrics_csv_path) ||
        std::filesystem::file_size(config_.abcycle_metrics_csv_path) == 0;

    std::ofstream out(config_.abcycle_metrics_csv_path, std::ios::app);
    if (!out.is_open()) {
        throw std::runtime_error("Failed to open abcycle metrics CSV file: " + config_.abcycle_metrics_csv_path);
    }

    if (need_header) {
        out << "abcycle_size,num_optimal_edges_decreased,num_optimal_edges_increased\n";
    }

    out << abcycle_size << "," << num_optimal_edges_decreased << "," << num_optimal_edges_increased << std::endl;
    out.close();
}

} // namespace eax::analysis