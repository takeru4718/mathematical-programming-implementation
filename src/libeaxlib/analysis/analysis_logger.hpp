#pragma once

#include <cstddef>
#include <string>

#include "analysis_config.hpp"

namespace eax::analysis {
class AnalysisLogger {
public:
    explicit AnalysisLogger(const AnalysisConfig& config)
        : config_(config){}

    bool enabled_abcycle_optimal_edges() const {
        return config_.enable_abcycle_optimal_edges_metrics;
    }

    const std::string& optimal_edges_path() const {
        return config_.optimal_edges_path;
    }

    void append_abcycle_optimal_edges(std::size_t abcycle_size,
                                    std::size_t num_optimal_edges_decreased,
                                    std::size_t num_optimal_edges_increased) const;
    
private:
    const AnalysisConfig& config_;
                                    
};

} // namespace eax::analysis