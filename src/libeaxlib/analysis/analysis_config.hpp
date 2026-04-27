#pragma once
#include <string>
#include <cstddef>

namespace eax::analysis {
struct AnalysisConfig {
    //ABサイクルの大きさに対する最適辺の減少数と増加数を計測するかどうか
    bool enable_abcycle_optimal_edges_metrics = false;
    std::string abcycle_metrics_csv_path = "../../analysis/abcycle_size_vs_optimal_edges_decreased.csv";
    std::string optimal_edges_path = "../../optimal_edges/rat575.tour";

    //ABサイクルの大きさを制御するかどうか
    bool enable_abcycle_size_control = false;
    size_t abcycle_size_max = 4;
};
}