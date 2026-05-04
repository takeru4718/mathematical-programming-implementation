#pragma once
#include <string>
#include <cstddef>

namespace eax::analysis {
struct AnalysisConfig {
    //ABサイクルの大きさに対する最適辺の減少数と増加数を計測するかどうか
    //analysis_logger.cppで使用
    bool enable_abcycle_optimal_edges_metrics = false;
    std::string abcycle_metrics_csv_path = "../../analysis/xpr2308_greedy_tabu/abcycle_size_vs_optimal_edges_decreased_with_generation.csv";
    std::string optimal_edges_directory = "../../result/xpr2308_opt_tour/";

    //ABサイクルの分布を取得するフラグ
    bool enable_abcycle_distribution_metrics = false;
    std::string abcycle_distribution_csv_path = "../../analysis/xpr2308_greedy_tabu/abcycle_distribution.csv";

    //ABサイクルの大きさを制御するかどうか
    bool enable_abcycle_size_control = false;
    size_t abcycle_size_max = 4;
};
}