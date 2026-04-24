#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <random>
#include <fstream>
#include <filesystem>
#include "analysis/optimal_edges_loader.hpp"
#include "analysis/analysis_logger.hpp"

namespace eax {

/**
 * @brief 1世代の世代交代を行う関数オブジェクト
 * @details 「局所的な交叉EAXを用いたGAの高速化とTSPへの適用」(DOI:10.1527/tjsai.22.542)
 *          に記述されている世代交代モデルの実装。
 * @tparam FitnessFunc 適応度を計算する関数オブジェクト
 * @tparam CrossOverFunc 交叉を行う関数オブジェクト
 */
template <typename FitnessFunc, typename CrossOverFunc>
class NagataGenerationChangeModel
{
public:
    NagataGenerationChangeModel(FitnessFunc fitness_func, CrossOverFunc cross_over_func)
        : fitness_func(std::move(fitness_func)), cross_over(std::move(cross_over_func)){}

    /**
     * @brief 世代交代モデルの１回の世代交代を実行する
     * @param population 集団
     * @param context 実行コンテキスト
     * @note contextはrandom_genメンバ変数を持ち、それはコンセプトstd::uniform_random_bit_generatorを満たす型である必要があります。
     */
    template <typename Individual, typename Context>
        requires(requires(std::vector<Individual> population, FitnessFunc fitness_func, CrossOverFunc cross_over, Context context) {
            { fitness_func(cross_over(population[0], population[1], context)[0], context) } -> std::convertible_to<double>;
            population[0] = cross_over(population[0], population[1], context)[0];
            context.random_gen;
        } && std::uniform_random_bit_generator<decltype(Context::random_gen)>)
    void operator()(std::vector<Individual>& population, Context& context)
    {
        //分析用のロガー
        analysis::AnalysisLogger local_logger(context.env.analysis_config);
        using Child = std::invoke_result_t<CrossOverFunc, Individual&, Individual&, Context&>::value_type;
        auto calc_all_fitness = [](const std::vector<Child>& children, Context& context, FitnessFunc& fitness_func) {
            std::vector<double> fitness_values(children.size());
            for (size_t i = 0; i < children.size(); ++i) {
                fitness_values[i] = fitness_func(children[i], context);
            }
            return fitness_values;
        };
        
        size_t population_size = population.size();

        if (population_size < 2) {
            return; // 集団サイズが2未満の場合は何もしない
        }

        std::vector<size_t> indices(population_size);
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), context.random_gen);

        for (size_t i = 0; i < population_size; ++i) {
            size_t parent_A_index = indices[i];
            size_t parent_B_index = indices[(i + 1) % population_size];
            Individual& parent_A = population[parent_A_index];
            Individual& parent_B = population[parent_B_index];
            std::vector<Child> children = cross_over(parent_A, parent_B, context);

            if (children.empty()) {
                continue; // 子供が生成されなかった場合はスキップ
            }

            children.emplace_back(parent_A);

            std::vector<double> children_fitness = calc_all_fitness(children, context, fitness_func);
            
            // 子供and親Aの中で最良の個体を選択
            size_t best_index = 0;
            double best_fitness = children_fitness[0];
            for (size_t j = 1; j < children_fitness.size(); ++j) {
                if (children_fitness[j] > best_fitness) {
                    best_fitness = children_fitness[j];
                    best_index = j;
                }
            }

            if (best_index != children.size() - 1 && local_logger.enabled_abcycle_optimal_edges()){
                //ここでabcycleのサイズを取得する
                const size_t ab_n = children[best_index].get_num_ab_cycle_modifications();
                //最適辺減少数を取得
                const auto& mods = children[best_index].get_modifications();
                const size_t prefix_n = std::min(ab_n, mods.size());//念のため
                //const size_t prefix_n = mods.size();
                //ここで最適辺読み込み(別ファイルから読み込むものとする)
                const auto& optimal_edges = analysis::load_optimal_edges(local_logger.optimal_edges_path());
                //最適辺が減少した数を保存
                size_t num_optimal_edges_decreased = 0;
                //最適辺が増加した数を保存
                size_t num_optimal_edges_increased = 0;
                //Modificationを逆に辿って消える辺を取得する
                for (size_t i = 0; i < prefix_n; ++i) {
                    const auto& m = mods[i];
                    auto [v1, v2] = m.edge1;
                    auto new_v2 = m.new_v2;
                    //modificationは順方向と逆方向の両方の情報を持っているため，v1 < v2のみを考える
                    if (v1 < v2 && std::find(optimal_edges.begin(), optimal_edges.end(), std::make_pair(v1, v2)) != optimal_edges.end()) {
                        num_optimal_edges_decreased++;
                    }
                    if (v1 < new_v2 && std::find(optimal_edges.begin(), optimal_edges.end(), std::make_pair(v1, new_v2)) != optimal_edges.end()) {
                        num_optimal_edges_increased++;
                    }
                }

                //ここでabcycleのサイズに対する最適辺の減少数をファイル書き出し
                local_logger.append_abcycle_optimal_edges(ab_n, num_optimal_edges_decreased, num_optimal_edges_increased);
            }


            parent_A = std::move(children[best_index]);
        }

    }
private:
    FitnessFunc fitness_func;
    CrossOverFunc cross_over;
};

}