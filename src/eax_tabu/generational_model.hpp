#pragma once

#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <optional>

#include "utils.hpp"
#include "genetic_algorithm.hpp"

namespace eax {
    /**
     * @brief GAを実行するオブジェクト
     * @tparam GenerationalStep 1世代の世代交代を行う関数オブジェクト
     * @tparam UpdateFunc 更新処理を行い、世代交代を継続する場合、NotTerminatedを返す関数オブジェクト
     * @tparam LoggingFunc ロギングを行う関数オブジェクト (デフォルトは何もしない関数)
     * @tparam PostProcess 後処理を行う関数オブジェクト (デフォルトは何もしない関数)
     */
    template <typename GenerationalStep, typename UpdateFunc, typename LoggingFunc = mpi::NOP_Function, typename PostProcess = mpi::NOP_Function> 
    class GenerationalModel {
    public:
        /**
         * @brief デフォルトコンストラクタ
         * @note GenerationalStepとPostProcessはデフォルト構築可能である必要があります。
         */
        GenerationalModel()
            requires(std::default_initializable<GenerationalStep> && std::default_initializable<PostProcess>)
            = default;
        
        /**
         * @brief コンストラクタ
         * @param generational_step 1世代の世代交代を行う関数オブジェクト
         * @param post_process 後処理を行う関数オブジェクト (デフォルトは何もしない関数)
         */
        GenerationalModel(GenerationalStep generational_step, UpdateFunc update_func, LoggingFunc logging, PostProcess post_process = {}) 
            : generational_step(std::move(generational_step)), update_func(std::move(update_func)), logging(std::move(logging)), post_process(std::move(post_process)) {}

        template <typename Individual, typename Context>
        constexpr std::pair<mpi::genetic_algorithm::TerminationReason, std::vector<Individual>>
            execute(std::vector<Individual> population, Context& context, size_t generation = 0)
        {
            mpi::genetic_algorithm::TerminationReason reason = update_func(population, context, generation);
            while(reason == mpi::genetic_algorithm::TerminationReason::NotTerminated) {
                logging(population, context, generation);
                generational_step(population, context);
                generation++;
                reason = update_func(population, context, generation);
            }
            
            post_process(population, context, generation, reason);
            return {reason, std::move(population)};
        }
    private:
        GenerationalStep generational_step;
        UpdateFunc update_func;
        LoggingFunc logging;
        PostProcess post_process;
    };

    /**
     * @brief 1世代の世代交代を行う関数オブジェクト
     * @tparam FitnessFunc 適応度を計算する関数オブジェクト
     * @tparam CrossOverFunc 交叉を行う関数オブジェクト
     * @tparam RandomGen 乱数生成器の型
     */
    template <typename FitnessFunc, typename CrossOverFunc>
    class GenerationalStep
    {
    public:
        GenerationalStep(FitnessFunc fitness_func, CrossOverFunc cross_over_func)
            : fitness_func(std::move(fitness_func)), cross_over(std::move(cross_over_func)){}

        /**
         * @brief 世代交代モデルの１回の世代交代を実行する
         * @param population 集団
         * @param context 実行コンテキスト
         * @note contextはrandom_genメンバ変数を持ち、それはコンセプトstd::uniform_random_bit_generatorを満たす型である必要があります。
         */
        template <typename Individual, typename Context>
        void operator()(std::vector<Individual>& population, Context& context)
        {
            using Child = std::invoke_result_t<CrossOverFunc, Individual&, Individual&, Context&>::value_type;
            auto calc_all_fitness = [](const std::vector<Child>& children, Context& context, FitnessFunc& fitness_func) {
                std::vector<double> fitness_values(children.size());
                for (size_t i = 0; i < children.size(); ++i) {
                    fitness_values[i] = fitness_func(children[i], context);
                }
                return fitness_values;
            };
            
            size_t population_size = population.size();

            std::vector<size_t> indices(population_size);
            std::iota(indices.begin(), indices.end(), 0);
            std::shuffle(indices.begin(), indices.end(), context.random_gen);

            population.push_back(population[indices[0]]); // 最初の個体を追加しておく
            indices.push_back(population.size() - 1);
            for (size_t i = 0; i < population_size; ++i) {
                size_t parent_A_index = indices[i];
                size_t parent_B_index = indices[(i + 1) % population_size];
                Individual& parent_A = population[parent_A_index];
                Individual& parent_B = population[parent_B_index];
                std::vector<Child> children = cross_over(parent_A, parent_B, context);

                if (children.empty()) {
                    continue; // 子供が生成されなかった場合はスキップ
                }

                std::vector<double> children_fitness = calc_all_fitness(children, context, fitness_func);
                
                // 子供の中で最良の個体を選択
                size_t best_index = 0;
                double best_fitness = children_fitness[0];
                for (size_t j = 1; j < children_fitness.size(); ++j) {
                    if (children_fitness[j] > best_fitness) {
                        best_fitness = children_fitness[j];
                        best_index = j;
                    }
                }
                
                if (best_fitness > 0.0) {
                    // 子供の適応度が0より大きい場合、親Aを子供に置き換える
                    parent_A = std::move(children[best_index]);
                } else {
                    // 子供の適応度が0以下の場合、親Aをそのままにする
                }
            }

            // 最後の個体を削除
            population.pop_back();
        }
    private:
        FitnessFunc fitness_func;
        CrossOverFunc cross_over;
    };

    //ER世代交代モデル
    /**
     * @brief 1世代の世代交代を行う関数オブジェクト
     * @tparam FitnessFunc 適応度を計算する関数オブジェクト
     * @tparam CrossOverFunc 交叉を行う関数オブジェクト
     * @tparam RandomGen 乱数生成器の型
     */
     template <typename FitnessFunc, typename CrossOverFunc>
     class GenerationalStepER
     {
     public:
         GenerationalStepER(FitnessFunc fitness_func, CrossOverFunc cross_over_func)
             : fitness_func(std::move(fitness_func)), cross_over(std::move(cross_over_func)){}
 
         /**
          * @brief 世代交代モデルの１回の世代交代を実行する
          * @param population 集団
          * @param context 実行コンテキスト
          * @note contextはrandom_genメンバ変数を持ち、それはコンセプトstd::uniform_random_bit_generatorを満たす型である必要があります。
          */
         template <typename Individual, typename Context>
         void operator()(std::vector<Individual>& population, Context& context)
         {
             using Child = std::invoke_result_t<CrossOverFunc, Individual&, Individual&, Context&>::value_type;
             auto calc_all_fitness = [](const std::vector<Child>& children, Context& context, FitnessFunc& fitness_func) {
                 std::vector<double> fitness_values(children.size());
                 for (size_t i = 0; i < children.size(); ++i) {
                     fitness_values[i] = fitness_func(children[i], context);
                 }
                 return fitness_values;
             };
             
             size_t population_size = population.size();
             
             std::vector<size_t> indices(population_size);
             std::iota(indices.begin(), indices.end(), 0);
             std::shuffle(indices.begin(), indices.end(), context.random_gen);

             std::vector<Individual> next_population;
             next_population.reserve(population_size);
 
             //population.push_back(population[indices[0]]); // 最初の個体を追加しておく
             //indices.push_back(population.size()-1);
             //TODO: 母集団のサイズが奇数の場合の処理の追加(エラー吐くなど)
             for (size_t i = 0; i < population_size/2; ++i) {
                 size_t parent_A_index = indices[i*2];
                 size_t parent_B_index = indices[i*2 + 1];
                 Individual& parent_A = population[parent_A_index];
                 Individual& parent_B = population[parent_B_index];
                 //新しい個体の作成
                 Individual new_individual_A = parent_A;
                 Individual new_individual_B = parent_B;
                 //親Aベースの交叉
                std::vector<Child> children_A = cross_over(parent_A, parent_B, context);
                //親Bベースの交叉
                std::vector<Child> children_B = cross_over(parent_B, parent_A, context);

                //子供の集団を作成
                std::vector<Child> children;
                children.reserve(children_A.size() + children_B.size());
                for (auto& child : children_A) {
                    children.emplace_back(std::move(child));
                }
                for (auto& child : children_B) {
                    children.emplace_back(std::move(child));
                }

                if (children.empty()) {
                    continue; // 子供が生成されなかった場合はスキップ
                }
 
                std::vector<double> children_fitness = calc_all_fitness(children, context, fitness_func);
                
                // 1つの子集団から2個体を選択（エリート選択 + ランキング選択）
                if (!children.empty()) {
                    // 1. エリート選択：最良の個体を選択
                    size_t best_index = 0;
                    double best_fitness = children_fitness[0];
                    for (size_t j = 1; j < children_fitness.size(); ++j) {
                        if (children_fitness[j] > best_fitness) {
                            best_fitness = children_fitness[j];
                            best_index = j;
                        }
                    }
                    
                    // エリート個体を取得（適応度が0より大きい場合のみ）
                    std::optional<Child> elite_child;
                    if (best_fitness > 0.0) {
                        elite_child = std::move(children[best_index]);
                    }
                    
                    // 2. ランキング選択：最良以外から1個体を選択
                    // 最良の個体を除外した残りから選択
                    if (children.size() > 1) {
                        // 最良の個体を除外
                        children.erase(children.begin() + best_index);
                        children_fitness.erase(children_fitness.begin() + best_index);
                        
                        if (!children.empty()) {
                            // インデックス配列を初期化
                            std::vector<size_t> child_indices(children.size());
                            std::iota(child_indices.begin(), child_indices.end(), 0);
                            
                            // 適応度に従って降順にソート（最上位がchild_indices[0]、最下位がchild_indices[children.size()-1]）
                            std::sort(child_indices.begin(), child_indices.end(), 
                                [&children_fitness](size_t a, size_t b) {
                                    return children_fitness[a] > children_fitness[b];
                                });
                            
                            // ランクに比例した重みをつける(最上位:children.size(), 最下位:1)
                            std::vector<size_t> rank_weights(children.size());
                            size_t total_rank = 0;
                            for(size_t k = 0; k < children.size(); k++){
                                rank_weights[k] = children.size() - k;  // 最上位に高い重み
                                total_rank += rank_weights[k];
                            }
                            
                            // ランキング選択で1個体を選択
                            size_t pick = context.random_gen() % total_rank;
                            size_t current = 0;
                            size_t selected_index = 0;
                            for(size_t k = 0; k < children.size(); k++){
                                current += rank_weights[k];
                                if(current > pick){
                                    selected_index = child_indices[k];
                                    break;
                                }
                            }
                            
                            // ランキング選択で選んだ個体（適応度が0より大きい場合のみ）
                            double selected_fitness = children_fitness[selected_index];
                            if (selected_fitness > 0.0) {
                                if (elite_child.has_value()) {
                                    // エリート個体でparent_A、ランキング選択個体でparent_Bを置き換え
                                    // 注意: parent_Bに適用する前に、parent_Aの状態を保存する必要がある
                                    // DeltaWithIndividual::apply_toが自動的にベース個体をコピーするので、先にparent_Bを適用
                                    new_individual_B = std::move(children[selected_index]);
                                    new_individual_A = std::move(*elite_child);
                                } else {
                                    // エリート個体がない場合は、ランキング選択個体でparent_Aを置き換え
                                    new_individual_A = std::move(children[selected_index]);
                                }
                            } else if (elite_child.has_value()) {
                                // ランキング選択個体が無効な場合、エリート個体でparent_Aを置き換え
                                new_individual_A = std::move(*elite_child);
                            }
                        } else if (elite_child.has_value()) {
                            // ランキング選択できる個体がない場合、エリート個体でparent_Aを置き換え
                            new_individual_A = std::move(*elite_child);
                        }
                    } else if (elite_child.has_value()) {
                        // 子が1個体しかない場合、エリート個体でparent_Aを置き換え
                        new_individual_A = std::move(*elite_child);
                    }
                }
                next_population.push_back(std::move(new_individual_A));
                next_population.push_back(std::move(new_individual_B));
             }
 
             population = std::move(next_population);
         }
     private:
         FitnessFunc fitness_func;
         CrossOverFunc cross_over;
     };
}