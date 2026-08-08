#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <random>
#include <stdexcept>

namespace eax {

/**
 * @brief 擬似MGGの奇数番生存選択方式
 */
enum class PseudoMggOddSelection {
    Roulette,
    Ranking,
};

/**
 * @brief 擬似MGG世代交代モデル
 * @details 親ペアリングは Nagata モデルと同じ（shuffle + 環状）。
 *          生存選択だけを変え、ループ番号 i が
 *          - 偶数: 家族（子 + 親A）からエリート選択
 *          - 奇数: 家族（子 + 親A）からルーレット or ランキング選択
 *          で parent A を置換する。
 *          適応度計算は FitnessFunc に委譲（実験では greedy を想定）。
 *          集団最良保護は未実装（保留）。
 */
template <typename FitnessFunc, typename CrossOverFunc>
class PseudoMggGenerationChangeModel
{
public:
    PseudoMggGenerationChangeModel(FitnessFunc fitness_func, CrossOverFunc cross_over_func,
                                   PseudoMggOddSelection odd_selection = PseudoMggOddSelection::Roulette)
        : fitness_func(std::move(fitness_func)),
          cross_over(std::move(cross_over_func)),
          odd_selection(odd_selection) {}

    template <typename Individual, typename Context>
        requires(requires(std::vector<Individual> population, FitnessFunc fitness_func, CrossOverFunc cross_over, Context context) {
            { fitness_func(cross_over(population[0], population[1], context)[0], context) } -> std::convertible_to<double>;
            population[0] = cross_over(population[0], population[1], context)[0];
            context.random_gen;
        } && std::uniform_random_bit_generator<decltype(Context::random_gen)>)
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
        if (population_size < 2) {
            return;
        }

        // 親ペアリングは従来Nagataと同じ: shuffle後の環状ペア
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
                continue;
            }

            // 家族 = 子 + 親A
            children.emplace_back(parent_A);
            std::vector<double> family_fitness = calc_all_fitness(children, context, fitness_func);

            size_t selected_index = 0;
            if (i % 2 == 0) {
                // 偶数番ループ: エリート選択
                selected_index = select_elite(family_fitness);
            } else if (odd_selection == PseudoMggOddSelection::Ranking) {
                // 奇数番ループ: 線形ランキング選択（最悪1 : 最良3）
                selected_index = select_ranking(family_fitness, context.random_gen);
            } else {
                // 奇数番ループ: ルーレット選択
                selected_index = select_roulette(family_fitness, context.random_gen);
            }

            parent_A = std::move(children[selected_index]);
        }
    }

private:
    static size_t select_elite(const std::vector<double>& fitness_values)
    {
        size_t best_index = 0;
        double best_fitness = fitness_values[0];
        for (size_t j = 1; j < fitness_values.size(); ++j) {
            if (fitness_values[j] > best_fitness) {
                best_fitness = fitness_values[j];
                best_index = j;
            }
        }
        return best_index;
    }

    template <typename URBG>
    static size_t select_roulette(const std::vector<double>& fitness_values, URBG& rng)
    {
        // greedy適応度は負になり得るので、minを引いて非負化してからルーレットする
        double min_fitness = *std::min_element(fitness_values.begin(), fitness_values.end());
        std::vector<double> weights(fitness_values.size());
        for (size_t j = 0; j < fitness_values.size(); ++j) {
            weights[j] = fitness_values[j] - min_fitness + 1e-12;
        }
        double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
        if (!(sum > 0.0)) {
            std::uniform_int_distribution<size_t> uni(0, weights.size() - 1);
            return uni(rng);
        }
        std::discrete_distribution<size_t> dist(weights.begin(), weights.end());
        return dist(rng);
    }

    /**
     * @brief 線形ランキング選択
     * @details 適応度降順の順位 r=0..n-1（0が最良）に対し
     *          weight = eta_max - (eta_max - eta_min) * r / (n-1)
     *          を使う。ここでは eta_min=1, eta_max=3（最悪1 : 最良3）。
     */
    template <typename URBG>
    static size_t select_ranking(const std::vector<double>& fitness_values, URBG& rng,
                                 double eta_min = 1.0, double eta_max = 3.0)
    {
        const size_t n = fitness_values.size();
        if (n == 1) {
            return 0;
        }

        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return fitness_values[a] > fitness_values[b];
        });

        std::vector<double> weights(n, 0.0);
        for (size_t rank = 0; rank < n; ++rank) {
            double weight = eta_max - (eta_max - eta_min) * static_cast<double>(rank) / static_cast<double>(n - 1);
            weights[order[rank]] = weight;
        }

        std::discrete_distribution<size_t> dist(weights.begin(), weights.end());
        return dist(rng);
    }

    FitnessFunc fitness_func;
    CrossOverFunc cross_over;
    PseudoMggOddSelection odd_selection;
};

}
