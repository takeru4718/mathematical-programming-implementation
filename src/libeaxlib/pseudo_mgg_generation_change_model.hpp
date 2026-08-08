#pragma once

#include <vector>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <random>
#include <stdexcept>

namespace eax {

/**
 * @brief 擬似MGG世代交代モデル
 * @details populationループの偶数番目ではエリート選択、奇数番目ではルーレット選択で
 *          交配相手(parent B)を選び、親子家族から最良個体で parent A を置換する。
 *          生存選択の適応度計算は FitnessFunc に委譲する（実験では greedy を想定）。
 */
template <typename FitnessFunc, typename CrossOverFunc>
class PseudoMggGenerationChangeModel
{
public:
    PseudoMggGenerationChangeModel(FitnessFunc fitness_func, CrossOverFunc cross_over_func)
        : fitness_func(std::move(fitness_func)), cross_over(std::move(cross_over_func)) {}

    template <typename Individual, typename Context>
        requires(requires(std::vector<Individual> population, FitnessFunc fitness_func, CrossOverFunc cross_over, Context context) {
            { fitness_func(cross_over(population[0], population[1], context)[0], context) } -> std::convertible_to<double>;
            population[0] = cross_over(population[0], population[1], context)[0];
            { population[0].get_distance() } -> std::convertible_to<double>;
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

        // ルーレット用: 経路長が短いほど選ばれやすいよう 1/L を使う
        std::vector<double> roulette_weights(population_size);
        for (size_t i = 0; i < population_size; ++i) {
            double length = static_cast<double>(population[i].get_distance());
            roulette_weights[i] = 1.0 / std::max(length, 1.0);
        }

        for (size_t i = 0; i < population_size; ++i) {
            size_t parent_A_index = i;
            size_t parent_B_index = 0;

            if (i % 2 == 0) {
                // 偶数番目: エリート選択（A以外の最良個体）
                parent_B_index = select_elite_mate(population, parent_A_index);
            } else {
                // 奇数番目: ルーレット選択（A以外）
                parent_B_index = select_roulette_mate(roulette_weights, parent_A_index, context.random_gen);
            }

            Individual& parent_A = population[parent_A_index];
            Individual& parent_B = population[parent_B_index];
            std::vector<Child> children = cross_over(parent_A, parent_B, context);

            if (children.empty()) {
                continue;
            }

            children.emplace_back(parent_A);

            std::vector<double> children_fitness = calc_all_fitness(children, context, fitness_func);

            size_t best_index = 0;
            double best_fitness = children_fitness[0];
            for (size_t j = 1; j < children_fitness.size(); ++j) {
                if (children_fitness[j] > best_fitness) {
                    best_fitness = children_fitness[j];
                    best_index = j;
                }
            }

            parent_A = std::move(children[best_index]);

            // A置換後にルーレット重みを更新（以降の奇数番選択で最新の集団を反映）
            double length = static_cast<double>(population[parent_A_index].get_distance());
            roulette_weights[parent_A_index] = 1.0 / std::max(length, 1.0);
        }
    }

private:
    template <typename Individual>
    static size_t select_elite_mate(const std::vector<Individual>& population, size_t excluded_index)
    {
        size_t best_index = excluded_index == 0 ? 1 : 0;
        auto best_length = population[best_index].get_distance();
        for (size_t j = 0; j < population.size(); ++j) {
            if (j == excluded_index) {
                continue;
            }
            auto length = population[j].get_distance();
            if (length < best_length) {
                best_length = length;
                best_index = j;
            }
        }
        return best_index;
    }

    template <typename URBG>
    static size_t select_roulette_mate(const std::vector<double>& weights, size_t excluded_index, URBG& rng)
    {
        std::vector<double> masked = weights;
        masked[excluded_index] = 0.0;
        double sum = std::accumulate(masked.begin(), masked.end(), 0.0);
        if (!(sum > 0.0)) {
            // フォールバック: 一様にA以外を選ぶ
            std::uniform_int_distribution<size_t> uni(0, masked.size() - 2);
            size_t pick = uni(rng);
            return pick >= excluded_index ? pick + 1 : pick;
        }
        std::discrete_distribution<size_t> dist(masked.begin(), masked.end());
        return dist(rng);
    }

    FitnessFunc fitness_func;
    CrossOverFunc cross_over;
};

}
