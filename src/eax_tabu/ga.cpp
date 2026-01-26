#include "ga.hpp"

#include <fstream>

#include "genetic_algorithm.hpp"

#include "object_pools.hpp"
#include "eax_tabu.hpp"
#include "greedy_evaluator.hpp"
#include "entropy_evaluator.hpp"
#include "distance_preserving_evaluator.hpp"

#include "generational_model.hpp"
#include "delta_with_individual.hpp"

namespace eax {
std::pair<mpi::genetic_algorithm::TerminationReason, std::vector<Individual>> execute_ga(
    std::vector<Individual>& population,
    Context& context,
    std::chrono::system_clock::time_point timeout_time,
    const std::string& log_file_name) {

    using namespace std;
    using Context = eax::Context;
    // オブジェクトプール
    eax::ObjectPools object_pools(context.env.tsp.city_count);
    
    // 交叉関数
    eax::EAX_tabu_Rand eax_tabu_rand(object_pools);
    eax::EAX_tabu_N_AB eax_tabu_n_ab(object_pools);
    eax::EAX_tabu_UNIFORM eax_tabu_uniform(object_pools);
    auto crossover_func = [&eax_tabu_rand, &eax_tabu_n_ab, &eax_tabu_uniform](Individual& parent1, const Individual& parent2,
                                Context& context) {
        auto& env = context.env;

        //ERモデルのとき親1と親2のタブーリストを統合
        std::vector<std::pair<size_t, size_t>> merged_tabu_edges;
        const std::vector<std::pair<size_t, size_t>>* tabu_edges_ptr;

        if (env.generational_model_type == GenerationalModelType::ER) {
            // 親1と親2のタブーリストを統合
            merged_tabu_edges.reserve(parent1.get_tabu_edges().size() + parent2.get_tabu_edges().size());
            merged_tabu_edges.insert(merged_tabu_edges.end(), 
                                     parent1.get_tabu_edges().begin(), 
                                     parent1.get_tabu_edges().end());
            merged_tabu_edges.insert(merged_tabu_edges.end(), 
                                     parent2.get_tabu_edges().begin(), 
                                     parent2.get_tabu_edges().end());
            tabu_edges_ptr = &merged_tabu_edges;
        } else {
            // MGGモデルのときは親1のタブーリストのみ
            tabu_edges_ptr = &parent1.get_tabu_edges();
        }

        struct {
            eax::EAX_tabu_Rand& eax_tabu_rand;
            eax::EAX_tabu_N_AB& eax_tabu_n_ab;
            eax::EAX_tabu_UNIFORM& eax_tabu_uniform;
            const Individual& parent1;
            const Individual& parent2;
            Context& context;
            const std::vector<std::pair<size_t, size_t>>* tabu_edges_ptr;

            auto operator()(const EAX_Rand_tag&) {
                return eax_tabu_rand(parent1, parent2, context.env.num_children, context.env.tsp, context.random_gen, {}, {}, std::forward_as_tuple(*tabu_edges_ptr));
            }
            
            auto operator()(const EAX_n_AB_tag& n_ab) {
                return eax_tabu_n_ab(parent1, parent2, context.env.num_children, context.env.tsp, context.random_gen, std::forward_as_tuple(n_ab.get_n()), {}, std::forward_as_tuple(*tabu_edges_ptr));
            }
            
            auto operator()(const EAX_UNIFORM_tag& uniform) {
                return eax_tabu_uniform(parent1, parent2, context.env.num_children, context.env.tsp, context.random_gen, std::forward_as_tuple(uniform.get_ratio()), {}, std::forward_as_tuple(*tabu_edges_ptr));
            }
        } visitor {eax_tabu_rand, eax_tabu_n_ab, eax_tabu_uniform, parent1, parent2, context, tabu_edges_ptr};
        
        auto deltas = std::visit(visitor, env.eax_type);

        std::vector<eax::DeltaWithIndividual<Individual>> children;
        children.reserve(deltas.size());
        for(auto& delta:deltas){
            children.emplace_back(parent1,delta);
        }
        return children;
    };

    // 適応度関数
    auto calc_fitness_lambda = [](const eax::DeltaWithIndividual<Individual>&child,Context& context) {
        return 1.0 / (child.delta.get_delta_distance() + child.individual_ptr->get_distance());
    };
    
    // 更新処理関数
    struct {
        std::chrono::system_clock::time_point timeout_time;
        mpi::genetic_algorithm::TerminationReason operator()(vector<Individual>& population, Context& context, size_t generation) {
            context.current_generation = generation;

            update_individual_and_edge_counts(population, context);
            
            if (std::chrono::system_clock::now() >= timeout_time) {
                return mpi::genetic_algorithm::TerminationReason::TimeLimit;
            }

            return continue_condition(population, context, generation);
        }
        
        void update_individual_and_edge_counts(vector<Individual>& population, Context& context) {
            for (auto& individual : population) {
                auto delta = individual.update(context);
                auto delta_H = eax::calc_delta_entropy(delta, context.pop_edge_counts, context.env.population_size);
                context.entropy += delta_H;
                for (const auto& mod : delta.get_modifications()) {
                    size_t v1 = mod.edge1.first;
                    size_t v2 = mod.edge1.second;
                    size_t new_v2 = mod.new_v2;
                    context.pop_edge_counts[v1][v2] -= 1;
                    context.pop_edge_counts[v1][new_v2] += 1;
                }
            }
        }

        mpi::genetic_algorithm::TerminationReason continue_condition(const vector<Individual>& population, Context& context, size_t generation) {
            double best_length = std::numeric_limits<double>::max();
            double average_length = 0.0;
            for (size_t i = 0; i < population.size(); ++i) {
                double length = population[i].get_distance();
                best_length = std::min(best_length, length);
                average_length += length;
            }
            average_length /= population.size();
            
            if (best_length < context.best_length) {
                context.best_length = best_length;
                context.generation_of_reached_best = generation;
                context.stagnation_generations = 0;
            }else {
                context.stagnation_generations += 1;
            }
            
            if (average_length - best_length < 0.001)
                return mpi::genetic_algorithm::TerminationReason::Converged; // 収束条件
            
            if (generation >= 3000)
                return mpi::genetic_algorithm::TerminationReason::MaxGenerations; // 最大世代数条件
            
            return mpi::genetic_algorithm::TerminationReason::NotTerminated;
        }
    } update_func {timeout_time};
    
    // ロガー
    std::ofstream log_ofs;
    if (!log_file_name.empty()) {
        log_ofs.open(log_file_name, std::ios::out);
        log_ofs << "Generation,BestLength,AverageLength,WorstLength,Entropy" << std::endl;
    }
    struct {
        std::ofstream& log_os;
        void operator()([[maybe_unused]]const vector<Individual>& population, Context& context, size_t generation) {
            if (context.start_time.time_since_epoch().count() == 0) {
                // 計測開始時刻が未設定なら、現在時刻を設定
                const_cast<Context&>(context).start_time = std::chrono::system_clock::now();
            } else {
                auto now = std::chrono::system_clock::now();
                context.elapsed_time += std::chrono::duration<double>(now - context.start_time).count();
                const_cast<Context&>(context).start_time = now;
            }

            if (!log_os.is_open()) {
                return;
            }
            std::vector<double> lengths(population.size());
            for (size_t i = 0; i < population.size(); ++i) {
                lengths[i] = population[i].get_distance();
            }
            double best_length = *std::min_element(lengths.begin(), lengths.end());
            double average_length = std::accumulate(lengths.begin(), lengths.end(), 0.0) / lengths.size();
            double worst_length = *std::max_element(lengths.begin(), lengths.end());
            
            log_os << generation << "," << best_length << "," << average_length << "," << worst_length << "," << context.entropy << std::endl;
        }
    } logging{log_ofs};
    
    struct {
        void operator()([[maybe_unused]]const vector<Individual>& population, Context& context, size_t generation, [[maybe_unused]]mpi::genetic_algorithm::TerminationReason reason) {
            context.final_generation = generation;
        }
    } post_process;

    // 世代交代処理
    if (context.env.generational_model_type == eax::GenerationalModelType::ER) {
        eax::GenerationalStepER generational_step_er(calc_fitness_lambda, crossover_func);
        eax::GenerationalModel genetic_algorithm(generational_step_er, update_func, logging, post_process);
        return genetic_algorithm.execute(population, context, context.current_generation);
    } else {
        // eax::GenerationalStep generational_step(calc_fitness_lambda, crossover_func);
        // eax::GenerationalModel genetic_algorithm(generational_step, update_func, logging, post_process);
        // return genetic_algorithm.execute(population, context, context.current_generation);
    }
}

Context create_context(const std::vector<Individual>& initial_population, Environment const& env) {
    Context context;
    context.env = env;

    context.set_initial_edge_counts(initial_population);
    context.random_gen = std::mt19937(env.random_seed);
    return context;
}

void serialize_population(const std::vector<Individual>& population, std::ostream& os) {
    os << "# Population" << std::endl;
    for (const auto& individual : population) {
        individual.serialize(os);
        os << std::endl;
    }
}

std::vector<Individual> deserialize_population(std::istream& is) {
    std::vector<Individual> population;
    std::string line;
    // # Population
    std::getline(is, line);
    if (line != "# Population") throw std::runtime_error("Expected '# Population'");
    while (std::getline(is, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        Individual individual = Individual::deserialize(iss);
        population.push_back(individual);
    }
    return population;
}
    
}