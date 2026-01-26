#pragma once

#include <vector>
#include <array>
#include <iostream>
#include <cmath>
#include <chrono>

#include "checksumed.hpp"
#include "crossover_delta.hpp"
#include "delta_with_individual.hpp"

namespace eax {

    class Context;

    class Individual : public Checksumed {
    public:
        Individual(const std::vector<size_t>& path, const std::vector<std::vector<int64_t>>& adjacency_matrix, size_t tabu_range);
        constexpr std::array<size_t, 2>& operator[](size_t index) {
            return doubly_linked_list[index];
        }

        constexpr std::array<size_t, 2> const& operator[](size_t index) const {
            return doubly_linked_list[index];
        }

        size_t size() const;
        size_t get_distance() const {
            return distance_;
        }
        std::vector<size_t> to_path() const;
        // Individual& operator=(const DeltaWithIndividual& child) {
        //     prev_diff = child;
        //     return *this;
        // }

        template <typename T = Individual>
            requires std::is_same_v<T, Individual>
        Individual& operator=(eax::DeltaWithIndividual<T>&& delta_view) {
            delta_view.apply_to(*this);
            prev_diff = delta_view.delta;
            return *this;
        }
        
        std::vector<std::pair<size_t, size_t>> const& get_tabu_edges() const {
            return tabu_edges[current_tabu_index];
        }
        
        CrossoverDelta update(eax::Context& context);
        
        int64_t& distance() {
            return distance_;
        }
        
        void serialize(std::ostream& os) const;
        static Individual deserialize(std::istream& is);
    private:
        Individual() = default;
        std::vector<std::array<size_t, 2>> doubly_linked_list;
        CrossoverDelta prev_diff;
        int64_t distance_;
        size_t tabu_range;
        size_t current_tabu_index = 0;
        std::vector<std::vector<std::pair<size_t, size_t>>> tabu_edges;
    };
}