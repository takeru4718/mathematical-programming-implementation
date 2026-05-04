#pragma once

#include "eaxdef.hpp"

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <cassert>

//回転 + 逆順同値を潰した代表系を返す
inline std::vector<size_t> canonicalize_tour(const std::vector<size_t>& path) {
    if (path.empty()) return {};

    const size_t n = path.size();

    //巡回路は最低3頂点以上であることを確認する
    assert(n >=  3);

    // 最小頂点を先頭に移動
    auto it_min = std::min_element(path.begin(), path.end());
    const size_t min_pos = static_cast<size_t>(std::distance(path.begin(), it_min));

    //辞書順で比較する=最小頂点につながる2つの頂点の大小で順方向，逆方向を判定する
    if (path[(min_pos + 1) % n] < path[(min_pos + n - 1) % n]) {
        //順方向
        std::vector<size_t> forward;
        forward.reserve(n);
        for (size_t k = 0; k < n; ++k) {
            forward.push_back(path[(min_pos + k) % n]);
        }
        return forward;
    } else {
        //逆方向
        std::vector<size_t> backward;
        backward.reserve(n);
        for (size_t k = 0; k < n; ++k) {
            backward.push_back(path[(min_pos + n - k) % n]);
        }
        return backward;
    }
}

namespace eax {

template <individual_readable Individual>
void print_best_solution(const std::vector<Individual>& population, std::ostream& os) {
    size_t best_index = 0;
    int64_t best_length = population[0].get_distance();
    for (size_t i = 1; i < population.size(); ++i) {
        int64_t length = population[i].get_distance();
        if (length < best_length) {
            best_length = length;
            best_index = i;
        }
    }
    
    auto& best_ind = population[best_index];

    std::vector<size_t> best_path;
    size_t prev = 0;
    size_t current = 0;
    for (size_t i = 0; i < best_ind.size(); ++i) {
        best_path.push_back(current);
        size_t next = best_ind[current][0];
        if (next == prev) {
            next = best_ind[current][1];
        }
        prev = current;
        current = next;
    }
    
    os << "Best Solution: " << best_length << std::endl;
    os << "Best Path: ";
    for (size_t city : best_path) {
        os << city << " ";
    }
    os << std::endl;
}

template <individual_readable Individual>
void print_best_solution_unique(const std::vector<Individual>& population, std::ostream& os) {
    int64_t best_length = population[0].get_distance();
    for (size_t i = 1; i < population.size(); ++i) {
        int64_t length = population[i].get_distance();
        if (length < best_length) {
            best_length = length;
        }
    }

    //ここで同じ最良解を持つ解を取得する
    std::vector<size_t> best_indices;
    for (size_t i = 0; i < population.size(); ++i) {
        if (population[i].get_distance() == best_length) {
            best_indices.push_back(i);
        }
    }

    //ここで重複ツアーは最初に出現した一回のみを記録するようにする
    std::set<std::vector<size_t>> seen;
    std::vector<size_t> unique_best_indices;

    for (size_t idx : best_indices) {
        auto& best_ind = population[idx];
        std::vector<size_t> best_path;
        size_t prev = 0;
        size_t current = 0;
        for (size_t i = 0; i < best_ind.size(); ++i) {
            best_path.push_back(current);
            size_t next = best_ind[current][0];
            if (next == prev) {
                next = best_ind[current][1];
            }
            prev = current;
            current = next;
        }
        
        auto canon = canonicalize_tour(best_path);
        //std::set::insertの戻り値の使用で，firstはその要素の位置，secondはboolでその要素が追加されたかどうかを示す
        if (seen.insert(canon).second) { // 初出のみtrue
            unique_best_indices.push_back(idx);
        }
    }

    //ここで最良解を保存する
    os << "Best Solution: " << best_length << std::endl;
    for (size_t i = 0; i < unique_best_indices.size(); ++i) {
        auto& best_ind = population[unique_best_indices[i]];
        std::vector<size_t> best_path;
        size_t prev = 0;
        size_t current = 0;
        for (size_t i = 0; i < best_ind.size(); ++i) {
            best_path.push_back(current);
            size_t next = best_ind[current][0];
            if (next == prev) {
                next = best_ind[current][1];
            }
            prev = current;
            current = next;
        }
        os << "Best Path: ";
        for (size_t city : best_path) {
            os << city << " ";
        }
        os << std::endl;
    }
}

}