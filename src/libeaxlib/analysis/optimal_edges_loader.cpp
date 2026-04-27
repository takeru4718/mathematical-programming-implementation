#include "optimal_edges_loader.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <cstddef>

namespace eax {
    namespace analysis {
        static std::vector<std::pair<size_t, size_t>> optimal_edges_cache;
        const std::vector<std::pair<size_t, size_t>>& load_optimal_edges(const std::string& file_name) {
            if (!optimal_edges_cache.empty()) {
                return optimal_edges_cache;
            }
    
            std::ifstream file(file_name);
            if (!file.is_open()) {
                throw std::runtime_error("Could not open the file: " + file_name);
            }
    
            std::vector<std::pair<size_t, size_t>> optimal_edges;
            std::vector<size_t> path;
            
            std::string line;
    
            //1行にpathがある
            std::getline(file, line);
            std::stringstream ss(line);
            std::string label;
            //Best Path: というラベルがあるので，読み飛ばす
            std::getline(ss, label, ':');
            size_t path_element;
            //実際のpathを読み込む．ここでは0から始まる番号であることに注意
            while(ss >> path_element){
                path.push_back(path_element);
            }
    
            if (path.size() < 2) {
                throw std::runtime_error("Path is too short: " + file_name);
            }
    
            //辺と都市数は同じ数だけある
            optimal_edges.reserve(path.size());
    
            //pathをoptimal_edgesに変換
            for (size_t i = 0; i < path.size() - 1; ++i) {
                //optimal_edges.first < optimal_edges.secondとなるようにする
                size_t u = std::min(path[i], path[i + 1]);
                size_t v = std::max(path[i], path[i + 1]);
                optimal_edges.emplace_back(u, v);
            }
            //最初と最後の辺を追加
            size_t u = path.front();
            size_t v = path.back();
            if (u > v) std::swap(u, v);
            optimal_edges.emplace_back(u, v);
    
            file.close();
    
            optimal_edges_cache = std::move(optimal_edges);
            return optimal_edges_cache;
        }
    }
}

