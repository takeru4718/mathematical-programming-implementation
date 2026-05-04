#include "optimal_edges_loader.hpp"
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <cstddef>
#include <set>
#include <filesystem>
#include <algorithm>

namespace eax {
    namespace analysis {
        static std::set<std::pair<size_t, size_t>> optimal_edges_cache;
        static std::string cached_dir;

        const std::set<std::pair<size_t, size_t>>& load_optimal_edges_from_directory(const std::string& dir) {
            if (!optimal_edges_cache.empty()) {
                //ディレクトリが異なる場合は想定していないため，エラーを投げる
                if (dir != cached_dir) {
                    throw std::runtime_error("optimal_edges cache is already initialized with: " + cached_dir);
                }
                return optimal_edges_cache;
            }

            //ディレクトリ内のファイルを取得
            std::vector<std::string> files = list_files_in_directory(dir);

            if (files.empty()) {
                throw std::runtime_error("No files found in the directory: " + dir);
            }

            std::set<std::pair<size_t, size_t>> unique_edges;
            //ファイルを読み込む
            for (const auto& file : files) {
                std::set<std::pair<size_t, size_t>> edges = load_optimal_edges_from_file(file);
                unique_edges.insert(edges.begin(), edges.end());
            }

            optimal_edges_cache = std::move(unique_edges);
            cached_dir = dir;
            return optimal_edges_cache;
        }

        std::vector<std::string> list_files_in_directory(const std::string& dir) {
            std::vector<std::string> files;

            for (const auto& entry : std::filesystem::directory_iterator(dir)) {
                if (!entry.is_regular_file()) continue;
                files.push_back(entry.path().string());
            }

            std::sort(files.begin(), files.end()); // 再現性のために並べる
            return files;
        }

        std::set<std::pair<size_t, size_t>> load_optimal_edges_from_file(const std::string& file_name) {
    
            std::ifstream file(file_name);
            if (!file.is_open()) {
                throw std::runtime_error("Could not open the file: " + file_name);
            }
    
            std::set<std::pair<size_t, size_t>> unique_edges;
            
            std::string line;

            //複数のpathがある場合にpathを読み込む
            while(std::getline(file, line)) {
                //Best Path: というラベルがない場合は読み飛ばす→Best Length:の行を読み飛ばす
                if (line.rfind("Best Path:", 0) != 0) continue; 
                std::vector<size_t> path;
                std::stringstream ss(line);
                std::string label;
                size_t path_element;
                //Best Path: というラベルがあるので，読み飛ばす
                std::getline(ss, label, ':');
                //実際のpathを読み込む．ここでは0から始まる番号であることに注意
                while(ss >> path_element){
                    path.push_back(path_element);
                }

                if (path.size() < 2) {
                    throw std::runtime_error("Path is too short: " + file_name);
                }

                for (size_t i = 0; i < path.size() - 1; ++i) {
                    size_t u = std::min(path[i], path[i + 1]);
                    size_t v = std::max(path[i], path[i + 1]);
                    unique_edges.emplace(u, v);
                }
                //最初と最後の辺を追加
                size_t u = path.front();
                size_t v = path.back();
                if (u > v) std::swap(u, v);
                unique_edges.emplace(u, v);
            }
    
            return unique_edges;
        }
    }
}

