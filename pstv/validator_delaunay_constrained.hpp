#ifndef VALIDATOR_CONSTRAINED_DELAUNAY_HPP
#define VALIDATOR_CONSTRAINED_DELAUNAY_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <pstv/dataset.hpp>
#include <pstv/validation.hpp>
#include <vector>
#include <string>
#include <unordered_set>

using std::vector;

namespace pstv {
class ValidatorConstrainedDelaunay {
  pstv::Validation validation;
  pstv::Dataset dataset;
  int count_non_delaunay = 0;
  int count_constrained_skipped = 0;
  std::unordered_set<std::string> constrained_edges;

public:
  explicit ValidatorConstrainedDelaunay(){};

  // Add a single constrained edge by vertex indices.
  // The edge is stored in normalized form "u-v" where u < v.
  void add_constrained_edge(int u, int v) {
      if (u > v) std::swap(u, v);
      constrained_edges.insert(std::to_string(u) + "-" + std::to_string(v));
  }

  // Load constrained edges from a CSV file.
  // Each line should contain two vertex indices separated by a comma: "u,v"
  // Lines starting with '#' are treated as comments and ignored.
  // Empty lines are ignored.
  // Returns the number of constrained edges loaded.
  int load_constrained_edges_csv(const std::string& filename) {
      std::ifstream ifs(filename);
      if (!ifs.is_open()) {
          std::cerr << "Warning: Could not open constrained edges file: " << filename << std::endl;
          return 0;
      }
      int count = 0;
      std::string line;
      while (std::getline(ifs, line)) {
          // Skip empty lines and comments
          if (line.empty() || line[0] == '#') continue;
          // Remove carriage return if present (Windows line endings)
          if (!line.empty() && line.back() == '\r') line.pop_back();
          if (line.empty()) continue;

          std::stringstream ss(line);
          std::string token;
          int u = -1, v = -1;
          if (std::getline(ss, token, ',')) {
              u = std::atoi(token.c_str());
          }
          if (std::getline(ss, token, ',')) {
              v = std::atoi(token.c_str());
          }
          if (u >= 0 && v >= 0 && u != v) {
              add_constrained_edge(u, v);
              count++;
          }
      }
      std::cout << "Loaded " << count << " constrained edges from: " << filename << std::endl;
      return count;
  }

  // Get number of registered constrained edges (including boundary edges after validate)
  size_t get_constrained_edge_count() const {
      return constrained_edges.size();
  }

  // Get number of constrained edges that were skipped during validation
  int get_constrained_skipped_count() const {
      return count_constrained_skipped;
  }

  bool validate(pstv::Dataset ds, bool mode = false) {
    dataset = ds;
    setup_data();
    
    // Automatically treat boundary edges as constrained
    // (boundary_edges is now properly populated by set_boundary_edges() in setup_data())
    for (const auto& edge : dataset.boundary_edges) {
        constrained_edges.insert(edge);
    }

    std::cout << "Constrained edges (including boundary): " << constrained_edges.size() << std::endl;

    while (0 < dataset.edges.size()) {
      auto edge_itr = dataset.edges.begin();
      std::string current_edge = *edge_itr;

      // Check if the edge is constrained -> skip Delaunay check and flipping
      if (constrained_edges.find(current_edge) != constrained_edges.end()) {
          dataset.edges.erase(current_edge);
          count_constrained_skipped++;
          continue; 
      }

      auto shared_triangles = dataset.edge_map[current_edge];
      if (shared_triangles.size() < 2) {
        dataset.edges.erase(current_edge);
      } else {
        auto triangle_itr = shared_triangles.begin();
        int triangle1_index = *triangle_itr;
        triangle_itr++;
        int triangle2_index = *triangle_itr;
        if (mode) {
          _validate_and_flip(triangle1_index, triangle2_index);
        } else {
          _validate(triangle1_index, triangle2_index);
        }
      }
    }

    std::cout << "Constrained edges skipped: " << count_constrained_skipped << std::endl;

    if (0 < count_non_delaunay) {
      std::cout << "Delaunay is NOT verified !" << std::endl;
      std::cout << "Non-Delaunay edges: " << count_non_delaunay << std::endl;
      return false;
    }
    std::cout << "Delaunay is verified !" << std::endl;
    return true;
  }

  void output_new_triangles_data(std::string filename) {
    std::ofstream ofs(filename);
    for (size_t i = 0; i < dataset.triangles.size(); ++i) {
      ofs << dataset.triangles[i][0] << "," << dataset.triangles[i][1] << "," << dataset.triangles[i][2] << std::endl;
    }
    std::cout << "Output New Triangles Data to: " << filename << std::endl;
  }

private:
  void setup_data() {
    dataset.set_edge_map();
    dataset.set_boundary_edges();   // FIX: initialize boundary_edges correctly
    dataset.set_edges();
    dataset.set_vertex_map();
    dataset.set_triangle_map();
  }

  // Helper: check if an edge should be re-inserted into the processing queue
  bool should_reinsert_edge(const std::string& edge) {
    // Don't re-insert if it's a boundary edge or a constrained edge
    if (dataset.boundary_edges.find(edge) != dataset.boundary_edges.end()) return false;
    if (constrained_edges.find(edge) != constrained_edges.end()) return false;
    return true;
  }

  void _validate(const int triangle_index1, const int triangle_index2) {
    vector<int> triangle1 = dataset.triangles[triangle_index1];
    vector<int> triangle2 = dataset.triangles[triangle_index2];
    vector<int> all = triangle1;
    all.insert(all.end(), triangle2.begin(), triangle2.end());
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    vector<int> non_shared1 = all;
    vector<int> non_shared2 = all;
    vector<int> shared = all;
    for (size_t i = 0; i < 3; ++i) {
      non_shared1.erase(std::remove(non_shared1.begin(), non_shared1.end(), triangle2[i]), non_shared1.end());
    }
    for (size_t i = 0; i < 3; ++i) {
      non_shared2.erase(std::remove(non_shared2.begin(), non_shared2.end(), triangle1[i]), non_shared2.end());
    }
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared1[0]), shared.end());
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared2[0]), shared.end());
    if (shared[0] < shared[1]) {
      dataset.edges.erase(std::to_string(shared[0]) + '-' + std::to_string(shared[1]));
    } else {
      dataset.edges.erase(std::to_string(shared[1]) + '-' + std::to_string(shared[0]));
    }
    if (validation.incircle(dataset.vertexes[triangle1[0]], dataset.vertexes[triangle1[1]], dataset.vertexes[triangle1[2]], dataset.vertexes[non_shared2[0]]) < 0) {
      count_non_delaunay++;
    }
  }

  void _validate_and_flip(const int triangle_index1, const int triangle_index2) {
    vector<int> triangle1 = dataset.triangles[triangle_index1];
    vector<int> triangle2 = dataset.triangles[triangle_index2];
    vector<int> all = triangle1;
    all.insert(all.end(), triangle2.begin(), triangle2.end());
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    vector<int> non_shared1 = all;
    vector<int> non_shared2 = all;
    vector<int> shared = all;
    for (size_t i = 0; i < 3; ++i) {
      non_shared1.erase(std::remove(non_shared1.begin(), non_shared1.end(), triangle2[i]), non_shared1.end());
    }
    for (size_t i = 0; i < 3; ++i) {
      non_shared2.erase(std::remove(non_shared2.begin(), non_shared2.end(), triangle1[i]), non_shared2.end());
    }
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared1[0]), shared.end());
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared2[0]), shared.end());
    if (shared[0] < shared[1]) {
      dataset.edges.erase(std::to_string(shared[0]) + '-' + std::to_string(shared[1]));
    } else {
      dataset.edges.erase(std::to_string(shared[1]) + '-' + std::to_string(shared[0]));
    }
    if (validation.incircle(dataset.vertexes[triangle1[0]], dataset.vertexes[triangle1[1]], dataset.vertexes[triangle1[2]], dataset.vertexes[non_shared2[0]]) < 0) {
      count_non_delaunay++;
      
      // ... Flip Logic ...
      vector<int> non_shared;
      non_shared.push_back(non_shared1[0]);
      non_shared.push_back(non_shared2[0]);
      vector<int> new_triangle1 = non_shared;
      vector<int> new_triangle2 = non_shared;
      new_triangle1.push_back(shared[0]);
      new_triangle2.push_back(shared[1]);
      std::sort(new_triangle1.begin(), new_triangle1.end());
      std::sort(new_triangle2.begin(), new_triangle2.end());
      dataset.triangles[triangle_index1] = new_triangle1;
      dataset.triangles[triangle_index2] = new_triangle2;

      // update edge_map
      std::string key_edge;
      if (shared[0] < shared[1]) {
        key_edge = std::to_string(shared[0]) + '-' + std::to_string(shared[1]);
      } else {
        key_edge = std::to_string(shared[1]) + '-' + std::to_string(shared[0]);
      }
      dataset.edge_map.erase(key_edge);
      std::string new_key_edge;
      if (non_shared1[0] < non_shared2[0]) {
        new_key_edge = std::to_string(non_shared1[0]) + '-' + std::to_string(non_shared2[0]);
      } else {
        new_key_edge = std::to_string(non_shared2[0]) + '-' + std::to_string(non_shared1[0]);
      }
      std::unordered_set<int> triangles_index;
      triangles_index.insert(triangle_index1);
      triangles_index.insert(triangle_index2);
      dataset.edge_map[new_key_edge] = triangles_index;

      // Re-add adjacent edges to processing queue
      // FIX: also check constrained_edges (not just boundary_edges)
      std::string edge;
      triangles_index.clear();
      for (int i = 0; i < 2; ++i) {
        if (non_shared1[0] < shared[i]) {
          edge = std::to_string(non_shared1[0]) + '-' + std::to_string(shared[i]);
          if (should_reinsert_edge(edge))
            dataset.edges.insert(edge);
          triangles_index = dataset.edge_map[edge];
          triangles_index.erase(triangle_index1);
          triangles_index.erase(triangle_index2);
          if (i == 0) {
            triangles_index.insert(triangle_index1);
          } else {
            triangles_index.insert(triangle_index2);
          }
          dataset.edge_map[edge] = triangles_index;
        } else {
          edge = std::to_string(shared[i]) + '-' + std::to_string(non_shared1[0]);
          if (should_reinsert_edge(edge))
            dataset.edges.insert(edge);
          triangles_index = dataset.edge_map[edge];
          triangles_index.erase(triangle_index1);
          triangles_index.erase(triangle_index2);
          if (i == 0) {
            triangles_index.insert(triangle_index1);
          } else {
            triangles_index.insert(triangle_index2);
          }
          dataset.edge_map[edge] = triangles_index;
        }
        if (non_shared2[0] < shared[i]) {
          edge = std::to_string(non_shared2[0]) + '-' + std::to_string(shared[i]);
          if (should_reinsert_edge(edge))
            dataset.edges.insert(edge);
          triangles_index = dataset.edge_map[edge];
          triangles_index.erase(triangle_index1);
          triangles_index.erase(triangle_index2);
          if (i == 0) {
            triangles_index.insert(triangle_index1);
          } else {
            triangles_index.insert(triangle_index2);
          }
          dataset.edge_map[edge] = triangles_index;
        } else {
          edge = std::to_string(shared[i]) + '-' + std::to_string(non_shared2[0]);
          if (should_reinsert_edge(edge))
            dataset.edges.insert(edge);
          triangles_index = dataset.edge_map[edge];
          triangles_index.erase(triangle_index1);
          triangles_index.erase(triangle_index2);
          if (i == 0) {
            triangles_index.insert(triangle_index1);
          } else {
            triangles_index.insert(triangle_index2);
          }
          dataset.edge_map[edge] = triangles_index;
        }
      }
    }
  }
};
} // namespace pstv

#endif // VALIDATOR_CONSTRAINED_DELAUNAY_HPP
