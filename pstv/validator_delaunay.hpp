#ifndef VALIDATOR_DELAUNAY_HPP
#define VALIDATOR_DELAUNAY_HPP

#include <pstv/dataset.hpp>
#include <pstv/validation.hpp>
#include <vector>

namespace pstv {
class ValidatorDelaunay {
  pstv::Validation validation;
  pstv::Dataset dataset;
  int count_non_delaunay = 0;

public:
  explicit ValidatorDelaunay(){};
  bool validate(pstv::Dataset ds, bool mode = false) {
    dataset = ds;
    setup_data();
    vector<int> triangle1;
    vector<int> triangle2;
    while (0 < dataset.edges.size()) {
      auto edge_itr = dataset.edges.begin();
      auto shared_triangles = dataset.edge_map[*edge_itr];
      if (shared_triangles.size() < 2) {
        dataset.edges.erase(*edge_itr);
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
    dataset.set_edges();
    dataset.set_vertex_map();
    dataset.set_triangle_map();
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

      // add edges to 'edges'
      std::string edge;
      triangles_index.clear();
      for (int i = 0; i < 2; ++i) {
        if (non_shared1[0] < shared[i]) {
          edge = std::to_string(non_shared1[0]) + '-' + std::to_string(shared[i]);
          if (dataset.boundary_edges.find(edge) == dataset.boundary_edges.end())
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
          if (dataset.boundary_edges.find(edge) == dataset.boundary_edges.end())
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
          if (dataset.boundary_edges.find(edge) == dataset.boundary_edges.end())
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
          if (dataset.boundary_edges.find(edge) == dataset.boundary_edges.end())
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

#endif // VALIDATOR_DELAUNAY_HPP