#ifndef DATASET_HPP
#define DATASET_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using std::vector;

namespace pstv {
class Dataset {
public:
  vector<vector<double>> vertexes;
  vector<vector<int>> triangles;
  vector<int> boundaries;
  // Trinagulation
  std::unordered_map<std::string, std::unordered_set<int>> edge_map;
  std::unordered_set<std::string> boundary_edges;
  std::unordered_set<std::string> edges;
  std::unordered_set<int> unprocessed_set;
  // Delaunay
  std::unordered_map<std::string, std::unordered_set<int>> vertex_map;
  vector<std::unordered_set<int>> triangle_map;

  explicit Dataset(){};
  explicit Dataset(std::string vertex_file, std::string triangle_file, std::string boundary_file) {
    set_vertex_from_csv(vertex_file);
    set_triangle_from_csv(triangle_file);
    set_boundary_from_csv(boundary_file);
  };

private:
  template <class _T>
  vector<vector<_T>> get_data(std::string filename) {
    std::ifstream ifs(filename);
    vector<vector<_T>> vv;
    for (std::string value; std::getline(ifs, value);) {
      vv.push_back(vector<_T>());
      for (std::stringstream ss(value); std::getline(ss, value, ',');) {
        _T num;
        if (std::is_integral<_T>::value) {
          num = atoi(value.c_str());
        } else {
          num = atof(value.c_str());
        }
        vv[vv.size() - 1].push_back(num);
      }
    }
    return vv;
  }
  void set_vertex_from_csv(std::string filename) {
    vertexes = get_data<double>(filename);
  }
  void set_triangle_from_csv(std::string filename) {
    triangles = get_data<int>(filename);
  }
  void set_boundary_from_csv(std::string filename) {
    vector<vector<int>> boundary_tmp = get_data<int>(filename);
    for (size_t i = 0; i < boundary_tmp.size(); ++i) {
      boundaries.push_back(boundary_tmp[i][0]);
    }
  }

public:
  bool check_overlap_vertexes() {
#pragma omp parallel for
    for (size_t i = 0; i < vertexes.size() - 1; ++i) {
      for (size_t j = i + 1; j < vertexes.size(); ++j) {
        if (vertexes[i] == vertexes[j]) {
          std::cout << "Vertexes Overlap!" << std::endl;
          std::cout << "index_set: (" << i << ", " << j << ")" << std::endl;
          return false; // overlap
        }
      }
    }
    return true;
  }

  void set_vertex_map() {
    vector<int> triangle;
    std::unordered_set<int> tmp;
    std::string key;
    for (int i = 0; i < triangles.size(); i++) {
      triangle = triangles[i];
      key = std::to_string(triangle[0]);
      if (vertex_map[key].size() == 0) {
        vertex_map[key] = {i};
      } else {
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
      key = std::to_string(triangle[1]);
      if (vertex_map[key].size() == 0) {
        vertex_map[key] = {i};
      } else {
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
      key = std::to_string(triangle[2]);
      if (vertex_map[key].size() == 0) {
        vertex_map[key] = {i};
      } else {
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
    }
  }

  void set_triangle_map() {
    vector<int> triangle;
    vector<int> _triangle;
    std::string key;
    std::string key_for_vertex_map;
    for (int i = 0, n = triangles.size(); i < n; ++i) {
      triangle = triangles[i];
      std::unordered_set<int> set;
      // set.clear();
      key_for_vertex_map = std::to_string(triangle[0]);
      for (const auto &triangle_index : vertex_map[key_for_vertex_map]) {
        _triangle = triangles[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      key_for_vertex_map = std::to_string(triangle[1]);
      for (const auto &triangle_index : vertex_map[key_for_vertex_map]) {
        _triangle = triangles[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      key_for_vertex_map = std::to_string(triangle[2]);
      for (const auto &triangle_index : vertex_map[key_for_vertex_map]) {
        _triangle = triangles[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      set.erase(triangle[0]);
      set.erase(triangle[1]);
      set.erase(triangle[2]);
      triangle_map.push_back(set);
    }
  }

  void set_edge_map() {
    vector<int> triangle;
    std::unordered_set<int> tmp;
    std::string key;
    // #pragma omp parallel for
    for (int i = 0, n = triangles.size(); i < n; ++i) {
      triangle = triangles[i];
      std::sort(triangle.begin(), triangle.end());
      key = std::to_string(triangle[0]) + '-' + std::to_string(triangle[1]);
      if (edge_map[key].size() == 0) {
        edge_map[key] = {i};
      } else {
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
      key = std::to_string(triangle[1]) + '-' + std::to_string(triangle[2]);
      if (edge_map[key].size() == 0) {
        edge_map[key] = {i};
      } else {
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
      key = std::to_string(triangle[0]) + '-' + std::to_string(triangle[2]);
      if (edge_map[key].size() == 0) {
        edge_map[key] = {i};
      } else {
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
    }
  }
  void set_boundary_edges() {
    for (int i = 0, n = boundaries.size(); i < n; ++i) {
      if (boundaries[i % n] < boundaries[(i + 1) % n]) {
        boundary_edges.insert(std::to_string(boundaries[i % n]) + '-' + std::to_string(boundaries[(i + 1) % n]));
      } else {
        boundary_edges.insert(std::to_string(boundaries[(i + 1) % n]) + '-' + std::to_string(boundaries[i % n]));
      }
    }
  }
  void set_edges() {
    vector<int> triangle;
    std::string edge;
    // #pragma omp parallel for
    for (int i = 0, n = triangles.size(); i < n; ++i) {
      triangle = triangles[i];
      std::sort(triangle.begin(), triangle.end());
      edge = std::to_string(triangle[0]) + '-' + std::to_string(triangle[1]);
      edges.insert(edge);
      edge = std::to_string(triangle[1]) + '-' + std::to_string(triangle[2]);
      edges.insert(edge);
      edge = std::to_string(triangle[0]) + '-' + std::to_string(triangle[2]);
      edges.insert(edge);
    }
    for (int i = 0, n = boundaries.size(); i < n; ++i) {
      if (boundaries[i % n] < boundaries[(i + 1) % n]) {
        edges.erase(std::to_string(boundaries[i % n]) + '-' + std::to_string(boundaries[(i + 1) % n]));
      } else {
        edges.erase(std::to_string(boundaries[(i + 1) % n]) + '-' + std::to_string(boundaries[i % n]));
      }
    }
  }
  void set_unprocessed_set() {
    for (int i = 0, n = triangles.size(); i < n; ++i) {
      unprocessed_set.insert(i);
    }
  }
};
} // namespace pstv
#endif // DATASET_HPP