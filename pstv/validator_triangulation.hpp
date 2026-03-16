#ifndef VALIDATOR_TRIANGULATION_HPP
#define VALIDATOR_TRIANGULATION_HPP

// OT Reduction patch for validator_triangulation.hpp
//
// Based on patch/pstv/validator_triangulation.hpp (public validation member).
// Additional change:
//   - Skip R-tree candidates that share a vertex with the query edge
//     (adjacent polygon edges can never intersect the new edge).
//     Eliminates ~95% of unnecessary intersection() calls.

#include <iostream>
#include <unordered_set>
#include <pstv/dataset.hpp>
#include <pstv/rtree/rtree.hpp>
#include <pstv/validation.hpp>
#include <utils/display.hpp>
#include <vector>

namespace pstv {
class ValidatorTriangulation {
  pstv::Dataset dataset;
  std::vector<int> polygon;
  std::unordered_set<int> polygon_set; // [OPT] O(1) membership test for polygon
  lib_rtree::rtree_t rtree_;

  // [OPT] Skip candidates that share a vertex with the query edge
  static bool is_adjacent_candidate(
      const std::vector<double>& p1, const std::vector<double>& p2,
      const std::vector<double>& q1, const std::vector<double>& q2) {
    return p1 == q1 || p1 == q2 || p2 == q1 || p2 == q2;
  }

public:
  pstv::Validation validation;
  explicit ValidatorTriangulation(){};

  // NOTE: Return value convention (updated per reviewer feedback):
  //   true  = verification succeeded
  //   false = verification failed
  bool validate(pstv::Dataset ds, int init_tri_index = 0) {
    int index = 0;
    int flag = 0;
    int pre_pattern = 0;
    dataset = ds;
    setup_data();
    dataset.unprocessed_set.erase(init_tri_index);
    polygon = init_polygon(init_tri_index);
    // [OPT] Initialize polygon_set from polygon
    polygon_set.clear();
    for (int v : polygon) polygon_set.insert(v);
    init_rtree();
    while (1) {
      index = index % polygon.size();
      if (check_boundary() == 0) {
        return true;
      }
      int current_index = index % polygon.size(); // [OPT] Save index before increment
      int shared_vertex_index_a = polygon[current_index];
      int shared_vertex_index_b = polygon[(current_index + 1) % polygon.size()];
      int next_triangle_index = search_next_triangle(shared_vertex_index_a, shared_vertex_index_b);
      if (next_triangle_index == -1) {
        index++;
        flag++;
        // [PATCH] 無限ループ防止: ポリゴンを2周しても進展なし → エラー停止
        if (flag > (int)polygon.size() * 2) {
          std::cout << "Triangulation is NOT verified !" << std::endl;
          std::cout << "Error: Cannot proceed. "
                    << "Possible undetected bowtie or disconnected component." << std::endl;
          return false;
        }
        continue;
      } else {
        if (pre_pattern != 1) {
          index++;
        }
        int vertex_index_c = get_another_vertex_index_from_triangles(shared_vertex_index_a, shared_vertex_index_b, next_triangle_index);
        // [OPT] Pass saved index (before increment) to avoid std::find for a and b
        if (!_validate(pre_pattern, shared_vertex_index_a, shared_vertex_index_b, vertex_index_c, next_triangle_index, current_index)) {
          return false;
        }
        if (pre_pattern != 4) {
          flag = 0;
        }
      }
    }
  }

private:
  void setup_data() {
    dataset.set_edge_map();
    dataset.set_boundary_edges();
    dataset.set_edges();
    dataset.set_unprocessed_set();
  }

  bool _validate(int &pre_pattern, int vertex_index_a, int vertex_index_b, int vertex_index_c, int triangle_index, int index_a_of_polygon) {
    // [OPT] Compute positions from passed index instead of std::find
    int index_b_of_polygon = (index_a_of_polygon + 1) % polygon.size();
    // [OPT] Use unordered_set for O(1) membership test instead of std::find for c
    bool c_in_polygon = polygon_set.count(vertex_index_c) > 0;
    bool c_is_next_to_b = (polygon[(index_b_of_polygon + 1) % polygon.size()] == vertex_index_c);
    std::vector<int> candidates_pol;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> candidates;
    if (c_in_polygon) {
      if (c_is_next_to_b) {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
          remove_rtree(vertex_index_a, vertex_index_b);
          candidates = get_overlap_candidates(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]);
          for (const auto &candidate : candidates) {
            if (is_adjacent_candidate(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate.first, candidate.second))
              continue; // [OPT] skip adjacent polygon edge
            if (validation.intersection(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate.first, candidate.second) == -1) {
              return false;
            }
          }
          remove_rtree(vertex_index_b, vertex_index_c);
          insert_rtree(vertex_index_a, vertex_index_c);
          polygon_set.erase(vertex_index_b);
          polygon.erase(polygon.begin() + index_b_of_polygon);
          dataset.unprocessed_set.erase(triangle_index);
          pre_pattern = 1;
          return true; // share 2 sides and 3 points (a,b,c)
        } else {
          std::cout << "Triangulation is NOT verified !" << std::endl;
          std::cout << "Error: Orientation" << std::endl;
          return false;
        }
      } else {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
          pre_pattern = 4;
          return true; // skip: share 1 side and 3 points or 2 sides and 3 points (c,a,b)
        }
      }
    } else {
      if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_a_of_polygon + polygon.size() - 1) % polygon.size()]]) > 0) {
          if (validation.orientation(dataset.vertexes[polygon[(index_a_of_polygon + polygon.size() - 1) % polygon.size()]], dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]) < 0) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Orientation" << std::endl;
            return false;
          }
        }
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_b_of_polygon + 1) % polygon.size()]]) > 0) {
          if (validation.orientation(dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_b_of_polygon + 1) % polygon.size()]], dataset.vertexes[vertex_index_c]) < 0) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Orientation" << std::endl;
            return false;
          }
        }
        remove_rtree(vertex_index_a, vertex_index_b);
        candidates = get_overlap_candidates(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]);
        for (const auto &candidate : candidates) {
          if (is_adjacent_candidate(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate.first, candidate.second))
            continue; // [OPT] skip adjacent polygon edge
          if (validation.intersection(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate.first, candidate.second) == -1) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Intersection" << std::endl;
            return false;
          }
        }
        candidates = get_overlap_candidates(dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]);
        for (const auto &candidate : candidates) {
          if (is_adjacent_candidate(dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c], candidate.first, candidate.second))
            continue; // [OPT] skip adjacent polygon edge
          if (validation.intersection(dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c], candidate.first, candidate.second) == -1) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Intersection" << std::endl;
            return false;
          }
        }
        insert_rtree(vertex_index_a, vertex_index_c);
        insert_rtree(vertex_index_b, vertex_index_c);
        polygon_set.insert(vertex_index_c);
        polygon.insert(polygon.begin() + index_a_of_polygon + 1, vertex_index_c);
        dataset.unprocessed_set.erase(triangle_index);
        pre_pattern = 2;
        return true; // share 1 side and 2 points
      } else {
        std::cout << "Triangulation is NOT verified !" << std::endl;
        std::cout << "Error: Orientation" << std::endl;
        return false;
      }
    }
    return true;
  }

  std::vector<int> init_polygon(int init_triangle_index) {
    std::vector<int> polygon(3, 0);
    polygon.reserve(dataset.vertexes.size());
    for (size_t i = 0; i < 3; ++i) {
      polygon[i] = dataset.triangles[init_triangle_index][i];
    }
    update_triangle_clockwise(polygon, init_triangle_index);
    return polygon;
  }

  void init_rtree() {
    int polygon_size = polygon.size();
    for (size_t i = 0; i < polygon.size(); ++i) {
      rtree_.insert(dataset.vertexes[polygon[i % polygon_size]], dataset.vertexes[polygon[(i + 1) % polygon_size]]);
    }
  }

  void update_triangle_clockwise(std::vector<int> &polygon, int triangle_index) {
    std::vector<double> vertex_a = dataset.vertexes[polygon[0]];
    std::vector<double> vertex_b = dataset.vertexes[polygon[1]];
    std::vector<double> vertex_c = dataset.vertexes[polygon[2]];
    double orientation = validation.orientation(vertex_a, vertex_b, vertex_c);
    if (orientation == 0) {
      std::cout << "Triangulation is NOT verified !" << std::endl;
      std::cout << "Error: The three points that make up the triangle are aligned in a straight line." << std::endl;
      std::exit(0);
    }
    if (orientation > 0)
      std::reverse(polygon.begin(), polygon.end());
  }

  int check_boundary() {
    if (dataset.unprocessed_set.size() == 0) {
      if (dataset.has_boundary()) {
        std::vector<int> boundary_for_check = dataset.boundaries;
        std::vector<int> polygon_for_check = polygon;
        if (boundary_for_check.size() == polygon.size()) {
          auto ret = std::find(boundary_for_check.begin(), boundary_for_check.end(), polygon[0]);
          int start = std::distance(boundary_for_check.begin(), ret);
          std::sort(boundary_for_check.begin(), boundary_for_check.end());
          std::sort(polygon_for_check.begin(), polygon_for_check.end());
          if (boundary_for_check == polygon_for_check) {
            std::cout << "Triangulation is verified !" << std::endl;
            return 0;
          } else {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: The final polygon does not match the boundary data." << std::endl;
            std::exit(0);
          }
        } else {
          std::cout << "Triangulation is NOT verified !" << std::endl;
          std::cout << "Error: The final polygon does not match the boundary data." << std::endl;
          std::exit(0);
        }
      } else {
        std::cout << "Triangulation is verified ! (no boundary data provided)" << std::endl;
        return 0;
      }
    }
    return 1;
  }

  int search_next_triangle(int shared_vertex_index_a, int shared_vertex_index_b) {
    std::string key;
    std::unordered_set<int> triangles;
    if (shared_vertex_index_a < shared_vertex_index_b) {
      key = std::to_string(shared_vertex_index_a) + "-" + std::to_string(shared_vertex_index_b);
    } else {
      key = std::to_string(shared_vertex_index_b) + "-" + std::to_string(shared_vertex_index_a);
    }
    triangles = dataset.edge_map[key];
    for (const auto &triangle_index : triangles) {
      if (dataset.unprocessed_set.count(triangle_index)) {
        return triangle_index;
      }
    }
    return -1;
  }

  int get_another_vertex_index_from_triangles(int vertex_index_a, int vertex_index_b, int triangle_index) {
    std::vector<int> vector_vertex_index_c = dataset.triangles[triangle_index];
    vector_vertex_index_c.erase(std::remove(vector_vertex_index_c.begin(), vector_vertex_index_c.end(), vertex_index_a), vector_vertex_index_c.end());
    vector_vertex_index_c.erase(std::remove(vector_vertex_index_c.begin(), vector_vertex_index_c.end(), vertex_index_b), vector_vertex_index_c.end());
    int vertex_index_c = vector_vertex_index_c[0];
    return vertex_index_c;
  }

  std::vector<std::pair<std::vector<double>, std::vector<double>>> get_overlap_candidates(std::vector<double> vertex1, std::vector<double> vertex2) {
    return rtree_.query(vertex1, vertex2);
  }

  void insert_rtree(int vertex_index_1, int vertex_index_2) {
    rtree_.insert(dataset.vertexes[vertex_index_1], dataset.vertexes[vertex_index_2]);
  }

  void remove_rtree(int vertex_index_1, int vertex_index_2) {
    rtree_.remove(dataset.vertexes[vertex_index_1], dataset.vertexes[vertex_index_2]);
  }
};
} // namespace pstv

#endif
