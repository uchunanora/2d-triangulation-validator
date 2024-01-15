#ifndef VALIDATOR_TRIANGULATION_HPP
#define VALIDATOR_TRIANGULATION_HPP

#include <iostream>
#include <pstv/dataset.hpp>
#include <pstv/interval-tree/interval_tree.hpp>
#include <pstv/validation.hpp>
#include <utils/display.hpp>
#include <vector>

namespace pstv {
class ValidatorTriangulation {
  pstv::Validation validation;
  pstv::Dataset dataset;
  std::vector<int> polygon;
  lib_interval_tree::interval_tree<lib_interval_tree::interval<double>> interval_tree_x;
  lib_interval_tree::interval_tree<lib_interval_tree::interval<double>> interval_tree_y;

public:
  explicit ValidatorTriangulation(){};

  bool validate(pstv::Dataset ds, int init_tri_index = 0) {
    int index = 0;
    int flag = 0;
    int pre_pattern = 0;
    dataset = ds;
    setup_data();
    dataset.unprocessed_set.erase(init_tri_index);
    polygon = init_polygon(init_tri_index);
    init_interval_tree();
    while (1) {
      index = index % polygon.size();
      if (check_boundary() == 0) {
        return true;
      }
      int shared_vertex_index_a = polygon[index % polygon.size()];
      int shared_vertex_index_b = polygon[(index + 1) % polygon.size()];
      int next_triangle_index = search_next_triangle(shared_vertex_index_a, shared_vertex_index_b);
      if (next_triangle_index == -1) {
        index++;
        flag++;
        continue;
      } else {
        if (pre_pattern != 1) {
          index++;
        }
        int vertex_index_c = get_another_vertex_index_from_triangles(shared_vertex_index_a, shared_vertex_index_b, next_triangle_index);
        if (_validate(pre_pattern, shared_vertex_index_a, shared_vertex_index_b, vertex_index_c, next_triangle_index)) {
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

  bool _validate(int &pre_pattern, int vertex_index_a, int vertex_index_b, int vertex_index_c, int triangle_index) {
    auto ret_a = std::find(polygon.begin(), polygon.end(), vertex_index_a);
    int index_a_of_polygon = std::distance(polygon.begin(), ret_a);
    auto ret_b = std::find(polygon.begin(), polygon.end(), vertex_index_b);
    int index_b_of_polygon = std::distance(polygon.begin(), ret_b);
    auto ret_c = std::find(polygon.begin(), polygon.end(), vertex_index_c);
    int index_c_of_polygon = std::distance(polygon.begin(), ret_c);
    std::vector<int> candidates_pol;
    std::vector<std::vector<std::vector<double>>> candidates;
    if (ret_c != polygon.end()) {
      if ((index_b_of_polygon + 1) % polygon.size() == index_c_of_polygon % polygon.size()) {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
          remove_interval_tree_x(vertex_index_a, vertex_index_b);
          remove_interval_tree_y(vertex_index_a, vertex_index_b);
          candidates = get_overlap_candidates(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]);
          for (const auto &candidate : candidates) {
            if (validation.intersection(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate[0], candidate[1]) == -1) {
              return true;
            }
          }
          remove_interval_tree_x(vertex_index_b, vertex_index_c);
          remove_interval_tree_y(vertex_index_b, vertex_index_c);
          insert_interval_tree_x(vertex_index_a, vertex_index_c);
          insert_interval_tree_y(vertex_index_a, vertex_index_c);
          polygon.erase(ret_b);
          dataset.unprocessed_set.erase(triangle_index);
          pre_pattern = 1;
          return false; // share 2 sides and 3 points (a,b,c)
        } else {
          std::cout << "Triangulation is NOT verified !" << std::endl;
          std::cout << "Error: Orientation" << std::endl;
          return true;
        }
      } else {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
          pre_pattern = 4;
          return false; // skip: share 1 side and 3 points or 2 sides and 3 points (c,a,b)
        }
      }
    } else {
      if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]) > 0) {
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_a_of_polygon + polygon.size() - 1) % polygon.size()]]) > 0) {
          if (validation.orientation(dataset.vertexes[polygon[(index_a_of_polygon + polygon.size() - 1) % polygon.size()]], dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]) < 0) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Orientation" << std::endl;
            return true;
          }
        }
        if (validation.orientation(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_b_of_polygon + 1) % polygon.size()]]) > 0) {
          if (validation.orientation(dataset.vertexes[vertex_index_b], dataset.vertexes[polygon[(index_b_of_polygon + 1) % polygon.size()]], dataset.vertexes[vertex_index_c]) < 0) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Orientation" << std::endl;
            return true;
          }
        }
        remove_interval_tree_x(vertex_index_a, vertex_index_b);
        remove_interval_tree_y(vertex_index_a, vertex_index_b);
        candidates = get_overlap_candidates(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c]);
        for (const auto &candidate : candidates) {
          if (validation.intersection(dataset.vertexes[vertex_index_a], dataset.vertexes[vertex_index_c], candidate[0], candidate[1]) == -1) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Intersection" << std::endl;
            return true;
          }
        }
        candidates = get_overlap_candidates(dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c]);
        for (const auto &candidate : candidates) {
          if (validation.intersection(dataset.vertexes[vertex_index_b], dataset.vertexes[vertex_index_c], candidate[0], candidate[1]) == -1) {
            std::cout << "Triangulation is NOT verified !" << std::endl;
            std::cout << "Error: Intersection" << std::endl;
            return true;
          }
        }
        insert_interval_tree_x(vertex_index_a, vertex_index_c);
        insert_interval_tree_x(vertex_index_b, vertex_index_c);
        insert_interval_tree_y(vertex_index_a, vertex_index_c);
        insert_interval_tree_y(vertex_index_b, vertex_index_c);
        std::vector<int>::iterator itr;
        itr = std::find(polygon.begin(), polygon.end(), vertex_index_a);
        polygon.insert(itr + 1, vertex_index_c);
        dataset.unprocessed_set.erase(triangle_index);
        pre_pattern = 2;
        return false; // share 1 side and 2 points
      } else {
        std::cout << "Triangulation is NOT verified !" << std::endl;
        std::cout << "Error: Orientation" << std::endl;
        return true;
      }
    }
    return 0;
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

  void init_interval_tree() {
    int polygon_size = polygon.size();
    for (size_t i = 0; i < polygon.size(); ++i) {
      interval_tree_x.insert(lib_interval_tree::make_safe_interval<double>(dataset.vertexes[polygon[i % polygon_size]][0], dataset.vertexes[polygon[(i + 1) % polygon_size]][0]), {dataset.vertexes[polygon[i % polygon_size]], dataset.vertexes[polygon[(i + 1) % polygon_size]]});
      interval_tree_y.insert(lib_interval_tree::make_safe_interval<double>(dataset.vertexes[polygon[i % polygon_size]][1], dataset.vertexes[polygon[(i + 1) % polygon_size]][1]), {dataset.vertexes[polygon[i % polygon_size]], dataset.vertexes[polygon[(i + 1) % polygon_size]]});
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
      std::vector<int> boundary_for_check = dataset.boundaries;
      std::vector<int> polygon_for_check = polygon;
      if (boundary_for_check.size() == polygon.size()) {
        auto ret = std::find(boundary_for_check.begin(), boundary_for_check.end(), polygon[0]);
        int start = std::distance(boundary_for_check.begin(), ret);
        std::sort(boundary_for_check.begin(), boundary_for_check.end());
        std::sort(polygon_for_check.begin(), polygon_for_check.end());
        if (boundary_for_check == polygon_for_check) {
          std::cout << "Triangulation is verified !" << std::endl;
          return 0; // Triangulation Verification
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
    return -1; // don't have next triangle.
  }

  int get_another_vertex_index_from_triangles(int vertex_index_a, int vertex_index_b, int triangle_index) {
    std::vector<int> vector_vertex_index_c = dataset.triangles[triangle_index];
    vector_vertex_index_c.erase(std::remove(vector_vertex_index_c.begin(), vector_vertex_index_c.end(), vertex_index_a), vector_vertex_index_c.end());
    vector_vertex_index_c.erase(std::remove(vector_vertex_index_c.begin(), vector_vertex_index_c.end(), vertex_index_b), vector_vertex_index_c.end());
    int vertex_index_c = vector_vertex_index_c[0];
    return vertex_index_c;
  }

  std::vector<std::vector<std::vector<double>>> get_overlap_candidates(std::vector<double> vertex1, std::vector<double> vertex2) {
    std::vector<std::vector<std::vector<double>>> overlap_intervals_x;
    std::vector<std::vector<std::vector<double>>> overlap_intervals_y;
    std::vector<std::vector<std::vector<double>>> candidates;
    overlap_intervals_x = get_overlap_intervals_x(vertex1, vertex2);
    overlap_intervals_y = get_overlap_intervals_y(vertex1, vertex2);
    for (int i = 0; i < overlap_intervals_x.size(); i++) {
      for (int j = 0; j < overlap_intervals_y.size(); j++) {
        if (overlap_intervals_x[i] == overlap_intervals_y[j]) {
          candidates.push_back(overlap_intervals_x[i]);
          break;
        }
      }
    }
    return candidates;
  }

  bool vv_equal(std::vector<double> p1, std::vector<double> p2, std::vector<double> q1, std::vector<double> q2) {
    if (p1 == q1) {
      if (p2 == q2) {
        return true;
      }
    } else if (p1 == q2) {
      if (p2 == q1) {
        return true;
      }
    }
    return false;
  }

  void insert_interval_tree_x(int vertex_index_1, int vertex_index_2) {
    interval_tree_x.insert(
        lib_interval_tree::make_safe_interval<double>(dataset.vertexes[vertex_index_1][0], dataset.vertexes[vertex_index_2][0]),
        {{dataset.vertexes[vertex_index_1][0], dataset.vertexes[vertex_index_1][1]}, {dataset.vertexes[vertex_index_2][0], dataset.vertexes[vertex_index_2][1]}});
  }

  void insert_interval_tree_y(int vertex_index_1, int vertex_index_2) {
    interval_tree_y.insert(
        lib_interval_tree::make_safe_interval<double>(dataset.vertexes[vertex_index_1][0], dataset.vertexes[vertex_index_2][0]),
        {{dataset.vertexes[vertex_index_1][0], dataset.vertexes[vertex_index_1][1]}, {dataset.vertexes[vertex_index_2][0], dataset.vertexes[vertex_index_2][1]}});
  }

  void remove_interval_tree_x(int vertex_index_1, int vertex_index_2) {
    interval_tree_x.find_all(lib_interval_tree::make_safe_interval<double>(dataset.vertexes[vertex_index_1][0], dataset.vertexes[vertex_index_2][0]), [this, &vertex_index_1, &vertex_index_2](auto const &iter) {
      if (vv_equal(dataset.vertexes[vertex_index_1], dataset.vertexes[vertex_index_2], iter.segment()[0], iter.segment()[1])) {
        interval_tree_x.erase(iter);
        return false;
      }
      return true;
    });
  }

  void remove_interval_tree_y(int vertex_index_1, int vertex_index_2) {
    interval_tree_y.find_all(lib_interval_tree::make_safe_interval<double>(dataset.vertexes[vertex_index_1][1], dataset.vertexes[vertex_index_2][1]), [this, &vertex_index_1, &vertex_index_2](auto const &iter) {
      if (vv_equal(dataset.vertexes[vertex_index_1], dataset.vertexes[vertex_index_2], iter.segment()[0], iter.segment()[1])) {
        interval_tree_y.erase(iter);
        return false;
      }
      return true;
    });
  }

  std::vector<std::vector<std::vector<double>>> get_overlap_intervals_x(std::vector<double> vertex1, std::vector<double> vertex2) {
    std::vector<std::vector<std::vector<double>>> overlap_intervals;
    interval_tree_x.overlap_find_all(lib_interval_tree::make_safe_interval<double>(vertex1[0], vertex2[0]), [&overlap_intervals](auto iter) {
      overlap_intervals.push_back(iter.segment());
      return true;
    });
    return overlap_intervals;
  }

  std::vector<std::vector<std::vector<double>>> get_overlap_intervals_y(std::vector<double> vertex1, std::vector<double> vertex2) {
    std::vector<std::vector<std::vector<double>>> overlap_intervals;
    interval_tree_y.overlap_find_all(lib_interval_tree::make_safe_interval<double>(vertex1[1], vertex2[1]), [&overlap_intervals](auto iter) {
      overlap_intervals.push_back(iter.segment());
      return true;
    });
    return overlap_intervals;
  }
};
} // namespace pstv

#endif