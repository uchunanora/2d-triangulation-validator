#ifndef DATASET_HPP
#define DATASET_HPP

#include <algorithm>
#include <fstream>
#include <iostream>
#include <queue>
#include <set>
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
  explicit Dataset(std::string vertex_file, std::string triangle_file) {
    set_vertex_from_csv(vertex_file);
    set_triangle_from_csv(triangle_file);
  };

  bool has_boundary() const { return !boundaries.empty(); }

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
  // ---- 旧 check_overlap_vertexes() は廃止 ----
  // 代わりに check_degenerate_triangles() を使用する。
  // 別IDの重複座標頂点（スリット/ボウタイの接触点）は許容される。

  /**
   * 退化三角形チェック（旧 check_overlap_vertexes の代替）
   *
   * 同一三角形内に同座標の頂点がある場合（面積ゼロ）を検出して REJECT する。
   * これが唯一の「重複座標による拒否」条件。
   * それ以外の重複座標（別IDのスリット/ボウタイ接触点）は許容する。
   *
   * @return true: 問題なし, false: 退化三角形が存在
   */
  bool check_degenerate_triangles() {
    for (size_t i = 0; i < triangles.size(); ++i) {
      int a = triangles[i][0];
      int b = triangles[i][1];
      int c = triangles[i][2];
      if (vertexes[a] == vertexes[b] ||
          vertexes[b] == vertexes[c] ||
          vertexes[a] == vertexes[c]) {
        std::cout << "Degenerate triangle detected (zero area)!" << std::endl;
        std::cout << "triangle_index: " << i
                  << ", vertices: (" << a << ", " << b << ", " << c << ")" << std::endl;
        return false;
      }
    }
    return true;
  }

  /**
   * ボウタイ頂点の検出と自動ID分離
   *
   * vertex_map（各頂点が属する三角形の集合）を用いて、
   * ある頂点 v の周囲三角形が「v 以外の頂点を共有するか」で
   * エッジ連結成分に分解する。連結成分が2つ以上ならボウタイ頂点。
   * 2番目以降の連結成分に新しい頂点IDを割り当てて分離する。
   *
   * 計算量: 典型的なメッシュでは O(N)（各頂点の隣接三角形数が定数）
   * 境界情報不要: triangles のみで判定可能
   *
   * @return 分離した回数
   */
  int detect_and_split_bowtie_vertices() {
    // vertex_map を構築（各頂点 → 所属三角形の集合）
    std::unordered_map<int, std::unordered_set<int>> v_map;
    for (int i = 0; i < (int)triangles.size(); ++i) {
      for (int k = 0; k < 3; ++k) {
        v_map[triangles[i][k]].insert(i);
      }
    }

    int split_count = 0;
    int next_id = (int)vertexes.size();

    // 元の頂点IDリストを先に取得（イテレーション中の変更を避ける）
    std::vector<int> original_vertex_ids;
    original_vertex_ids.reserve(v_map.size());
    for (const auto& kv : v_map) {
      original_vertex_ids.push_back(kv.first);
    }

    for (int v : original_vertex_ids) {
      const auto& tri_set = v_map[v];
      if (tri_set.size() <= 1) continue;

      std::vector<int> tri_list(tri_set.begin(), tri_set.end());

      // BFS で連結成分分解
      // 2つの三角形が「v 以外の頂点を共有する」場合に連結とみなす
      std::vector<int> component(tri_list.size(), -1);
      int num_components = 0;

      for (size_t start = 0; start < tri_list.size(); ++start) {
        if (component[start] != -1) continue;
        component[start] = num_components;
        std::queue<size_t> q;
        q.push(start);
        while (!q.empty()) {
          size_t cur = q.front(); q.pop();
          int ti = tri_list[cur];
          // v 以外の頂点を収集
          std::set<int> cur_others;
          for (int k = 0; k < 3; ++k) {
            if (triangles[ti][k] != v) cur_others.insert(triangles[ti][k]);
          }
          for (size_t nb = 0; nb < tri_list.size(); ++nb) {
            if (component[nb] != -1) continue;
            int tj = tri_list[nb];
            bool shares_other = false;
            for (int k = 0; k < 3; ++k) {
              if (triangles[tj][k] != v && cur_others.count(triangles[tj][k])) {
                shares_other = true;
                break;
              }
            }
            if (shares_other) {
              component[nb] = num_components;
              q.push(nb);
            }
          }
        }
        num_components++;
      }

      if (num_components <= 1) continue;

      // 2番目以降の連結成分に新IDを割り当て
      for (int comp = 1; comp < num_components; ++comp) {
        int new_id = next_id++;
        vertexes.push_back(vertexes[v]); // 同座標の頂点を追加

        // 三角形の参照を新IDに書き換え
        for (size_t idx = 0; idx < tri_list.size(); ++idx) {
          if (component[idx] != comp) continue;
          int ti = tri_list[idx];
          for (int k = 0; k < 3; ++k) {
            if (triangles[ti][k] == v) {
              triangles[ti][k] = new_id;
            }
          }
        }

        // boundary 内の参照も更新（境界情報がある場合）
        if (has_boundary()) {
          for (size_t bi = 0; bi < boundaries.size(); ++bi) {
            if (boundaries[bi] == v) {
              // この boundary 位置の前後の頂点が、新IDの連結成分に属するか確認
              int prev_bi = (int)((bi + boundaries.size() - 1) % boundaries.size());
              int next_bi = (int)((bi + 1) % boundaries.size());
              bool prev_in_comp = false, next_in_comp = false;
              for (size_t idx = 0; idx < tri_list.size(); ++idx) {
                if (component[idx] != comp) continue;
                int ti = tri_list[idx];
                for (int k = 0; k < 3; ++k) {
                  if (triangles[ti][k] == boundaries[prev_bi]) prev_in_comp = true;
                  if (triangles[ti][k] == boundaries[next_bi]) next_in_comp = true;
                }
              }
              if (prev_in_comp || next_in_comp) {
                boundaries[bi] = new_id;
              }
            }
          }
        }
        split_count++;
      }
    }
    if (split_count > 0) {
      std::cout << "Bowtie vertices detected and split: " << split_count << " split(s)" << std::endl;
    }
    return split_count;
  }

  /**
   * スリット頂点の検出と自動ID分離（boundary情報ありの場合のみ）
   *
   * boundaries 配列内で同一IDが2回以上出現する頂点を検出し、
   * 2回目以降の出現に新しい頂点IDを割り当てて分離する。
   *
   * boundary情報なしの場合は何もしない（スリットは内部エッジ扱いで無害）。
   *
   * @return 分離した回数
   */
  int detect_and_split_slit_vertices() {
    if (!has_boundary()) return 0;

    std::unordered_map<int, std::vector<size_t>> id_positions;
    for (size_t i = 0; i < boundaries.size(); ++i) {
      id_positions[boundaries[i]].push_back(i);
    }

    struct SplitPlan {
      int v;
      int new_id;
      size_t boundary_pos;
      std::vector<int> target_tris;
    };
    std::vector<SplitPlan> plans;
    int next_id = (int)vertexes.size();

    for (const auto& kv : id_positions) {
      if (kv.second.size() <= 1) continue;
      int v = kv.first;
      const auto& positions = kv.second;

      // v を含む三角形を収集
      std::vector<int> v_tris;
      for (int ti = 0; ti < (int)triangles.size(); ++ti) {
        for (int k = 0; k < 3; ++k) {
          if (triangles[ti][k] == v) { v_tris.push_back(ti); break; }
        }
      }
      if (v_tris.empty()) continue;

      // 辺 v-u を共有する三角形マップ
      std::unordered_map<int, std::vector<int>> edge_to_tris;
      for (int ti : v_tris) {
        for (int k = 0; k < 3; ++k) {
          int u = triangles[ti][k];
          if (u != v) edge_to_tris[u].push_back(ti);
        }
      }

      // 各出現の prev_b / next_b
      std::vector<int> prev_bs, next_bs;
      for (size_t p : positions) {
        prev_bs.push_back(boundaries[(p + boundaries.size() - 1) % boundaries.size()]);
        next_bs.push_back(boundaries[(p + 1) % boundaries.size()]);
      }

      // 各出現の三角形を、v 周りの辺巡回で収集
      // 出現 oi: 辺 v-next_b[oi] を含む三角形から開始し、
      //          辺 v-prev_b[oi] を含む三角形まで巡回
      // 巡回方向: 三角形 T(v, a, b) で辺 v-a から辺 v-b へ進む
      //           （a が「入口辺」、b が「出口辺」）

      for (size_t occ = 1; occ < positions.size(); ++occ) {
        int pb = prev_bs[occ], nb = next_bs[occ];

        auto& nb_tris = edge_to_tris[nb];

        // 辺 v-nb を共有する三角形のうち、辺 v→nb の「左側」にあるものを選ぶ
        // （CCW 三角形で v=T[k], nb=T[(k+1)%3] となる三角形）
        // boundary パスで ...→pb→v→nb→... の方向で、
        // この出現の三角形は辺 v→nb の左側（内側）にある
        int start_tri = -1;
        if (nb_tris.size() == 1) {
          start_tri = nb_tris[0];
        } else {
          for (int ti : nb_tris) {
            for (int k = 0; k < 3; ++k) {
              if (triangles[ti][k] == v && triangles[ti][(k + 1) % 3] == nb) {
                start_tri = ti;
                break;
              }
            }
            if (start_tri != -1) break;
          }
        }

        if (start_tri == -1) continue;

        // 他の出現の boundary 辺の頂点集合（この出現の pb/nb を除く）
        std::set<int> other_boundary_verts;
        for (size_t oi2 = 0; oi2 < positions.size(); ++oi2) {
          if (oi2 == occ) continue;
          other_boundary_verts.insert(prev_bs[oi2]);
          other_boundary_verts.insert(next_bs[oi2]);
        }
        other_boundary_verts.erase(pb);
        other_boundary_verts.erase(nb);

        // start_tri から辺巡回で pb に到達するまで三角形を収集
        // 他の出現の boundary 辺に到達したら停止
        std::vector<int> target;
        target.push_back(start_tri);
        int cur_tri = start_tri;
        int prev_vertex = nb;
        int safety = (int)v_tris.size() + 2;
        while (safety-- > 0) {
          int next_vertex = -1;
          for (int k = 0; k < 3; ++k) {
            int u = triangles[cur_tri][k];
            if (u != v && u != prev_vertex) { next_vertex = u; break; }
          }
          if (next_vertex == -1 || next_vertex == pb) break;
          if (other_boundary_verts.count(next_vertex)) break;
          auto& adj = edge_to_tris[next_vertex];
          int next_tri = -1;
          for (int t : adj) {
            if (t != cur_tri) { next_tri = t; break; }
          }
          if (next_tri == -1) break;
          target.push_back(next_tri);
          prev_vertex = next_vertex;
          cur_tri = next_tri;
        }

        SplitPlan plan;
        plan.v = v;
        plan.new_id = next_id++;
        plan.boundary_pos = positions[occ];
        plan.target_tris = target;
        plans.push_back(std::move(plan));
      }
    }

    // 一括書き換え（三角形の頂点 + boundary パスの該当位置）
    for (auto& plan : plans) {
      vertexes.push_back(vertexes[plan.v]);
      for (int ti : plan.target_tris) {
        for (int k = 0; k < 3; ++k) {
          if (triangles[ti][k] == plan.v) {
            triangles[ti][k] = plan.new_id;
            break;
          }
        }
      }
      boundaries[plan.boundary_pos] = plan.new_id;
    }

    int split_count = (int)plans.size();
    if (split_count > 0) {
      std::cout << "Slit vertices detected and split: " << split_count << " split(s)" << std::endl;
    }
    return split_count;
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
    if (!has_boundary()) return;
    for (int i = 0, n = (int)boundaries.size(); i < n; ++i) {
      int a = boundaries[i];
      int b = boundaries[(i + 1) % n];
      if (a < b) {
        boundary_edges.insert(std::to_string(a) + '-' + std::to_string(b));
      } else if (a > b) {
        boundary_edges.insert(std::to_string(b) + '-' + std::to_string(a));
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
    if (has_boundary()) {
      for (int i = 0, n = (int)boundaries.size(); i < n; ++i) {
        int a = boundaries[i];
        int b = boundaries[(i + 1) % n];
        if (a < b) {
          edges.erase(std::to_string(a) + '-' + std::to_string(b));
        } else if (a > b) {
          edges.erase(std::to_string(b) + '-' + std::to_string(a));
        }
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
