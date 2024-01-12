#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <set>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <gmpxx.h>
#include <interval-tree/interval_tree.hpp>

using std::cout; using std::endl; using std::vector; using std::pair;
using namespace lib_interval_tree;

//for filter
double u,theta,const6u,u_n,iccerrboundA;

//COUNT
int cnt_checkConvexHull = 0;
int cnt_searchNextTriangle = 0;
int cnt_orientation = 0;
int cnt_orientation_gmp = 0;
int cnt_nonDelaunay = 0;
int cnt_incircle = 0;

int cnt_tree_insert = 0;
int cnt_tree_erase = 0;
int cnt_tree_search = 0;
int CNT = 0;

#pragma region Display
void disp_delaunay(){
  cout << "△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△" << endl;
  cout << endl;
  cout << "       /  ____/  /          /   / /     /     /   / /" << endl;
  cout << "   /  /    _/   /       /  /  _/ /  /  /  /  /___  / " << endl;
  cout << "_____/______/_______/__/__/_____/__/__/__/__/_____/  " << endl;
  cout << endl;
  cout << "▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽△▽" << endl;
}
template<class T> void disp(std::string var_name, T &var){
  cout << var_name << ": " << var << endl;
}
template<class T> void disp_v(vector<T> &vec){
  // cout << std::setprecision(16);
  cout << "{ ";
  for(const auto &item : vec){
    cout << item << " ";
  }
  cout << " }" << endl;
}
template<class T> void disp_vv(vector<vector<T>> &vv){
  // cout << std::setprecision(16);
  cout << "{ \n";
  for(const auto &v : vv){
    for(const auto &item : v){
      cout << item << " ";
    }
    cout << "\n";
  }
  cout << " }" << endl;
}
void disp_time(timespec start,timespec end,std::string name){
  printf("%s time = ",name.c_str());
  if (end.tv_nsec < start.tv_nsec) {
    printf("%5ld.%09ld", end.tv_sec - start.tv_sec - 1,
           end.tv_nsec + (long int)1.0e+9 - start.tv_nsec);
  } else {
    printf("%5ld.%09ld", end.tv_sec - start.tv_sec,
           end.tv_nsec - start.tv_nsec);
  }
  printf("(sec)\n");
}
#pragma endregion Display

namespace pstv{
  vector<vector<double>> vert;
  vector<vector<int>> tri;
  vector<int> boundary;
  std::unordered_map<std::string,std::unordered_set<int>> edge_map;
  std::unordered_set<std::string> boundary_edges;
  std::unordered_set<std::string> edges; //外周辺以外の辺
  std::unordered_set<int> unprocessed_set; //triのindex群
  std::unordered_map<std::string,std::unordered_set<int>> vertex_map; //vertexを含む三角形群
  vector<std::unordered_set<int>> triangle_map; //三角形の周囲にあるvertex群の配列
  vector<vector<double>> x_sorted_vert;
  vector<vector<double>> y_sorted_vert;
  vector<pair<vector<double>, vector<double>>> edge_intervals; //辺のx座標の区間
  interval_tree<interval<double>> intervalTreeX;
  interval_tree<interval<double>> intervalTreeY;
  size_t tree_flag = 0;
  int hyper_polygon_size;
  int hyper_end_polygon_size;

#pragma region GetData
  template<class _T> vector<vector<_T>> get_data(std::string filename){
    //intかそれ以外かでしか対応できていない...
    std::ifstream ifs(filename);
    vector<vector<_T>> vv;
    for (std::string value; std::getline(ifs, value);) {
      vv.push_back(vector<_T>());
      for (std::stringstream ss(value); std::getline(ss, value, ',');) {
        _T num;
        if(std::is_integral<_T>::value){
          num = atoi(value.c_str());
        }else{
          num = atof(value.c_str());
          // cout << std::setprecision(16);
        }
        vv[vv.size()-1].push_back(num);
      }
    }
    return vv;
  }
  std::unordered_map<std::string,std::unordered_set<int>> get_edge_map(){
    vector<int> triangle;
    std::unordered_set<int> tmp;
    std::string key;
    // #pragma omp parallel for
    for(int i=0,n=tri.size();i<n;++i){
      triangle = tri[i];
      std::sort(triangle.begin(),triangle.end());
      key = std::to_string(triangle[0])+'-'+std::to_string(triangle[1]);
      if(edge_map[key].size() == 0){
        edge_map[key] = {i};
      }else{
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
      key = std::to_string(triangle[1])+'-'+std::to_string(triangle[2]);
      if(edge_map[key].size() == 0){
        edge_map[key] = {i};
      }else{
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
      key = std::to_string(triangle[0])+'-'+std::to_string(triangle[2]);
      if(edge_map[key].size() == 0){
        edge_map[key] = {i};
      }else{
        tmp = edge_map[key];
        tmp.insert(i);
        edge_map[key] = tmp;
      }
    }
    return edge_map;
  }

  std::unordered_set<std::string> get_boundary_edges(){
    for(int i=0,n=boundary.size();i<n;++i){
      if(boundary[i%n] < boundary[(i+1)%n]){
        boundary_edges.insert(std::to_string(boundary[i%n])+'-'+std::to_string(boundary[(i+1)%n]));
      }else{
        boundary_edges.insert(std::to_string(boundary[(i+1)%n])+'-'+std::to_string(boundary[i%n]));
      }
    }
    return boundary_edges;
  }

  std::unordered_set<std::string> get_edges(){
    vector<int> triangle;
    std::string edge;
    // #pragma omp parallel for
    for(int i=0,n=tri.size();i<n;++i){
      triangle = tri[i];
      std::sort(triangle.begin(),triangle.end());
      edge = std::to_string(triangle[0])+'-'+std::to_string(triangle[1]);
      edges.insert(edge);
      edge = std::to_string(triangle[1])+'-'+std::to_string(triangle[2]);
      edges.insert(edge);
      edge = std::to_string(triangle[0])+'-'+std::to_string(triangle[2]);
      edges.insert(edge);
    }
    for(int i=0,n=boundary.size();i<n;++i){
      if(boundary[i%n] < boundary[(i+1)%n]){
        edges.erase(std::to_string(boundary[i%n])+'-'+std::to_string(boundary[(i+1)%n]));
      }else{
        edges.erase(std::to_string(boundary[(i+1)%n])+'-'+std::to_string(boundary[i%n]));
      }
    }
    return edges;
  }


  std::unordered_set<int> get_unprocessed_set(){
    for(int i=0,n=tri.size();i<n;++i){
      unprocessed_set.insert(i);
    }
    return unprocessed_set;
  }
  std::unordered_map<std::string,std::unordered_set<int>> get_vertex_map(){
    vector<int> triangle;
    std::unordered_set<int> tmp;
    std::string key;
    for(int i=0,n=tri.size();i<n;++i){
      triangle = tri[i];
      key = std::to_string(triangle[0]);
      if(vertex_map[key].size() == 0){
        vertex_map[key] = {i};
      }else{
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
      key = std::to_string(triangle[1]);
      if(vertex_map[key].size() == 0){
        vertex_map[key] = {i};
      }else{
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
      key = std::to_string(triangle[2]);
      if(vertex_map[key].size() == 0){
        vertex_map[key] = {i};
      }else{
        tmp = vertex_map[key];
        tmp.insert(i);
        vertex_map[key] = tmp;
      }
    }
    return vertex_map;
  }
  vector<std::unordered_set<int>> get_triangle_map(){
    vector<int> triangle;
    vector<int> _triangle;
    std::string key;
    std::string key_for_vertex_map;
    for(int i=0,n=tri.size();i<n;++i){
      triangle = tri[i];
      std::unordered_set<int> set;
      //set.clear();
      key_for_vertex_map = std::to_string(triangle[0]);
      for(const auto& triangle_index : vertex_map[key_for_vertex_map]){
        _triangle = tri[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      key_for_vertex_map = std::to_string(triangle[1]);
      for(const auto& triangle_index : vertex_map[key_for_vertex_map]){
        _triangle = tri[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      key_for_vertex_map = std::to_string(triangle[2]);
      for(const auto& triangle_index : vertex_map[key_for_vertex_map]){
        _triangle = tri[triangle_index];
        set.insert(_triangle[0]);
        set.insert(_triangle[1]);
        set.insert(_triangle[2]);
      }
      set.erase(triangle[0]);
      set.erase(triangle[1]);
      set.erase(triangle[2]);
      triangle_map.push_back(set);
    }
    return triangle_map;
  }

  void get_xy_sort(){
    //disp_vv(vert);
    x_sorted_vert = vert; //コピーするよりvert作るときに一緒に作ったほうがはやい
    //disp_vv(x_sorted_vert);
    y_sorted_vert = vert; //コピーするよりvert作るときに一緒に作ったほうがはやい
    sort(x_sorted_vert.begin(),x_sorted_vert.end(),[](const vector<double> &alpha,const vector<double> &beta){return alpha[0] < beta[0];});
    sort(y_sorted_vert.begin(),y_sorted_vert.end(),[](const vector<double> &alpha,const vector<double> &beta){return alpha[1] < beta[1];});
  }
#pragma endregion GetData
#pragma region Tool
  void update_tri_clockwise(vector<vector<int>>& tri, int tri_index, int mode=0){
    //mode=0:clockwise sort, mode=1:unclockwise sort
    vector<double> vert_A = vert[tri[tri_index][0]];
    vector<double> vert_B = vert[tri[tri_index][1]];
    vector<double> vert_C = vert[tri[tri_index][2]];
    vector<double> vector_AB = {vert_B[0]-vert_A[0], vert_B[1]-vert_A[1]};
    vector<double> vector_AC = {vert_C[0]-vert_A[0], vert_C[1]-vert_A[1]};
    double cross_product = vector_AB[0]*vector_AC[1] - vector_AB[1]*vector_AC[0];
    if(fabs(cross_product) <= theta*(fabs(vector_AB[0]*vector_AC[1] + vector_AB[1]*vector_AC[0]))){
      mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_cross_product;
      mpq_ax = vert[tri[tri_index][0]][0];
      mpq_ay = vert[tri[tri_index][0]][1];
      mpq_bx = vert[tri[tri_index][1]][0];
      mpq_by = vert[tri[tri_index][1]][1];
      mpq_cx = vert[tri[tri_index][2]][0];
      mpq_cy = vert[tri[tri_index][2]][1];
      mpq_cross_product = (mpq_bx-mpq_ax)*(mpq_cy-mpq_ay) - (mpq_by-mpq_ay)*(mpq_cx-mpq_ax);
      int mpq_sgn = sgn(mpq_cross_product);
      if(mode == 0){
        if(mpq_sgn > 0){
          std::reverse(tri[tri_index].begin(), tri[tri_index].end());
        }
      }else if(mode == 1){
        if(mpq_sgn > 0){
          std::reverse(tri[tri_index].begin(), tri[tri_index].end());
        }
      }
    }else{
      if(mode == 0){
        if(cross_product > 0){
          std::reverse(tri[tri_index].begin(), tri[tri_index].end());
        }
      }else if(mode == 1){
        if(cross_product < 0){
          std::reverse(tri[tri_index].begin(), tri[tri_index].end());
        }
      }
    }
  }
  void flip(const int tri_index1, const int tri_index2){
    cout << "triangle1: "; disp_v(tri[tri_index1]);
    cout << "triangle2: "; disp_v(tri[tri_index2]);
    vector<int> triangle1 = tri[tri_index1];
    vector<int> triangle2 = tri[tri_index2];
    vector<int> all = triangle1;
    all.insert(all.end(), triangle2.begin(), triangle2.end());
    std::sort(all.begin(), all.end());
    all.erase(std::unique(all.begin(), all.end()), all.end());
    vector<int> non_shared1 = all;
    vector<int> non_shared2 = all;
    vector<int> shared = all;
    for(size_t i=0; i<3; ++i){
      non_shared1.erase(std::remove(non_shared1.begin(), non_shared1.end(), triangle2[i]), non_shared1.end());
    }
    for(size_t i=0; i<3; ++i){
      non_shared2.erase(std::remove(non_shared2.begin(), non_shared2.end(), triangle1[i]), non_shared2.end());
    }
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared1[0]), shared.end());
    shared.erase(std::remove(shared.begin(), shared.end(), non_shared2[0]), shared.end());

    vector<int> non_shared;
    non_shared.push_back(non_shared1[0]);
    non_shared.push_back(non_shared2[0]);
    vector<int> new_triangle1 = non_shared;
    vector<int> new_triangle2 = non_shared;
    new_triangle1.push_back(shared[0]);
    new_triangle2.push_back(shared[1]);
    std::sort(new_triangle1.begin(), new_triangle1.end());
    std::sort(new_triangle2.begin(), new_triangle2.end());
    tri[tri_index1] = new_triangle1;
    tri[tri_index2] = new_triangle2;
    cout << "new_triangle1: "; disp_v(tri[tri_index1]);
    cout << "new_triangle2: "; disp_v(tri[tri_index2]);
  }
  void output_tri_data(std::string filename){
    std::ofstream ofs(filename);
    for(size_t i=0; i<tri.size(); ++i){
      ofs << tri[i][0] << "," << tri[i][1] << "," << tri[i][2] << endl;
    }
  }

#pragma endregion Tool
#pragma region Judgements
  //check_overlap_vertもっとはやくできるはず...
  int check_overlap_vert(){
    size_t flag = 1;
    #pragma omp parallel for
    for(size_t i=0;i<vert.size()-1;++i){
      for(size_t j=i+1;j<vert.size();++j){
        if(vert[i] == vert[j]){
          cout << "index_set: (" << i << ", " << j << ")" << endl;
          flag = 0; //overlap
        }
      }
    }
    return flag;
  }
  // int check_overlap_vert(){
  //   std::unordered_set<double> set;
  //   for(size_t i=0;i<vert.size();++i){
  //     if(!set.insert(vert[i][0]).second){
  //       if(vert[i] == vert[j]){
  //         return 0;
  //       }
  //     }
  //   }
  //   return 1;
  // }

  /*
  int check_overlap_vert(){
    //cout << "check_overlap_vertex function." << endl;
    vector<double> sort1;
    vector<double> sort2;
    for(size_t i=0; i<vert.size()-1;++i){
      for(size_t j=i+1; j<vert.size();++j){
        sort1 = vert[i];
        sort2 = vert[j];
        std::sort(sort1.begin(),sort1.end());
        std::sort(sort2.begin(),sort2.end());
        if(sort1[0] == sort2[0]){
          if(sort1[1] == sort2[1]){
            return false; //overlap
          }
        }
      }
    }
    return true; //non overlap
  }
  */
  void filter(){
  	int i;
  	u = 1.;
  	for(i = 1; i <= 53; i++){
  		u /= 2;
  	}
  	theta = 3 * u;
    const6u = 1 + 6 * u;
  	u_n = 1.;
  	for(i = 1; i <= 1022; i++){
  		u_n /= 2;
  	}
  	iccerrboundA = 10 * u + 144 * u * u + u_n;
  }
  double verifyorientation(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c, int mode=0){
    //cout << "mode: " << mode << endl;
    cnt_orientation++;
    const double axby = (vert_a[0] - vert_c[0]) * (vert_b[1] - vert_c[1]);
    const double aybx = (vert_a[1] - vert_c[1]) * (vert_b[0] - vert_c[0]);
    const double det = axby - aybx;
    if(fabs(det) <= theta*(fabs(axby + aybx) + u_n)){
      // return 0;
      cnt_orientation_gmp++;
      mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_det;
      mpq_ax = vert_a[0];
      mpq_ay = vert_a[1];
      mpq_bx = vert_b[0];
      mpq_by = vert_b[1];
      mpq_cx = vert_c[0];
      mpq_cy = vert_c[1];
      // cout << "[ " << (mpq_ax - mpq_cx)*(mpq_by - mpq_cy) << ", " << (mpq_ay - mpq_cy)*(mpq_bx - mpq_cx) << " ]" << endl;
      mpq_det = (mpq_ax - mpq_cx)*(mpq_by - mpq_cy) - (mpq_ay - mpq_cy)*(mpq_bx - mpq_cx);
      int mpq_sgn = sgn(mpq_det); // -1 or 0 or 1
      // cout << "mpq_det: " << mpq_det << endl;
      return mpq_sgn;
    }
    return det;
  }
  double verifyincircle(vector<double> vertex_a, vector<double> vertex_b, vector<double> vertex_c, vector<double> vertex_d){
    //det>0:点dはC(a,b,c)外, det<0:点dはC(a,b,c)内, det=0:点dはC(a,b,c)上
    //vertex_a,b,cを反時計回りに.
    double vector_ABx = vertex_b[0]-vertex_a[0];
    double vector_ABy = vertex_b[1]-vertex_a[1];
    double vector_ACx = vertex_c[0]-vertex_a[0];
    double vector_ACy = vertex_c[1]-vertex_a[1];
    double cross_product = vector_ABx*vector_ACy - vector_ABy*vector_ACx;
    if(fabs(cross_product) <= theta*(fabs(vector_ABx*vector_ACy + vector_ABy*vector_ACx) + u_n)){
      mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_cross_product;
      mpq_ax = vertex_a[0];
      mpq_ay = vertex_a[1];
      mpq_bx = vertex_b[0];
      mpq_by = vertex_b[1];
      mpq_cx = vertex_c[0];
      mpq_cy = vertex_c[1];
      mpq_cross_product = (mpq_bx-mpq_ax)*(mpq_cy-mpq_ay) - (mpq_by-mpq_ay)*(mpq_cx-mpq_ax);
      int mpq_sgn = sgn(mpq_cross_product);
      if(mpq_sgn > 0){
        vector<double> tmp = vertex_b;
        vertex_b = vertex_c;
        vertex_c = tmp;
      }
    }else{
      if(cross_product > 0){
        vector<double> tmp = vertex_b;
        vertex_b = vertex_c;
        vertex_c = tmp;
      }
    }

    cnt_incircle += 1;

    double adx = vertex_a[0] - vertex_d[0];
    double bdx = vertex_b[0] - vertex_d[0];
    double cdx = vertex_c[0] - vertex_d[0];
    double ady = vertex_a[1] - vertex_d[1];
    double bdy = vertex_b[1] - vertex_d[1];
    double cdy = vertex_c[1] - vertex_d[1];

    double bdxcdy = bdx * cdy;
    double cdxbdy = cdx * bdy;
    double alift = adx * adx + ady * ady;

    double cdxady = cdx * ady;
    double adxcdy = adx * cdy;
    double blift = bdx * bdx + bdy * bdy;

    double adxbdy = adx * bdy;
    double bdxady = bdx * ady;
    double clift = cdx * cdx + cdy * cdy;

    double det = alift * (bdxcdy - cdxbdy)
    			+ blift * (cdxady - adxcdy)
    			+ clift * (adxbdy - bdxady);
    double permanent = (fabs(bdxcdy) + fabs(cdxbdy)) * alift
    				+ (fabs(cdxady) + fabs(adxcdy)) * blift
    				+ (fabs(adxbdy) + fabs(bdxady)) * clift;
    double errbound = iccerrboundA * permanent;

    if(fabs(det) <= errbound){
      mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_dx, mpq_dy;
      mpq_ax = vertex_a[0];
      mpq_ay = vertex_a[1];
      mpq_bx = vertex_b[0];
      mpq_by = vertex_b[1];
      mpq_cx = vertex_c[0];
      mpq_cy = vertex_c[1];
      mpq_dx = vertex_d[0];
      mpq_dy = vertex_d[1];
      mpq_class mpq_adx = mpq_ax - mpq_dx;
      mpq_class mpq_bdx = mpq_bx - mpq_dx;
      mpq_class mpq_cdx = mpq_cx - mpq_dx;
      mpq_class mpq_ady = mpq_ay - mpq_dy;
      mpq_class mpq_bdy = mpq_by - mpq_dy;
      mpq_class mpq_cdy = mpq_cy - mpq_dy;
      mpq_class mpq_bdxcdy = mpq_bdx * mpq_cdy;
      mpq_class mpq_cdxbdy = mpq_cdx * mpq_bdy;
      mpq_class mpq_alift = mpq_adx * mpq_adx + mpq_ady * mpq_ady;
      mpq_class mpq_cdxady = mpq_cdx * mpq_ady;
      mpq_class mpq_adxcdy = mpq_adx * mpq_cdy;
      mpq_class mpq_blift = mpq_bdx * mpq_bdx + mpq_bdy * mpq_bdy;
      mpq_class mpq_adxbdy = mpq_adx * mpq_bdy;
      mpq_class mpq_bdxady = mpq_bdx * mpq_ady;
      mpq_class mpq_clift = mpq_cdx * mpq_cdx + mpq_cdy * mpq_cdy;
      mpq_class mpq_det = mpq_alift * (mpq_bdxcdy - mpq_cdxbdy) + mpq_blift * (mpq_cdxady - mpq_adxcdy) + mpq_clift * (mpq_adxbdy - mpq_bdxady);
      int mpq_sgn = sgn(mpq_det);
      det = mpq_sgn;
    }

    return det;
  }
  int online(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c){
    //Return -1 if C is on Line[AB], otherwise return 1.
    vector<double> vector_x = {vert_a[0],vert_b[0],vert_c[0]};
    vector<double> vector_y = {vert_a[1],vert_b[1],vert_c[1]};
    std::sort(vector_x.begin(),vector_x.end());
    std::sort(vector_y.begin(),vector_y.end());
    if(verifyorientation(vert_a, vert_b, vert_c,-1) == 0){
      if(vert_a[0] == vert_b[0]){
        if(vector_y[1] == vert_c[1]){
          return -1;
        }
      }else{
        if(vector_x[1] == vert_c[0]){
          if(vert_a[1] == vert_b[1]){
            return -1;
          }else{
            if(vector_y[1] == vert_c[1]){
              return -1;
            }
          }
        }
      }
    }
    return 1;
  }
  
  // int intersection(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c, vector<double> vert_d){
  //   //Return -1 if Line-Segment[AB] and Line-Segment[CD] intersect, otherwise return 1.
  //   if(verifyorientation(vert_a,vert_b,vert_c,-2)*verifyorientation(vert_a,vert_b,vert_d,-3) < 0){
  //     if(verifyorientation(vert_c,vert_d,vert_a,-4)*verifyorientation(vert_c,vert_d,vert_b,-5) < 0){
  //       return -1;
  //     }
  //   }
  //   return 1;
  // }

  inline int onSegment(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c){
    if(std::min(vert_a[0],vert_b[0]) < vert_c[0] && vert_c[0] < std::max(vert_a[0],vert_b[0]) && std::min(vert_a[1],vert_b[1]) < vert_c[1] && vert_c[1] < std::max(vert_a[1],vert_b[1])){
      return true;
    }
    return false;
  }

  inline int intersection(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c, vector<double> vert_d){
    //Return -1 if Line-Segment[AB] and Line-Segment[CD] intersect, otherwise return 1.
    double o1,o2,o3,o4;
    o1 = verifyorientation(vert_a,vert_b,vert_c,-2);
    o2 = verifyorientation(vert_a,vert_b,vert_d,-3);
    o3 = verifyorientation(vert_c,vert_d,vert_a,-4);
    o4 = verifyorientation(vert_c,vert_d,vert_b,-5);
    if(o1*o2 < 0 && o3*o4 < 0) return -1;
    if(o1==0 && onSegment(vert_a,vert_b,vert_c)) return -1;
    if(o2==0 && onSegment(vert_a,vert_b,vert_d)) return -1;
    if(o3==0 && onSegment(vert_c,vert_d,vert_a)) return -1;
    if(o4==0 && onSegment(vert_c,vert_d,vert_b)) return -1;
    return 1;
  }


  double radius_upper(double ax, double ay, double bx, double by, double cx, double cy){
    double a2,b2,c2,a,b,c;
    double au,bu,cu;
    double v,r;
    
    //辺の長さの上限
    a2 = (bx-cx)*(bx-cx) + (by-cy)*(by-cy);
    b2 = (ax-cx)*(ax-cx) + (ay-cy)*(ay-cy);
    c2 = (ax-bx)*(ax-bx) + (ay-by)*(ay-by);
    a = std::sqrt(a2);
    b = std::sqrt(b2);
    c = std::sqrt(c2);
    au = a * const6u;
    bu = b * const6u;
    cu = c * const6u;

    double x = (1-6*u)*(a+b+c)*0.5;
    if(x <= au || x <= bu || x <= cu){
      r = std::numeric_limits<double>::infinity();
      return r;
    }
    v = x * (x - au) * (x - bu) * (x - cu);
    r = (1+16*u)*std::sqrt(a2*b2*c2/v)*0.25;
    if(!std::isfinite(r)){
      r = std::numeric_limits<double>::infinity();
      return r;
    }
    return r;
  }
#pragma endregion Judgements
#pragma region Division
  vector<int> init_polygon(int init_tri_index){
    update_tri_clockwise(tri, init_tri_index);
    vector<int> polygon(3,0);
    polygon.reserve(vert.size());//はやいのか？
    for(size_t i=0;i<3;++i){
      polygon[i] = tri[init_tri_index][i];
      // intervalTreeX.insertX({vert[tri[init_tri_index][i]][0],vert[tri[init_tri_index][i]][1]}, {vert[tri[init_tri_index][(i+1)%3]][0],vert[tri[init_tri_index][(i+1)%3]][1]});
      // intervalTreeY.insertY({vert[tri[init_tri_index][i]][0],vert[tri[init_tri_index][i]][1]}, {vert[tri[init_tri_index][(i+1)%3]][0],vert[tri[init_tri_index][(i+1)%3]][1]});
    }
    return polygon;
  }
  int check_convex_hull(vector<int> polygon){
    //操作している多角形の凸包と三角形分割全体の凸包が一致していて、未処理の三角形がないかcheck
    //コピーしてるけどポインタでする方が速くね.
    cnt_checkConvexHull++;
    if(unprocessed_set.size() == 0){
      vector<int> boundary_for_check = boundary;
      vector<int> polygon_for_check = polygon;
      if(boundary_for_check.size() == polygon.size()){
        auto ret = std::find(boundary_for_check.begin(), boundary_for_check.end(), polygon[0]);
        size_t start = std::distance(boundary_for_check.begin(), ret);
        // //順周りcheck
        // for(size_t i=0,n=boundary_for_check.size(); i<n; ++i){
        //   if(boundary_for_check[(start+i) % boundary_for_check.size()] != polygon_for_check[i]){
        //     break;
        //   }
        //   if(i == boundary_for_check.size()-1){
        //     cout << "unprocessed_set: " << unprocessed_set.size() << endl;
        //     cout << "boundary_for_check size: " << boundary_for_check.size() << endl;
        //     cout << "boundary_for_check: "; disp_v(boundary_for_check);
        //     cout << "polygon size: " << polygon.size() << endl;
        //     cout << "polygon: "; disp_v(polygon);
        //     return 0; //Division Guaranteed
        //   }
        // }
        // //逆周りcheck
        // for(size_t i=0,n=boundary_for_check.size(); i<n; ++i){
        //   if(boundary_for_check[(start-i+boundary_for_check.size()) % boundary_for_check.size()] != polygon_for_check[i]){
        //     break;
        //   }
        //   if(i == boundary_for_check.size()-1){
        //     cout << "unprocessed_set: " << unprocessed_set.size() << endl;
        //     cout << "boundary_for_check size: " << boundary_for_check.size() << endl;
        //     cout << "boundary_for_check: "; disp_v(boundary_for_check);
        //     cout << "polygon size: " << polygon.size() << endl;
        //     cout << "polygon: "; disp_v(polygon);
        //     return 0; //Division Guaranteed
        //   }
        // }
        std::sort(boundary_for_check.begin(), boundary_for_check.end());
        std::sort(polygon_for_check.begin(), polygon_for_check.end());
        if(boundary_for_check == polygon_for_check){
          return 0; //Division Guaranteed
        }else{
          return -1; //Division Error
        }
      }else{
        cout << "unprocessed_set: " << unprocessed_set.size() << endl;
        cout << "boundary_for_check size: " << boundary_for_check.size() << endl;
        cout << "boundary_for_check: "; disp_v(boundary_for_check);
        cout << "polygon size: " << polygon.size() << endl;
        cout << "polygon: "; disp_v(polygon);
        return -1; //Division Error
      }
    }
    return 1;
  }
  int search_next_triangle(int shared_vert_index_a, int shared_vert_index_b){
    cnt_searchNextTriangle++;
    std::string key;
    std::unordered_set<int> triangles;
    if(shared_vert_index_a < shared_vert_index_b){
      key = std::to_string(shared_vert_index_a) + "-" + std::to_string(shared_vert_index_b);
    }else{
      key = std::to_string(shared_vert_index_b) + "-" + std::to_string(shared_vert_index_a);
    }
    triangles = edge_map[key];
    for(const auto& triangle_index : triangles){
      if(unprocessed_set.count(triangle_index)){
        return triangle_index;
      }
    }
    return -1; //don't have next_triangle
  }
  int get_another_vertex_index_from_triangle(int vert_index_a, int vert_index_b, int triangle_index){
    //三角形１つと頂点２つ与えられたときにもう１つの頂点を返す関数
    vector<int> vector_vert_index_c = tri[triangle_index];
    vector_vert_index_c.erase(std::remove(vector_vert_index_c.begin(),vector_vert_index_c.end(),vert_index_a), vector_vert_index_c.end());
    vector_vert_index_c.erase(std::remove(vector_vert_index_c.begin(),vector_vert_index_c.end(),vert_index_b), vector_vert_index_c.end());
    int vert_index_c = vector_vert_index_c[0];
    return vert_index_c;
  }
  /*
  int get_another_vertex_index_from_triangle(int vert_index_a, int vert_index_b, int triangle_index){
    //三角形１つと頂点２つ与えられたときにもう１つの頂点を返す関数
    vector<int> triangle = tri[triangle_index];
    int vert_index_c = -1;
    for(size_t i=0;i<3;++i){
      vert_index_c = triangle[i];
      if(vert_index_c != vert_index_a){
        if(vert_index_c != vert_index_b){
          return vert_index_c;
        }
      }
    }
    return vert_index_c;
  }
  */
  /*
  int get_another_vertex_index_from_triangle(int vert_index_a, int vert_index_b, int triangle_index){
    //三角形１つと頂点２つ与えられたときにもう１つの頂点を返す関数
    vector<int> triangle = tri[triangle_index];
    std::unordered_set<int> set(triangle.begin(),triangle.end());
    set.erase(vert_index_a);
    set.erase(vert_index_b);
    int vert_index_c = -1;
    for(const auto& index : set){
      vert_index_c = index;
    }
    return vert_index_c;
  }
  */

  // vector<pair<vector<double>, vector<double>>> get_orientation_candidates(vector<pair<vector<double>, vector<double>>> overlappingIntervals, vector<double> p1, vector<double> p2){
  //   double x_min = p1[0];
  //   double y_min = p1[1];
  //   double x_max = p2[0];
  //   double y_max = p2[1];
  //   vector<pair<vector<double>, vector<double>>> candidates;
  //   if(x_min > x_max){
  //     double x_tmp = x_min;
  //     x_min = x_max;
  //     x_max = x_tmp;
  //   }
  //   if(y_min > y_max){
  //     double y_tmp = y_min;
  //     y_min = y_max;
  //     y_max = y_tmp;
  //   }
  //   double x1,y1,x2,y2;
  //   vector<double> vertex1, vertex2;
  //   for(const auto& segment : overlappingIntervals){
  //     vertex1 = segment.first;
  //     vertex2 = segment.second;
  //     x1 = vertex1[0];
  //     y1 = vertex1[1];
  //     x2 = vertex2[0];
  //     y2 = vertex2[1];
  //     if(x1 <= x_min && x2 <= x_min) continue;
  //     if(x_max <= x1 && x_max <= x2) continue;
  //     if(y1 <= y_min && y2 <= y_min) continue;
  //     if(y_max <= y1 && y_max <= y2) continue;
  //     candidates.push_back(segment);
  //   }
  //   return candidates;
  // }

  vector<int> get_orientation_candidates(vector<int> polygon,vector<double> p1,vector<double> p2, int index_a){
    double x_min = p1[0];
    double y_min = p1[1];
    double x_max = p2[0];
    double y_max = p2[1];
    vector<int> candidates;
    if(x_min > x_max){
      double x_tmp = x_min;
      x_min = x_max;
      x_max = x_tmp;
    }
    if(y_min > y_max){
      double y_tmp = y_min;
      y_min = y_max;
      y_max = y_tmp;
    }
    //(index_a-1,index_a),(index_a,index_b),(index_b,index_b+1)の3辺は確実に必要ないので飛ばす
    double x1,y1,x2,y2;
    int index;
    int PS = polygon.size();
    for(int i=0;i<PS-3;++i){
      index = (index_a + 2 + i) % PS;
      x1 = vert[polygon[index]][0];
      y1 = vert[polygon[index]][1];
      x2 = vert[polygon[(index+1)%PS]][0];
      y2 = vert[polygon[(index+1)%PS]][1];
      if(x1 <= x_min && x2 <= x_min) continue;
      if(x_max <= x1 && x_max <= x2) continue;
      if(y1 <= y_min && y2 <= y_min) continue;
      if(y_max <= y1 && y_max <= y2) continue;
      candidates.push_back(index);
    }
    // disp_v(candidates);
    return candidates;
  }

  // vector<pair<vector<double>, vector<double>>> get_orientation_candidates(vector<double> p1, vector<double> p2){
  //   vector<pair<vector<double>, vector<double>>> overlappingIntervalsX = intervalTreeX.searchOverlapX({p1[0]},{p2[0]});
  //   vector<pair<vector<double>, vector<double>>> overlappingIntervalsY = intervalTreeY.searchOverlapY({p1[1]},{p2[1]});
  //   vector<pair<vector<double>, vector<double>>> candidates;
  //   for (const auto& element1 : overlappingIntervalsX) {
  //       for (const auto& element2 : overlappingIntervalsY) {
  //           if (element1 == element2) {
  //               candidates.push_back(element1);
  //               break; // 重複している要素が見つかったので、次の要素を調べる
  //           }
  //       }
  //   }
  //   return candidates;
  // }

  vector<vector<vector<double>>> get_orientation_candidates2(vector<double> p1, vector<double> p2){
    vector<vector<vector<double>>> overlappingIntervalsX;
    vector<vector<vector<double>>> overlappingIntervalsY;
    vector<vector<vector<double>>> candidates;
    // cout << "=====start=====" << endl;
    // disp_v(p1);
    // disp_v(p2);
    cnt_tree_search++;
    intervalTreeX.overlap_find_all(make_safe_interval<double>(p1[0],p2[0]), [&overlappingIntervalsX](auto iter){
      overlappingIntervalsX.push_back(iter.segment());
      return true;
    });
    intervalTreeY.overlap_find_all(make_safe_interval<double>(p1[1],p2[1]), [&overlappingIntervalsY](auto iter){
      overlappingIntervalsY.push_back(iter.segment());
      return true;
    });
    // for(int i=0;i<overlappingIntervalsX.size();i++){
    //   // disp_vv(overlappingIntervalsX[i]);
    //   for(int j=0;j<overlappingIntervalsY.size();j++){
    //     if(overlappingIntervalsX[i] == overlappingIntervalsY[j]){
    //       candidates.push_back(overlappingIntervalsX[i]);
    //       break;
    //     }
    //   }
    // }
    // cout << "overlappingIntervalsX size: " << overlappingIntervalsX.size() << endl;
    // cout << candidates.size() << endl;
    return candidates;
  }

  bool vv_equal(vector<double> p1, vector<double> p2, vector<double> q1, vector<double> q2){
    if(p1 == q1){
      if(p2 == q2){
        return true;
      }
    }else if(p1 == q2){
      if(p2 == q1){
        return true;
      }
    }
    return false;
  }

  int guarantee_division_for_non_convex(int& pre_pattern, int vertex_index_a, int vertex_index_b, int vertex_index_c, int triangle_index , vector<int>& polygon){
    auto ret_a = std::find(polygon.begin(), polygon.end(), vertex_index_a);
    int index_a_of_polygon = std::distance(polygon.begin(), ret_a);
    auto ret_b = std::find(polygon.begin(), polygon.end(), vertex_index_b);
    int index_b_of_polygon = std::distance(polygon.begin(), ret_b);
    auto ret_c = std::find(polygon.begin(), polygon.end(), vertex_index_c);
    int index_c_of_polygon = std::distance(polygon.begin(), ret_c);
    vector<int> candidates_pol;
    vector<vector<vector<double>>> candidates;
    vector<vector<double>> vv;
    // pair<vector<double>, vector<double>> segment;
    // vector<pair<vector<double>, vector<double>>> overlappingIntervals;
    // cout << intervalTreeX.size() << endl;
    double edge_pow2;
    if(tree_flag == 0){
      if(hyper_polygon_size < polygon.size()){
        tree_flag = 1;
        //intervalTree初期化
        int polygon_size = polygon.size();
        for(int i=0;i<polygon.size();i++){
          cnt_tree_insert += 2;
          intervalTreeX.insert(make_safe_interval<double>(vert[polygon[i % polygon_size]][0], vert[polygon[(i+1) % polygon_size]][0]), {vert[polygon[i % polygon_size]] , vert[polygon[(i+1) % polygon_size]]});
          intervalTreeY.insert(make_safe_interval<double>(vert[polygon[i % polygon_size]][1], vert[polygon[(i+1) % polygon_size]][1]), {vert[polygon[i % polygon_size]] , vert[polygon[(i+1) % polygon_size]]});
        }
      }
    }else if(tree_flag == 1){
      if(polygon.size() < hyper_end_polygon_size){
        tree_flag = 0;
        cout << "tree flag end! : cnt_checkConvexHull=" << cnt_checkConvexHull << endl;
        intervalTreeX.clear();
        intervalTreeY.clear();
      }
    }
    if(ret_c != polygon.end()){
      if( (index_b_of_polygon + 1) % polygon.size() == index_c_of_polygon % polygon.size()){
        if(verifyorientation(vert[vertex_index_a], vert[vertex_index_b], vert[vertex_index_c],1) > 0){
          if(tree_flag == 1){
            int N = 0;
            // CNT++;
            intervalTreeX.find_all(make_safe_interval<double>(vert[vertex_index_a][0], vert[vertex_index_b][0]), [&N, &vv, &vertex_index_a, &vertex_index_b](auto const& iter){
              // cout << "=====start=====" << endl;
              // vv = {vert[vertex_index_a],vert[vertex_index_b]};
              // disp_vv(vv);
              // disp_v(iter.segment()[0]);
              // disp_v(iter.segment()[1]);
              // if(vv_equal(vert[vertex_index_a], vert[vertex_index_b], iter.segment()[0], iter.segment()[1])){
                N++;
                intervalTreeX.erase(iter);
                cnt_tree_erase++;
                return false;
              // }
              // disp_vv(vv);
              // disp_v(iter.segment()[0]);
              // disp_v(iter.segment()[1]);
              // cout << "======end======" << endl;
              return true;
            });
            if(N == 0){cout << "!!!: " << CNT << endl;CNT++;}
            intervalTreeY.find_all(make_safe_interval<double>(vert[vertex_index_a][1], vert[vertex_index_b][1]), [&vv, &vertex_index_a, &vertex_index_b](auto const& iter){
              // vv = {vert[vertex_index_a],vert[vertex_index_b]};
              // if(vv_equal(vert[vertex_index_a], vert[vertex_index_b], iter.segment()[0], iter.segment()[1])){
                intervalTreeY.erase(iter);
                cnt_tree_erase++;
                return false;
              // }
              return true;
            });

            candidates = get_orientation_candidates2(vert[vertex_index_a], vert[vertex_index_c]);
            for(const auto& candidate : candidates){
              if(intersection(vert[vertex_index_a], vert[vertex_index_c], candidate[0], candidate[1]) == -1){
                return true;
              }
            }

            intervalTreeX.find_all(make_safe_interval<double>(vert[vertex_index_b][0], vert[vertex_index_c][0]), [&vv, &vertex_index_b, &vertex_index_c](auto const& iter){
              // vv = {vert[vertex_index_b],vert[vertex_index_c]};
              // if(vv_equal(vert[vertex_index_b], vert[vertex_index_c], iter.segment()[0], iter.segment()[1])){
                intervalTreeX.erase(iter);
                cnt_tree_erase++;
                return false;
              // }
              return true;
            });
            intervalTreeX.insert(make_safe_interval<double>(vert[vertex_index_a][0], vert[vertex_index_c][0]), {{vert[vertex_index_a][0],vert[vertex_index_a][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
            cnt_tree_insert++;
            intervalTreeY.find_all(make_safe_interval<double>(vert[vertex_index_b][1], vert[vertex_index_c][1]), [&vv, &vertex_index_b, &vertex_index_c](auto const& iter){
              // vv = {vert[vertex_index_b],vert[vertex_index_c]};
              // if(vv_equal(vert[vertex_index_b], vert[vertex_index_c], iter.segment()[0], iter.segment()[1])){
                intervalTreeY.erase(iter);
                cnt_tree_erase;
                return false;
              // }
              return true;
            });
            intervalTreeY.insert(make_safe_interval<double>(vert[vertex_index_a][1], vert[vertex_index_c][1]), {{vert[vertex_index_a][0],vert[vertex_index_a][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
            cnt_tree_insert++;
          }else{
            candidates_pol = get_orientation_candidates(polygon,vert[vertex_index_a], vert[vertex_index_c],index_a_of_polygon);
            for(int i=0; i<candidates_pol.size(); ++i){
              if(i != index_a_of_polygon){
                if(intersection(vert[vertex_index_a], vert[vertex_index_c], vert[polygon[candidates_pol[i]]], vert[polygon[(candidates_pol[i]+1)%polygon.size()]]) == -1){
                  return true;
                }
              }
            }
          }

          polygon.erase(ret_b);
          unprocessed_set.erase(triangle_index);
          pre_pattern = 1; //1個進む
          return false; //share 2 sides and 3 points
        }else{
          return true;
        }
      }else{
        if(verifyorientation(vert[vertex_index_a], vert[vertex_index_b], vert[vertex_index_c],2) > 0){
          pre_pattern = 4;
          return false; //skip : share 1 sides and 3 points または 2辺3点共有だけどc,a,bの順になっている場合
        }
      }
    }else{
      // cout << "a:" << vertex_index_a << ",b:" << vertex_index_b << ",c:" << vertex_index_c << endl;
      if(verifyorientation(vert[vertex_index_a], vert[vertex_index_b], vert[vertex_index_c],3) > 0){
        if(verifyorientation(vert[vertex_index_a], vert[vertex_index_b], vert[polygon[(index_a_of_polygon+polygon.size()-1)%polygon.size()]],4) > 0){
          if(verifyorientation(vert[polygon[(index_a_of_polygon+polygon.size()-1)%polygon.size()]], vert[vertex_index_a], vert[vertex_index_c],5) < 0){
            return true;
          }
        }
        if(verifyorientation(vert[vertex_index_a], vert[vertex_index_b], vert[polygon[(index_b_of_polygon+1)%polygon.size()]],6) > 0){
          if(verifyorientation(vert[vertex_index_b], vert[polygon[(index_b_of_polygon+1)%polygon.size()]], vert[vertex_index_c],7) < 0){
            return true;
          }
        }
        if(tree_flag == 1){
          intervalTreeX.find_all(make_safe_interval<double>(vert[vertex_index_a][0], vert[vertex_index_b][0]), [&vv, &vertex_index_a, &vertex_index_b](auto const& iter){
            // vv = {vert[vertex_index_a],vert[vertex_index_b]};
            // if(vv_equal(vert[vertex_index_a], vert[vertex_index_b], iter.segment()[0], iter.segment()[1])){
              intervalTreeX.erase(iter);
              cnt_tree_erase++;
              return false;
            // }
            return true;
          });
          intervalTreeY.find_all(make_safe_interval<double>(vert[vertex_index_a][1], vert[vertex_index_b][1]), [&vv, &vertex_index_a, &vertex_index_b](auto const& iter){
            // vv = {vert[vertex_index_a],vert[vertex_index_b]};
            // if(vv_equal(vert[vertex_index_a], vert[vertex_index_b], iter.segment()[0], iter.segment()[1])){
              intervalTreeY.erase(iter);
              cnt_tree_erase++;
              return false;
            // }
            return true;
          });
          candidates = get_orientation_candidates2(vert[vertex_index_a], vert[vertex_index_c]);
          for(const auto& candidate : candidates){
            if(intersection(vert[vertex_index_a], vert[vertex_index_c], candidate[0], candidate[1]) == -1){
              return true;
            }
          }
          candidates = get_orientation_candidates2(vert[vertex_index_b], vert[vertex_index_c]);
          for(const auto& candidate : candidates){
            if(intersection(vert[vertex_index_b], vert[vertex_index_c], candidate[0], candidate[1]) == -1){
              return true;
            }
          }
          intervalTreeX.insert(make_safe_interval<double>(vert[vertex_index_a][0], vert[vertex_index_c][0]), {{vert[vertex_index_a][0],vert[vertex_index_a][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
          intervalTreeX.insert(make_safe_interval<double>(vert[vertex_index_b][0], vert[vertex_index_c][0]), {{vert[vertex_index_b][0],vert[vertex_index_b][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
          intervalTreeY.insert(make_safe_interval<double>(vert[vertex_index_a][1], vert[vertex_index_c][1]), {{vert[vertex_index_a][0],vert[vertex_index_a][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
          intervalTreeY.insert(make_safe_interval<double>(vert[vertex_index_b][1], vert[vertex_index_c][1]), {{vert[vertex_index_b][0],vert[vertex_index_b][1]},{vert[vertex_index_c][0],vert[vertex_index_c][1]}});
          cnt_tree_insert += 4;
        }else{
          candidates_pol = get_orientation_candidates(polygon,vert[vertex_index_a], vert[vertex_index_c],index_a_of_polygon);
          for(int i=0; i<candidates_pol.size(); ++i){
            if(intersection(vert[vertex_index_a], vert[vertex_index_c], vert[polygon[candidates_pol[i]]], vert[polygon[(candidates_pol[i]+1)%polygon.size()]]) == -1){
              return true;
            }
          }
          candidates_pol = get_orientation_candidates(polygon,vert[vertex_index_b], vert[vertex_index_c],index_a_of_polygon);
          for(int i=0; i<candidates_pol.size(); ++i){
            if(intersection(vert[vertex_index_b], vert[vertex_index_c], vert[polygon[candidates_pol[i]]], vert[polygon[(candidates_pol[i]+1)%polygon.size()]]) == -1){
              return true;
            }
          }
        }

        vector<int>::iterator itr;
        itr = std::find(polygon.begin(), polygon.end(), vertex_index_a);
        polygon.insert(itr+1, vertex_index_c);
        unprocessed_set.erase(triangle_index);
        pre_pattern = 2;
        return false; //share 1 sides and 2 points
      }else{
        return true;
      }
    }
    return 0;
  }
  int guarantee_division(int init_tri_index){
    int index = 0;
    int flag = 0;
    int flag_hole = 0;
    int pre_pattern = 0;// 0(初期値) or 1(非凸2辺3点共有のとき) or 2(凸または非凸1辺2点共有のとき) or 4(非凸1辺3点共有のとき(mod3で1))
    unprocessed_set.erase(init_tri_index);
    vector<int> polygon = init_polygon(init_tri_index);
    while(1){
      index = index % polygon.size();
      // cout << "polygon size: " << polygon.size() << endl;
      // cout << "index: " << index << endl;
      // cout << "pre_pattern: " << pre_pattern << endl;
      int cch = check_convex_hull(polygon);
      if(cch == 0){
        return 0;
      }else if(cch == -1){
        return -1;
      }
      int shared_vert_index_a = polygon[index % polygon.size()];
      int shared_vert_index_b = polygon[(index+1) % polygon.size()];
      int next_triangle_index = search_next_triangle(shared_vert_index_a, shared_vert_index_b);
      // cout << "nti: " << next_triangle_index << endl;
      if(next_triangle_index == -1){
        index++;
        flag++;
        continue;
      }else{
        if(pre_pattern != 1){
          // cout << "!!!!!1" << endl;
          index++; //2辺3点共有以外は+1する
        }
        int vert_index_c = get_another_vertex_index_from_triangle(shared_vert_index_a, shared_vert_index_b, next_triangle_index);
        // cout << "vert_index_c: " << vert_index_c << endl;
        if(guarantee_division_for_non_convex(pre_pattern, shared_vert_index_a, shared_vert_index_b, vert_index_c, next_triangle_index, polygon)){
          return -1;
        }
        if(pre_pattern != 4){
          // cout << "!!!!!2" << endl;
          flag = 0;
        }
      }
      // return -1;
    }
    return 0;
  }
#pragma endregion Divison
#pragma region Delaunay
void get_incircle_candidates(int& lower, int& upper,double x_min, double x_max){
  //x軸に関して途中まで二分探索風にしてあとは細かく探索
  int left = 0;
  int right = x_sorted_vert.size()-1;
  int mid;
  while(right - left > 1){
    mid = (left+right)/2;
    if(x_sorted_vert[mid][0] <= x_min){
      left = mid;
    }
    else if(x_max <= x_sorted_vert[mid][0]){
      right = mid;
    }
    else{
      break;
    }
  }
  lower = mid;
  upper = mid;
  while( lower != 0 && x_min < x_sorted_vert[lower-1][0]){
    lower--;
  }
  while( upper != x_sorted_vert.size()-1 && x_sorted_vert[upper+1][0] < x_max){
    upper++;
  }
}

//違反辺かどうかをverifyincircleで確認して、違反ならFlipする
void edge_flip(const int tri_index1, const int tri_index2){
  // cout << "triangle1: "; disp_v(tri[tri_index1]);
  // cout << "triangle2: "; disp_v(tri[tri_index2]);
  vector<int> triangle1 = tri[tri_index1];
  vector<int> triangle2 = tri[tri_index2];
  vector<int> all = triangle1;
  all.insert(all.end(), triangle2.begin(), triangle2.end());
  std::sort(all.begin(), all.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());
  vector<int> non_shared1 = all;
  vector<int> non_shared2 = all;
  vector<int> shared = all;
  for(size_t i=0; i<3; ++i){
    non_shared1.erase(std::remove(non_shared1.begin(), non_shared1.end(), triangle2[i]), non_shared1.end());
  }
  for(size_t i=0; i<3; ++i){
    non_shared2.erase(std::remove(non_shared2.begin(), non_shared2.end(), triangle1[i]), non_shared2.end());
  }
  shared.erase(std::remove(shared.begin(), shared.end(), non_shared1[0]), shared.end());
  shared.erase(std::remove(shared.begin(), shared.end(), non_shared2[0]), shared.end());

  if(shared[0] < shared[1]){
    edges.erase(std::to_string(shared[0])+'-'+std::to_string(shared[1]));
  }else{
    edges.erase(std::to_string(shared[1])+'-'+std::to_string(shared[0]));
  }
  if(verifyincircle(vert[triangle1[0]], vert[triangle1[1]], vert[triangle1[2]], vert[non_shared2[0]]) < 0){
    // cout << "incircle out!!!" << endl;
    //違反辺
    cnt_nonDelaunay++;
    //Flip
    vector<int> non_shared;
    non_shared.push_back(non_shared1[0]);
    non_shared.push_back(non_shared2[0]);
    vector<int> new_triangle1 = non_shared;
    vector<int> new_triangle2 = non_shared;
    new_triangle1.push_back(shared[0]);
    new_triangle2.push_back(shared[1]);
    std::sort(new_triangle1.begin(), new_triangle1.end());
    std::sort(new_triangle2.begin(), new_triangle2.end());
    tri[tri_index1] = new_triangle1;
    tri[tri_index2] = new_triangle2;
    // cout << "triangle1: "; disp_v(tri[tri_index1]);
    // cout << "triangle2: "; disp_v(tri[tri_index2]);
    cout << "new_triangle1: "; disp_v(tri[tri_index1]);
    cout << "new_triangle2: "; disp_v(tri[tri_index2]);
    //edge_map更新
    std::string key_edge;
    if(shared[0] < shared[1]){
      key_edge = std::to_string(shared[0])+'-'+std::to_string(shared[1]);
    }else{
      key_edge = std::to_string(shared[1])+'-'+std::to_string(shared[0]);
    }
    edge_map.erase(key_edge);
    std::string new_key_edge;
    if(non_shared1[0] < non_shared2[0]){
      new_key_edge = std::to_string(non_shared1[0])+'-'+std::to_string(non_shared2[0]);
    }else{
      new_key_edge = std::to_string(non_shared2[0])+'-'+std::to_string(non_shared1[0]);
    }
    std::unordered_set<int> triangles_index;
    triangles_index.insert(tri_index1);
    triangles_index.insert(tri_index2);
    edge_map[new_key_edge] = triangles_index;
    
    //edgesに周りの辺追加
    std::string edge;
    triangles_index.clear();
    for(int i=0; i<2; ++i){
      if(non_shared1[0] < shared[i]){
        edge = std::to_string(non_shared1[0])+'-'+std::to_string(shared[i]);
        if(boundary_edges.find(edge) == boundary_edges.end()) edges.insert(edge);
        triangles_index = edge_map[edge];
        triangles_index.erase(tri_index1);
        triangles_index.erase(tri_index2);
        if(i==0){
          triangles_index.insert(tri_index1);
        }else{
          triangles_index.insert(tri_index2);
        }
        edge_map[edge] = triangles_index;
      }else{
        edge = std::to_string(shared[i])+'-'+std::to_string(non_shared1[0]);
        if(boundary_edges.find(edge) == boundary_edges.end()) edges.insert(edge);
        triangles_index = edge_map[edge];
        triangles_index.erase(tri_index1);
        triangles_index.erase(tri_index2);
        if(i==0){
          triangles_index.insert(tri_index1);
        }else{
          triangles_index.insert(tri_index2);
        }
        edge_map[edge] = triangles_index;
      }
      if(non_shared2[0] < shared[i]){
        edge = std::to_string(non_shared2[0])+'-'+std::to_string(shared[i]);
        if(boundary_edges.find(edge) == boundary_edges.end()) edges.insert(edge);
        triangles_index = edge_map[edge];
        triangles_index.erase(tri_index1);
        triangles_index.erase(tri_index2);
        if(i==0){
          triangles_index.insert(tri_index1);
        }else{
          triangles_index.insert(tri_index2);
        }
        edge_map[edge] = triangles_index;
      }else{
        edge = std::to_string(shared[i])+'-'+std::to_string(non_shared2[0]);
        if(boundary_edges.find(edge) == boundary_edges.end()) edges.insert(edge);
        triangles_index = edge_map[edge];
        triangles_index.erase(tri_index1);
        triangles_index.erase(tri_index2);
        if(i==0){
          triangles_index.insert(tri_index1);
        }else{
          triangles_index.insert(tri_index2);
        }
        edge_map[edge] = triangles_index;
      }
    }
  }
}

void guarantee_delaunay(){
  vector<int> triangle1;
  vector<int> triangle2;
  // cout << "!!!!!!" << endl;
  // disp_v(boundary);
  while(0 < edges.size()){
    // cout << edges.size() << endl;
    auto edge_itr = edges.begin();
    // cout << "!1" << endl;
    auto triangles = edge_map[*edge_itr];
    // cout << triangles.size() << endl;
    if(triangles.size()<2){
      // cout << *edge_itr << endl;
      edges.erase(*edge_itr);
    }else{
      auto triangle_itr = triangles.begin();
      int triangle1_index = *triangle_itr;
      triangle_itr++;
      // cout << "!2" << endl;
      // cout << edges.size() << endl;
      // cout << *triangle_itr << endl;
      int triangle2_index = *triangle_itr;
      // cout << "!3" << endl;
      // cout << "!!!!!!" << endl;
      edge_flip(triangle1_index,triangle2_index);
      // cout << "!4" << endl;
    }

    // auto triangle_itr = triangles.begin();
    // int triangle1_index = *triangle_itr;
    // triangle_itr++;
    // // cout << "!2" << endl;
    // // cout << edges.size() << endl;
    // // cout << *triangle_itr << endl;
    // int triangle2_index = *triangle_itr;
    // // cout << "!3" << endl;
    // // cout << "!!!!!!" << endl;
    // edge_flip(triangle1_index,triangle2_index);
    // // cout << "!4" << endl;
  }
}
#pragma endregion Delaunay
}

int main(){
  using pstv::vert; using pstv::tri; using pstv::boundary;
  using pstv::x_sorted_vert; using pstv::y_sorted_vert;
  using pstv::edge_map; using pstv::boundary_edges; using pstv::edges; using pstv::unprocessed_set;
  using pstv::vertex_map; using pstv::triangle_map;
  using pstv::intervalTreeX; using pstv::intervalTreeY;
  using pstv::hyper_polygon_size; using pstv::hyper_end_polygon_size; using pstv::tree_flag;

  std::string dir = "../data/1000/1/";
  vert = pstv::get_data<double>(dir+"vert.csv");
  pstv::get_xy_sort();
  tri = pstv::get_data<int>(dir+"tri.csv");
  vector<vector<int>> boundary_tmp = pstv::get_data<int>(dir+"boundary.csv");
  for(size_t i=0; i<boundary_tmp.size(); ++i){
    boundary.push_back(boundary_tmp[i][0]);
  }
  edge_map = pstv::get_edge_map();
  boundary_edges = pstv::get_boundary_edges();
  edges = pstv::get_edges();
  unprocessed_set = pstv::get_unprocessed_set();
  pstv::filter();
  size_t overlap = pstv::check_overlap_vert();
  if(overlap){
    if(pstv::guarantee_division(1) == 0){
      cout << "Division Guaranteed." << endl;
      vertex_map = pstv::get_vertex_map();
      triangle_map = pstv::get_triangle_map();
      pstv::guarantee_delaunay();
      if(0 < cnt_nonDelaunay){
        cout << "Not Fullfilling Delaunay. nonDelaunay Edge number: " << cnt_nonDelaunay << endl;
        pstv::output_tri_data(dir+"new_tri.csv");
      }else{
        cout << "Delaunay Fullfilled." << endl;
      }
    }else{
      cout << "Triangulation has Error!!!!." << endl;
    }
  }else{
    cout << "Triangulation Overlap!!!!" << endl;
  }

}