#ifndef DATASET_HPP
#define DATASET_HPP

#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

using std::vector;

namespace pstv {
  class Dataset {
  public:
    vector<vector<double>> vertexes;
    vector<vector<int>> triangles;
    vector<int> boundaries;
    std::unordered_map<std::string,std::unordered_set<int>> edge_map;
    std::unordered_set<std::string> boundary_edges;
    std::unordered_set<std::string> edges;
    std::unordered_set<int> unprocessed_set;
  
    Dataset(std::string vertex_file, std::string triangle_file, std::string boundary_file){
      set_vertex_from_csv(vertex_file);
      set_triangle_from_csv(triangle_file);
      set_boundary_from_csv(boundary_file);
      set_edge_map();
      set_boundary_edges();
      set_edges();
      set_unprocessed_set();
    };
    
  private:
    template<class _T> vector<vector<_T>> get_data(std::string filename){
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
          }
          vv[vv.size()-1].push_back(num);
        }
      }
      return vv;
    } 
    void set_vertex_from_csv(std::string filename){
      vertexes = get_data<double>(filename);
    }
    void set_triangle_from_csv(std::string filename){
      triangles = get_data<int>(filename);
    }
    void set_boundary_from_csv(std::string filename){
      vector<vector<int>> boundary_tmp = get_data<int>(filename);
      for(size_t i=0; i<boundary_tmp.size(); ++i){
        boundaries.push_back(boundary_tmp[i][0]);
      }
    }

    void set_edge_map(){
      vector<int> triangle;
      std::unordered_set<int> tmp;
      std::string key;
      // #pragma omp parallel for
      for(int i=0,n=triangles.size();i<n;++i){
        triangle = triangles[i];
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
    }
    void set_boundary_edges(){
      for(int i=0,n=boundaries.size();i<n;++i){
        if(boundaries[i%n] < boundaries[(i+1)%n]){
          boundary_edges.insert(std::to_string(boundaries[i%n])+'-'+std::to_string(boundaries[(i+1)%n]));
        }else{
          boundary_edges.insert(std::to_string(boundaries[(i+1)%n])+'-'+std::to_string(boundaries[i%n]));
        }
      }
    }
    void set_edges(){
      vector<int> triangle;
      std::string edge;
      // #pragma omp parallel for
      for(int i=0,n=triangles.size();i<n;++i){
        triangle = triangles[i];
        std::sort(triangle.begin(),triangle.end());
        edge = std::to_string(triangle[0])+'-'+std::to_string(triangle[1]);
        edges.insert(edge);
        edge = std::to_string(triangle[1])+'-'+std::to_string(triangle[2]);
        edges.insert(edge);
        edge = std::to_string(triangle[0])+'-'+std::to_string(triangle[2]);
        edges.insert(edge);
      }
      for(int i=0,n=boundaries.size();i<n;++i){
        if(boundaries[i%n] < boundaries[(i+1)%n]){
          edges.erase(std::to_string(boundaries[i%n])+'-'+std::to_string(boundaries[(i+1)%n]));
        }else{
          edges.erase(std::to_string(boundaries[(i+1)%n])+'-'+std::to_string(boundaries[i%n]));
        }
      }
    }
    void set_unprocessed_set(){
      for(int i=0,n=triangles.size();i<n;++i){
        unprocessed_set.insert(i);
      }
    }
  };
}
#endif // DATASET_HPP