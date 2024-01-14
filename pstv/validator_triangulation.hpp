#ifndef VALIDATOR_TRIANGULATION_HPP
#define VALIDATOR_TRIANGULATION_HPP

#include <iostream>
#include <vector>
#include <pstv/dataset.hpp>
#include <pstv/varidation.hpp>

namespace pstv
{
  class ValidatorTriangulation
  {
    pstv::Validation varidation;
    pstv::Dataset dataset;

  public:
    explicit ValidatorTriangulation(){};

    bool validate(pstv::Dataset ds, int init_tri_index = 0)
    {
      int index = 0;
      int flag = 0;
      int pre_pattern = 0;
      dataset = ds;
      std::cout << dataset.unprocessed_set.size() << std::endl;
      dataset.unprocessed_set.erase(init_tri_index);
      std::cout << dataset.unprocessed_set.size() << std::endl;
      vector<int> polygon = init_polygon(init_tri_index);
      return 0;
    }

    vector<int> init_polygon(int init_tri_index)
    {
      update_tri_clockwise(dataset.triangles, init_tri_index);
      vector<int> polygon(3, 0);
      polygon.reserve(dataset.vertexes.size());
      for (size_t i = 0; i < 3; ++i)
      {
        polygon[i] = dataset.triangles[init_tri_index][i];
      }
      return polygon;
    }

    void update_tri_clockwise(vector<vector<int>> &tri, int tri_index, int mode = 0)
    {
      // mode=0:clockwise sort, mode=1:unclockwise sort
      vector<double> vert_A = vert[tri[tri_index][0]];
      vector<double> vert_B = vert[tri[tri_index][1]];
      vector<double> vert_C = vert[tri[tri_index][2]];
      vector<double> vector_AB = {vert_B[0] - vert_A[0], vert_B[1] - vert_A[1]};
      vector<double> vector_AC = {vert_C[0] - vert_A[0], vert_C[1] - vert_A[1]};
      double cross_product = vector_AB[0] * vector_AC[1] - vector_AB[1] * vector_AC[0];
      if (fabs(cross_product) <= theta * (fabs(vector_AB[0] * vector_AC[1] + vector_AB[1] * vector_AC[0])))
      {
        mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_cross_product;
        mpq_ax = vert[tri[tri_index][0]][0];
        mpq_ay = vert[tri[tri_index][0]][1];
        mpq_bx = vert[tri[tri_index][1]][0];
        mpq_by = vert[tri[tri_index][1]][1];
        mpq_cx = vert[tri[tri_index][2]][0];
        mpq_cy = vert[tri[tri_index][2]][1];
        mpq_cross_product = (mpq_bx - mpq_ax) * (mpq_cy - mpq_ay) - (mpq_by - mpq_ay) * (mpq_cx - mpq_ax);
        int mpq_sgn = sgn(mpq_cross_product);
        if (mode == 0)
        {
          if (mpq_sgn > 0)
          {
            std::reverse(tri[tri_index].begin(), tri[tri_index].end());
          }
        }
        else if (mode == 1)
        {
          if (mpq_sgn > 0)
          {
            std::reverse(tri[tri_index].begin(), tri[tri_index].end());
          }
        }
      }
      else
      {
        if (mode == 0)
        {
          if (cross_product > 0)
          {
            std::reverse(tri[tri_index].begin(), tri[tri_index].end());
          }
        }
        else if (mode == 1)
        {
          if (cross_product < 0)
          {
            std::reverse(tri[tri_index].begin(), tri[tri_index].end());
          }
        }
      }
    }
  };
}

#endif