#include <iostream>
#include <pstv/dataset.hpp>
#include <pstv/validation.hpp>
#include <pstv/validator_triangulation.hpp>
#include <utils/display.hpp>

int main() {
  std::string dir = "../data/1000/1/";
  pstv::Dataset dataset(dir + "vert.csv", dir + "tri.csv", dir + "boundary.csv");
  if (dataset.check_overlap_vertexes()) {
    pstv::ValidatorTriangulation vt;
    bool v = vt.validate(dataset);
  }
}