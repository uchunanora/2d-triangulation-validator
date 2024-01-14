#include <iostream>
#include <pstv/dataset.hpp>
#include <pstv/varidation.hpp>
#include <pstv/validator_triangulation.hpp>
#include <utils/display.hpp>

int main()
{
  std::string dir = "../data/1000/1/";
  pstv::Dataset dataset(dir + "vert.csv", dir + "tri.csv", dir + "boundary.csv");
  pstv::Validation validation;
  double a = validation.orientation(dataset.vertexes[0], dataset.vertexes[1], dataset.vertexes[2]);
  std::cout << a << std::endl;
  bool b = dataset.check_overlap_vertexes();
  std::cout << b << std::endl;
  pstv::ValidatorTriangulation VT;
  bool v = VT.validate(dataset);
}