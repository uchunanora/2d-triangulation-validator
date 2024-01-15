#include <pstv/dataset.hpp>
#include <pstv/validation.hpp>
#include <pstv/validator_delaunay.hpp>
#include <pstv/validator_triangulation.hpp>
#include <utils/display.hpp>

int main() {
  std::string dir = "./dataset/nonVerifiedDelaunay/";
  pstv::Dataset dataset(dir + "vertexes.csv", dir + "triangles.csv", dir + "boundaries.csv");
  if (dataset.check_overlap_vertexes()) {
    pstv::ValidatorTriangulation vt;
    bool is_verified_triangulation = vt.validate(dataset);
    if (is_verified_triangulation) {
      pstv::ValidatorDelaunay vd;
      bool is_verified_delaunay = vd.validate(dataset, true);
      if (!is_verified_delaunay) {
        std::string outpit_dir = "./new_tri.csv";
        vd.output_new_triangles_data(outpit_dir);
      }
    }
  }
}