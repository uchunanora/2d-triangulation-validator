#ifndef VARIDATION_HPP
#define VARIDATION_HPP

#include <vector>
#include <math.h>
#include <gmpxx.h>

using std::vector;

namespace pstv {
  class Validation {
    double u,theta,const6u,u_n,iccerrboundA;
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
  public:
    Validation(){
      filter();
    }
    double orientation(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c){
      const double axby = (vert_a[0] - vert_c[0]) * (vert_b[1] - vert_c[1]);
      const double aybx = (vert_a[1] - vert_c[1]) * (vert_b[0] - vert_c[0]);
      const double det = axby - aybx;
      if(fabs(det) <= theta*(fabs(axby + aybx) + u_n)){
        mpq_class mpq_ax, mpq_ay, mpq_bx, mpq_by, mpq_cx, mpq_cy, mpq_det;
        mpq_ax = vert_a[0];
        mpq_ay = vert_a[1];
        mpq_bx = vert_b[0];
        mpq_by = vert_b[1];
        mpq_cx = vert_c[0];
        mpq_cy = vert_c[1];
        mpq_det = (mpq_ax - mpq_cx)*(mpq_by - mpq_cy) - (mpq_ay - mpq_cy)*(mpq_bx - mpq_cx);
        int mpq_sgn = sgn(mpq_det); // -1 or 0 or 1
        return mpq_sgn;
      }
      return det;
    }

    double incircle(vector<double> vertex_a, vector<double> vertex_b, vector<double> vertex_c, vector<double> vertex_d){
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
      vector<double> vector_x = {vert_a[0],vert_b[0],vert_c[0]};
      vector<double> vector_y = {vert_a[1],vert_b[1],vert_c[1]};
      std::sort(vector_x.begin(),vector_x.end());
      std::sort(vector_y.begin(),vector_y.end());
      if(orientation(vert_a, vert_b, vert_c) == 0){
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

    int onsegment(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c){
      if(std::min(vert_a[0],vert_b[0]) < vert_c[0] && vert_c[0] < std::max(vert_a[0],vert_b[0]) && std::min(vert_a[1],vert_b[1]) < vert_c[1] && vert_c[1] < std::max(vert_a[1],vert_b[1])){
        return true;
      }
      return false;
    }

    int intersection(vector<double> vert_a, vector<double> vert_b, vector<double> vert_c, vector<double> vert_d){
      double o1,o2,o3,o4;
      o1 = orientation(vert_a,vert_b,vert_c);
      o2 = orientation(vert_a,vert_b,vert_d);
      o3 = orientation(vert_c,vert_d,vert_a);
      o4 = orientation(vert_c,vert_d,vert_b);
      if(o1*o2 < 0 && o3*o4 < 0) return -1;
      if(o1==0 && onsegment(vert_a,vert_b,vert_c)) return -1;
      if(o2==0 && onsegment(vert_a,vert_b,vert_d)) return -1;
      if(o3==0 && onsegment(vert_c,vert_d,vert_a)) return -1;
      if(o4==0 && onsegment(vert_c,vert_d,vert_b)) return -1;
      return 1;
    }

  };
}

#endif //VALIDATION_HPP