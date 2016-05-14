#ifndef _flowmap_
#define _flowmap_

#include "ftle.h"
#include "parameter.h"
#include "matrix.h"
#include "rk45.h"

using namespace std;

class flowmap {
 public:
  flowmap(const shared_ptr<parameter> &objpara);
  void calculate_trajectories();
  virtual ~flowmap();
  vector3d<double> phi_x() const { return phi_x_; };
  vector3d<double> phi_y() const { return phi_y_; };
  vector3d<double> phi_z() const { return phi_z_; };
  vector3d<double> phi_dudx() const { return phi_dudx_; };
  vector3d<double> phi_dudy() const { return phi_dudy_; };
  vector3d<double> phi_dudz() const { return phi_dudz_; };
  vector3d<double> phi_dvdx() const { return phi_dvdx_; };
  vector3d<double> phi_dvdy() const { return phi_dvdy_; };
  vector3d<double> phi_dvdz() const { return phi_dvdz_; };
  vector3d<double> phi_dwdx() const { return phi_dwdx_; };
  vector3d<double> phi_dwdy() const { return phi_dwdy_; };
  vector3d<double> phi_dwdz() const { return phi_dwdz_; };
  vector<double> x() const { return x_; };
  vector<double> y() const { return y_; };
  vector<double> z() const { return z_; };

 private:
  int nx;
  int ny;
  int nz;
  void initialize_expression(int number_of_threads);
  vector<double> velocity(double t, vector<double> &y);
  void rk45(double t0, double t1, double tol, double hmin, double hmax, size_t maxiter, vector<double> &d);

  exprtk::expression<double> *exp_u;
  exprtk::expression<double> *exp_v;
  exprtk::expression<double> *exp_w;
  exprtk::expression<double> *exp_dudx;
  exprtk::expression<double> *exp_dudy;
  exprtk::expression<double> *exp_dudz;
  exprtk::expression<double> *exp_dvdx;
  exprtk::expression<double> *exp_dvdy;
  exprtk::expression<double> *exp_dvdz;
  exprtk::expression<double> *exp_dwdx;
  exprtk::expression<double> *exp_dwdy;
  exprtk::expression<double> *exp_dwdz;
  vector<double> var_x;
  vector<double> var_y;
  vector<double> var_z;
  vector<double> var_t;
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  double zmin;
  double zmax;
  string function_u;
  string function_v;
  string function_w;
  string function_dudx;
  string function_dudy;
  string function_dudz;
  string function_dvdx;
  string function_dvdy;
  string function_dvdz;
  string function_dwdx;
  string function_dwdy;
  string function_dwdz;
  double t;
  double t0;
  double tol;
  double hmin;
  double hmax;
  vector3d<double> phi_x_;
  vector3d<double> phi_y_;
  vector3d<double> phi_z_;
  vector3d<double> phi_dudx_;
  vector3d<double> phi_dudy_;
  vector3d<double> phi_dudz_;
  vector3d<double> phi_dvdx_;
  vector3d<double> phi_dvdy_;
  vector3d<double> phi_dvdz_;
  vector3d<double> phi_dwdx_;
  vector3d<double> phi_dwdy_;
  vector3d<double> phi_dwdz_;
  vector<double> x_;
  vector<double> y_;
  vector<double> z_;
};

#endif
