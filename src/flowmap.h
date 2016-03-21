#ifndef _flowmap_
#define _flowmap_

#include "ftle.h"
#include "parameter.h"
#include "matrix.h"

using namespace std;

class flowmap {
 public:
  flowmap(shared_ptr<parameter> &objpara);
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
  void ode45(double t0, double t1, double tol, double hmin, double hmax, size_t maxiter, vector<double> &d);

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

  // ode45 coefficients
  const double a21 = 1.0/5.0;
  const double a31 = 3.0/40.0;
  const double a32 = 9.0/40.0;
  const double a41 = 44.0/45.0;
  const double a42 =-56.0/15.0;
  const double a43 = 32.0/9.0;
  const double a51 = 19372.0/6561.0;
  const double a52 =-25360.0/2187.0;
  const double a53 = 64448.0/6561.0;
  const double a54 =-212.0/729.0;
  const double a61 = 9017.0/3168.0;
  const double a62 =-355.0/33.0;
  const double a63 = 46732.0/5247.0;
  const double a64 = 49.0/176.0;
  const double a65 =-5103.0/18656.0;
  const double a71 = 35.0/384.0;
  //const double a72 = 0.0
  const double a73 = 500.0/1113.0;
  const double a74 = 125.0/192.0;
  const double a75 =-2187.0/6784.0;
  const double a76 = 11.0/84.0;

  //const double c1 = 0.0
  const double c2 = 1.0/5.0;
  const double c3 = 3.0/10.0;
  const double c4 = 4.0/5.0;
  const double c5 = 8.0/9.0;
  //const double c6 = 1.0;
  //const double c7 = 1.0;

  const double b1 = 35.0/384.0;
  //const double b2 = 0.0;
  const double b3 = 500.0/1113.0;
  const double b4 = 125.0/192.0;
  const double b5 =-2187.0/6784.0;
  const double b6 = 11.0/84.0;
  //const double b7 = 0.0

  const double b1p = 5179.0/57600.0;
  //const double b2p = 0.0
  const double b3p = 7571.0/16695.0;
  const double b4p = 393.0/640.0;
  const double b5p =-92097.0/339200.0;
  const double b6p = 187.0/2100.0;
  const double b7p = 1.0/40.0;
};

#endif
