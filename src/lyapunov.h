#ifndef _lyapunov_
#define _lyapunov_

#include "ftle.h"
#include "flowmap.h"
#include "parameter.h"

using namespace std;

class lyapunov {
 public:
  lyapunov(shared_ptr<parameter> &objpara, shared_ptr<flowmap> &objphi);
  void calculate_ftle();
  vector3d<double> ftle() const { return ftle_; };
  vector3d<double> eig1() const { return eig1_; };
  vector3d<double> eig2() const { return eig2_; };
  vector3d<double> eig3() const { return eig3_; };
  vector4d<double> v1() const { return v1_; };
  vector4d<double> v2() const { return v2_; };
  vector4d<double> v3() const { return v3_; };

 private:
  double t;
  int nx;
  int ny;
  int nz;
  vector3d<double> phi_x;
  vector3d<double> phi_y;
  vector3d<double> phi_z;
  vector3d<double> phi_dudx;
  vector3d<double> phi_dudy;
  vector3d<double> phi_dudz;
  vector3d<double> phi_dvdx;
  vector3d<double> phi_dvdy;
  vector3d<double> phi_dvdz;
  vector3d<double> phi_dwdx;
  vector3d<double> phi_dwdy;
  vector3d<double> phi_dwdz;
  vector3d<double> ftle_;
  vector3d<double> eig1_;
  vector3d<double> eig2_;
  vector3d<double> eig3_;
  vector4d<double> v1_;
  vector4d<double> v2_;
  vector4d<double> v3_;
};

#endif
