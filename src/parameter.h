#ifndef _parameter_
#define _parameter_

#include "ftle.h"

using namespace std;

class parameter {
 public:
  parameter(const char *fichier);
  double t() const { return t_; };
  double t0() const { return t0_; };
  double tolerance() const {return tolerance_;};
  double hmin() const { return hmin_; };
  double hmax() const { return hmax_; };

  int nx() const { return nx_; };
  int ny() const { return ny_; };
  int nz() const { return nz_; };
  double xmin() const { return xmin_; };
  double xmax() const { return xmax_; };
  double ymin() const { return ymin_; };
  double ymax() const { return ymax_; };
  double zmin() const { return zmin_; };
  double zmax() const { return zmax_; };
  string fu() const { return fu_; };
  string fv() const { return fv_; };
  string fw() const { return fw_; };
  string fdudx() const { return fdudx_; };
  string fdudy() const { return fdudy_; };
  string fdudz() const { return fdudz_; };
  string fdvdx() const { return fdvdx_; };
  string fdvdy() const { return fdvdy_; };
  string fdvdz() const { return fdvdz_; };
  string fdwdx() const { return fdwdx_; };
  string fdwdy() const { return fdwdy_; };
  string fdwdz() const { return fdwdz_; };

 private:
  double t0_;
  double t_;
  double tolerance_;
  double hmin_;
  double hmax_;
  int nx_;
  int ny_;
  int nz_;
  double xmin_;
  double xmax_;
  double ymin_;
  double ymax_;
  double zmin_;
  double zmax_;
  string fu_;
  string fv_;
  string fw_;
  string fdudx_;
  string fdudy_;
  string fdudz_;
  string fdvdx_;
  string fdvdy_;
  string fdvdz_;
  string fdwdx_;
  string fdwdy_;
  string fdwdz_;
};
#endif
