#ifndef _flowmap_
#define _flowmap_

#include "ftle.h"
#include "parameter.h"
#include "matrix.h"

using namespace std;

class flowmap
{
public:
	flowmap(shared_ptr<parameter>& objpara); 
	virtual ~flowmap();
	double*** phi_x() const { return phi_x_; };
	double*** phi_y() const { return phi_y_; };
	double*** phi_z() const { return phi_z_; };
	double*** phi_dudx() const { return phi_dudx_; };
	double*** phi_dudy() const { return phi_dudy_; };
	double*** phi_dudz() const { return phi_dudz_; };
	double*** phi_dvdx() const { return phi_dvdx_; };
	double*** phi_dvdy() const { return phi_dvdy_; };
	double*** phi_dvdz() const { return phi_dvdz_; };
	double*** phi_dwdx() const { return phi_dwdx_; };
	double*** phi_dwdy() const { return phi_dwdy_; };
	double*** phi_dwdz() const { return phi_dwdz_; };
	vector<double> x() const { return x_; };
	vector<double> y() const { return y_; };
	vector<double> z() const { return z_; };

private:
	_MthOde ode;
    int count_particules_ext;
	int nx;
	int ny;
	int nz;
	void initialize_expression(int number_of_threads);
	void interpolation (double t, double* y, double* g);
	void (flowmap::*ptr_ode) (double t, int npt, double tini, vector<double>& d);
	void int_euler (double t, int npt, double tini, vector<double>& d);
	void int_rk5 (double t, int npt, double tini, vector<double>& d);
	exprtk::expression<double>* exp_u;
  exprtk::expression<double>* exp_v;
  exprtk::expression<double>* exp_w;
  exprtk::expression<double>* exp_dudx;
  exprtk::expression<double>* exp_dudy;
  exprtk::expression<double>* exp_dudz;
  exprtk::expression<double>* exp_dvdx;
  exprtk::expression<double>* exp_dvdy;
  exprtk::expression<double>* exp_dvdz;
  exprtk::expression<double>* exp_dwdx;
  exprtk::expression<double>* exp_dwdy;
  exprtk::expression<double>* exp_dwdz;
	vector<double> var_x;
  vector<double> var_y;
  vector<double> var_z;
  vector<double> var_t;
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
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	double t;
	int   npt;
	int   t0;

	double*** phi_x_;
	double*** phi_y_;
	double*** phi_z_;
	double*** phi_dudx_;
	double*** phi_dudy_;
	double*** phi_dudz_;
	double*** phi_dvdx_;
	double*** phi_dvdy_;
	double*** phi_dvdz_;
	double*** phi_dwdx_;
	double*** phi_dwdy_;
	double*** phi_dwdz_;
	vector<double> x_;
	vector<double> y_;
	vector<double> z_;
};

#endif
