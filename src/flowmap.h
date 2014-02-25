#ifndef _flowmap_
#define _flowmap_

#include "ftle.h"
#include "parameter.h"
#include "matrix.h"

using namespace std;

class flowmap
{
public:
	flowmap (parameter* objpara);
	virtual ~flowmap(void);
	double*** get_phi_x() const { return phi_x; };
	double*** get_phi_y() const { return phi_x; };
	double*** get_phi_z() const { return phi_x; };
	double*** get_phi_dudx() const { return phi_dudx; };
	double*** get_phi_dudy() const { return phi_dudy; };
	double*** get_phi_dudz() const { return phi_dudz; };
	double*** get_phi_dvdx() const { return phi_dvdx; };
	double*** get_phi_dvdy() const { return phi_dvdy; };
	double*** get_phi_dvdz() const { return phi_dvdz; };
	double*** get_phi_dwdx() const { return phi_dwdx; };
	double*** get_phi_dwdy() const { return phi_dwdy; };
	double*** get_phi_dwdz() const { return phi_dwdz; };
	std::vector<double> get_x() const { return x; };
	std::vector<double> get_y() const { return y; };
	std::vector<double> get_z() const { return z; };

private:
	_MthOde ode;
    int count_particules_ext;
	int Nx;
	int Ny;
	int Nz;
	void initialize_expression(int number_of_threads);
	void interpolation (double t, double* y, double* g);
	void (flowmap::*ptr_ode) (double T, int Npt, double TIni, vector<double>& D);
	void int_euler (double T, int Npt, double TIni, vector<double>& D);
	void int_rk5 (double T, int Npt, double TIni, vector<double>& D);
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
    std::string function_u;
    std::string function_v;
    std::string function_w;
    std::string function_dudx;
    std::string function_dudy;
    std::string function_dudz;
    std::string function_dvdx;
    std::string function_dvdy;
    std::string function_dvdz;
    std::string function_dwdx;
    std::string function_dwdy;
    std::string function_dwdz;
	
	double*** phi_x;
	double*** phi_y;
	double*** phi_z;
	double*** phi_dudx;
	double*** phi_dudy;
	double*** phi_dudz;
	double*** phi_dvdx;
	double*** phi_dvdy;
	double*** phi_dvdz;
	double*** phi_dwdx;
	double*** phi_dwdy;
	double*** phi_dwdz;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	double T;
	int   Npt;
	int   T0;
};

#endif
