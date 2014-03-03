#ifndef _lyapunov_
#define _lyapunov_

#include "ftle.h"
#include "flowmap.h"
#include "parameter.h"

using namespace std;

class lyapunov
{
public:
	lyapunov (parameter* objpara, flowmap* objphi);
	virtual ~lyapunov(void);
	double***  ftle() const { return ftle_; };
	double***  eig1() const { return eig1_; };
	double***  eig2() const { return eig2_; };
	double***  eig3() const { return eig3_; };
	double**** v1() const { return v1_; };
	double**** v2() const { return v2_; };
	double**** v3() const { return v3_; };
	
private:
        void Ecrire_Tecplot(const char* fichier);
        void Ecrire_Tecplot_Hessian(const char* fichier);
        void Ecrire_Tecplot_Binary(const char* fichier);
        FILE* fp;
        double t;
        int nnode;
        int nbelm;
        int nx;
        int ny;
        int nz;
        double   xmin;
        double   xmax;
        double   ymax;
        double   ymin;
        double   zmax;
        double   zmin;
        double*  coo;
        int* 	 cnc;
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
        double*** ftle_;
        double*** eig1_;
        double*** eig2_;
        double*** eig3_;
        double**** v1_;
        double**** v2_;
        double**** v3_;
};

#endif
