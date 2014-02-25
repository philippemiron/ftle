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
	double***  get_Lyapunov() const { return ftle; };
	double***  get_Eig1() const { return eig1; };
	double***  get_Eig2() const { return eig2; };
	double***  get_Eig3() const { return eig3; };
	double**** get_V1() const { return V1; };
	double**** get_V2() const { return V2; };
	double**** get_V3() const { return V3; };
	
private:
        void Ecrire_Tecplot(const char* fichier);
        void Ecrire_Tecplot_Hessian(const char* fichier);
        void Ecrire_Tecplot_Binary(const char* fichier);
        FILE* fp;
        double T;
        int nnode;
        int nbelm;
        int NbChamps;
        int Nx;
        int Ny;
        int Nz;
        double   XMin;
        double   XMax;
        double   YMax;
        double   YMin;
        double   ZMax;
        double   ZMin;
        double* coo;
        int* cnc;
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
        double*** ftle;
        double*** eig1;
        double*** eig2;
        double*** eig3;
        double**** V1;
        double**** V2;
        double**** V3;
};

#endif
