#ifndef _parameter_
#define _parameter_

#include "ftle.h"
#include "enum.h"

class parameter
{
public:
    parameter(const char* fichier);
    virtual ~parameter(void);
    double	get_T();
    int	get_Npt();
    int	get_T0();
    int	get_TFin();
    int	get_Nx();
    int	get_Ny();
    int	get_Nz();
    _MthOde  get_MthOde();
    double  getXmin();
    double  getXmax();
    double  getYmin();
    double  getYmax();
    double  getZmin();
    double  getZmax();    
    std::string	getOutput_File_CG();
    std::string	getFunctionU();
    std::string	getFunctionV();
    std::string	getFunctionW();  
    std::string	getFunctionDudx();
    std::string	getFunctionDudy();
    std::string	getFunctionDudz();  
    std::string	getFunctionDvdx();
    std::string	getFunctionDvdy();
    std::string	getFunctionDvdz();  
    std::string	getFunctionDwdx();
    std::string	getFunctionDwdy();
    std::string	getFunctionDwdz();  

private:
    double  T;
    int   Npt;
    int   T0;
    int   Nx;
    int   Ny;
    int   Nz;
    double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
    _MthOde  MthOde;
    char fichier_piv[80];
    std::string output_file_cg;
    std::string fu;
    std::string fv;
    std::string fw;
    std::string fdudx;
    std::string fdudy;
    std::string fdudz;
    std::string fdvdx;
    std::string fdvdy;
    std::string fdvdz;
    std::string fdwdx;
    std::string fdwdy;
    std::string fdwdz;
};
#endif
