#ifndef _parameter_
#define _parameter_

#include "ftle.h"
#include "enum.h"

class parameter
{
public:
    parameter(const char* fichier);
    virtual ~parameter(void);
    double	get_T() const { return T; };
    int	get_Npt() const { return Npt; };
    int	get_T0() const { return T0; };
    int	get_Nx() const { return Nx; };
    int	get_Ny() const { return Ny; };
    int	get_Nz() const { return Nz; };
    _MthOde  get_MthOde() const { return MthOde; };
    double  getXmin() const { return xmin; };
    double  getXmax() const { return xmax; };
    double  getYmin() const { return ymin; };
    double  getYmax() const { return ymax; };
    double  getZmin() const { return zmin; };
    double  getZmax() const { return zmax; };    
    std::string	getOutput_File_CG() const { return output_file_cg; };
    std::string	getFunctionU() const { return fu; };
    std::string	getFunctionV() const { return fv; };
    std::string	getFunctionW() const { return fw; };  
    std::string	getFunctionDudx() const { return fdudx; };
    std::string	getFunctionDudy() const { return fdudy; };
    std::string	getFunctionDudz() const { return fdudz; };  
    std::string	getFunctionDvdx() const { return fdvdx; };
    std::string	getFunctionDvdy() const { return fdvdy; };
    std::string	getFunctionDvdz() const { return fdvdz; };  
    std::string	getFunctionDwdx() const { return fdwdx; };
    std::string	getFunctionDwdy() const { return fdwdy; };
    std::string	getFunctionDwdz() const { return fdwdz; };  

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
