#ifndef _parameter_
#define _parameter_

#include "ftle.h"
#include "enum.h"

using namespace std;

class parameter
{
public:
  parameter(const char* fichier);
  double	t() const { return t_; };
  int	npt() const { return npt_; };
  int	t0() const { return t0_; };
  int	nx() const { return nx_; };
  int	ny() const { return ny_; };
  int	nz() const { return nz_; };
  _MthOde mthode() const { return mthode_; };
  double  xmin() const { return xmin_; };
  double  xmax() const { return xmax_; };
  double  ymin() const { return ymin_; };
  double  ymax() const { return ymax_; };
  double  zmin() const { return zmin_; };
  double  zmax() const { return zmax_; };    
  std::string	filecg() const { return filecg_; };
  std::string	fu() const { return fu_; };
  std::string	fv() const { return fv_; };
  std::string	fw() const { return fw_; };  
  std::string	fdudx() const { return fdudx_; };
  std::string	fdudy() const { return fdudy_; };
  std::string	fdudz() const { return fdudz_; };  
  std::string	fdvdx() const { return fdvdx_; };
  std::string	fdvdy() const { return fdvdy_; };
  std::string	fdvdz() const { return fdvdz_; };  
  std::string	fdwdx() const { return fdwdx_; };
  std::string	fdwdy() const { return fdwdy_; };
  std::string	fdwdz() const { return fdwdz_; };  

private:
  double  t_;
  int   npt_;
  int   t0_;
  int   nx_;
  int   ny_;
  int   nz_;
  double xmin_;
  double xmax_;
  double ymin_;
  double ymax_;
  double zmin_;
  double zmax_;
  _MthOde  mthode_;
  std::string filecg_;
  std::string fu_;
  std::string fv_;
  std::string fw_;
  std::string fdudx_;
  std::string fdudy_;
  std::string fdudz_;
  std::string fdvdx_;
  std::string fdvdy_;
  std::string fdvdz_;
  std::string fdwdx_;
  std::string fdwdy_;
  std::string fdwdz_;
};
#endif
