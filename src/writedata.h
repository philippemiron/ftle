#ifndef _writedata_
#define _writedata_

#include "ftle.h"
#include "parameter.h"
#include "flowmap.h"
#include "lyapunov.h"

class writedata
{
public:
   writedata(parameter* objpara, flowmap* objfm, lyapunov* objftle);

private:
  void write_ftle(std::string file);
  int nx;
	int ny;	
	int nz;
	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> z;
	double*** ftle;
};

#endif
