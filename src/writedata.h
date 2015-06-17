#ifndef _writedata_
#define _writedata_

#include "ftle.h"
#include "parameter.h"
#include "flowmap.h"
#include "lyapunov.h"

using namespace std;

class writedata
{
public:
   writedata(shared_ptr<parameter>& objpara, shared_ptr<flowmap>& objfm, shared_ptr<lyapunov>& objftle);
private:
  void write_ftle(std::string file);
  int nx;
	int ny;	
	int nz;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	double*** ftle;
};

#endif
