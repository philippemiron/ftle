#include "writedata.h"

writedata::writedata(shared_ptr<parameter>& objpara, shared_ptr<flowmap>& objfm, shared_ptr<lyapunov>& objftle):
  nx(objpara->nx()),
	ny(objpara->ny()),
	nz(objpara->nz()),
	x(objfm->x()),
	y(objfm->y()),
	z(objfm->z()),
	ftle(objftle->ftle())
{
	write_ftle("ftle.dat"); 
};

void writedata::write_ftle(std::string file) 
{
		
//Opening the file
std::ofstream myfile(file.c_str());
if (myfile.is_open()) {
	// TP Headlines
  myfile << "TITLE = \"ftle\"" << std::endl << "VARIABLES = x, y, z, ftle" << std::endl;
  myfile << "zone I=" << nx << ", J=" << ny << ", K=" << nz << ", DATAPACKING=POINT" << std::endl;

  for(int k(0); k<nz; k++) {
    for(int j(0); j<ny; j++) {
      for(int i(0); i<nx; i++) {
        myfile << x[i] << " " << y[j] << " " << z[k] << " " << ftle[i][j][k] << std::endl;
      }
    }
  }
	// Closing the file
	myfile.close();
} else {
	std::cout << "Can't open output file name: " << file << std::endl;
	exit(0);
}

};


