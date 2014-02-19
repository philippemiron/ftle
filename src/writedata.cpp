#include "writedata.h"

writedata::writedata(parameter* objpara, flowmap* objfm, lyapunov* objftle):
    nx(objpara->get_Nx()),
	ny(objpara->get_Ny()),
	nz(objpara->get_Nz()),
	x(objfm->get_x()),
	y(objfm->get_y()),
	z(objfm->get_z()),
	ftle(objftle->get_Lyapunov())
{
	write_ftle("ftle.dat"); 
};
writedata::~writedata(void)
{
};

void writedata::write_ftle(std::string file) 
{
		
//Opening the file
std::ofstream myfile(file.c_str());
if (myfile.is_open())
{

	// TP Headlines
    myfile << "TITLE = \"ftle\"" << std::endl << "VARIABLES = x, y, z, ftle" << std::endl;
    myfile << "zone I=" << nx << ", J=" << ny << ", K=" << nz << ", DATAPACKING=POINT" << std::endl;

    for(int k(0); k<nz; k++) {
        for(int j(0); j<ny; j++) {
            for(int i(0); i<nx; i++) {
                myfile << x[i] << " " << y[j] << " " << z[k] << ftle[k][j][i] << std::endl;
            }
        }
    }
	// Closing the file
	myfile.close();
} 
else 
{
	std::cout << "Can't open output file name: " << file << std::endl;
	exit(0);
}

};


