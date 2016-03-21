#include "writedata.h"

writedata::writedata(shared_ptr<parameter> &objpara, shared_ptr<flowmap> &objfm, shared_ptr<lyapunov> &objftle) :
    nx(objpara->nx()),
    ny(objpara->ny()),
    nz(objpara->nz()),
    x(objfm->x()),
    y(objfm->y()),
    z(objfm->z()),
    eig1(objftle->eig1()),
    eig2(objftle->eig2()),
    eig3(objftle->eig3()),
    v1(objftle->v1()),
    v2(objftle->v2()),
    v3(objftle->v3()),
    ftle(objftle->ftle()) {
  write_ftle("ftle.dat");
};

void writedata::write_ftle(std::string file) {

//Opening the file
  std::ofstream myfile(file.c_str());
  if (myfile.is_open()) {
    // TP Headlines
    myfile << "TITLE = \"ftle\"" << endl;
    myfile << "VARIABLES = x, y, z, xi1_x, xi1_y, xi1_z, xi2_x, xi2_y, xi2_z, xi3_x, xi3_y, xi3_z, L1, L2, L3, ftle"
        << endl;
    myfile << "zone I=" << nx << ", J=" << ny << ", K=" << nz << ", DATAPACKING=POINT" << endl;

    for (int k(0); k < nz; k++) {
      for (int j(0); j < ny; j++) {
        for (int i(0); i < nx; i++) {
          myfile << x[i] << " " << y[j] << " " << z[k] << " ";
          myfile << v1[i][j][k][0] << " " << v1[i][j][k][1] << " " << v1[i][j][k][2] << " ";
          myfile << v2[i][j][k][0] << " " << v2[i][j][k][1] << " " << v2[i][j][k][2] << " ";
          myfile << v3[i][j][k][0] << " " << v3[i][j][k][1] << " " << v3[i][j][k][2] << " ";
          myfile << eig1[i][j][k] << " " << eig2[i][j][k] << " " << eig3[i][j][k] << " " << ftle[i][j][k] << endl;
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


