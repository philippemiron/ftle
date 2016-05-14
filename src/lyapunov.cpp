#include "lyapunov.h"
extern "C" {
void dsyevd_(char &jobz,
             char &uplo,
             int &n,
             double *a,
             int &lda,
             double *w,
             double *work,
             int &lwork,
             int *iwork,
             int &liwork,
             int &info);
}

lyapunov::lyapunov(const shared_ptr<parameter> &objpara, const shared_ptr<flowmap> &objfm) :
    t(objpara->t()),
    nx(objpara->nx()),
    ny(objpara->ny()),
    nz(objpara->nz()),
    x(objfm->x()),
    y(objfm->y()),
    z(objfm->z()),
    phi_x(objfm->phi_x()),
    phi_y(objfm->phi_y()),
    phi_z(objfm->phi_z()),
    phi_dudx(objfm->phi_dudx()),
    phi_dudy(objfm->phi_dudy()),
    phi_dudz(objfm->phi_dudz()),
    phi_dvdx(objfm->phi_dvdx()),
    phi_dvdy(objfm->phi_dvdy()),
    phi_dvdz(objfm->phi_dvdz()),
    phi_dwdx(objfm->phi_dwdx()),
    phi_dwdy(objfm->phi_dwdy()),
    phi_dwdz(objfm->phi_dwdz()) {

  vecResize(ftle_, nx, ny, nz);
  vecResize(eig1_, nx, ny, nz);
  vecResize(eig2_, nx, ny, nz);
  vecResize(eig3_, nx, ny, nz);
  vecResize(v1_, nx, ny, nz, 3);
  vecResize(v2_, nx, ny, nz, 3);
  vecResize(v3_, nx, ny, nz, 3);
};

void lyapunov::calculate_ftle() {
  // Lapack parameters
  int n(3), lda(n), lwork(1 + 6 * n + 2 * n * n), liwork(3 + 5 * n), info(0);
  char jobz = 'V';
  char uplo = 'U';

  // Calculate Cauchy-Green matrix at every node
  // Solve eingenvalue problems
  // Obtain Lyapunov exponent
#pragma omp parallel for
  for (int i = 0; i < nx; i++) {
    for (int j(0); j < ny; j++) {
      for (int k(0); k < nz; k++) {
        //  D1  D2  D3
        //  D2  D4  D5
        //  D3  D5  D6
        double D1 = phi_dudx[i][j][k] * phi_dudx[i][j][k] + phi_dvdx[i][j][k] * phi_dvdx[i][j][k]
            + phi_dwdx[i][j][k] * phi_dwdx[i][j][k];
        double D2 = phi_dudx[i][j][k] * phi_dudy[i][j][k] + phi_dvdx[i][j][k] * phi_dvdy[i][j][k]
            + phi_dwdx[i][j][k] * phi_dwdy[i][j][k];
        double D3 = phi_dudx[i][j][k] * phi_dudz[i][j][k] + phi_dvdx[i][j][k] * phi_dvdz[i][j][k]
            + phi_dwdx[i][j][k] * phi_dwdz[i][j][k];
        double D4 = phi_dudy[i][j][k] * phi_dudy[i][j][k] + phi_dvdy[i][j][k] * phi_dvdy[i][j][k]
            + phi_dwdy[i][j][k] * phi_dwdy[i][j][k];
        double D5 = phi_dudy[i][j][k] * phi_dudz[i][j][k] + phi_dvdy[i][j][k] * phi_dvdz[i][j][k]
            + phi_dwdy[i][j][k] * phi_dwdz[i][j][k];
        double D6 = phi_dudz[i][j][k] * phi_dudz[i][j][k] + phi_dvdz[i][j][k] * phi_dvdz[i][j][k]
            + phi_dwdz[i][j][k] * phi_dwdz[i][j][k];

        // Lapack parameter
        // and variable to store
        // output of function
        int iwork[liwork];
        double work[lwork];
        double w[n]; // eigen value
        // Store collumn-wise!!!! this will
        // store the eigenvectors matrix
        double a[9] = {
            D1, 0.00, 0.00,
            D2, D4, 0.00,
            D3, D5, D6
        };

        // Solve eigenproblem
        dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
        if (info != 0)
          cout << "Error using Lapack dsyevd function" << endl;

        for (int m(0); m < 3; m++) {
          v1_[i][j][k][m] = a[m];
          v2_[i][j][k][m] = a[n + m];
          v3_[i][j][k][m] = a[2 * n + m];
        }

        eig1_[i][j][k] = w[0];
        eig2_[i][j][k] = w[1];
        eig3_[i][j][k] = w[2];

        // Lyapunov exponent
        double eig_max = max(max(w[0], w[1]), w[2]);
        if (eig_max > numeric_limits<double>::epsilon())
          ftle_[i][j][k] = log(eig_max) / 2.0 / fabs(t);
        else
          ftle_[i][j][k] = 0;
      }
    }
  }
};

void lyapunov::output(string file) {
  ofstream myfile(file);
  if (myfile.is_open()) {
    myfile << "TITLE = \"ftle\"" << endl;
    myfile << "VARIABLES = x, y, z, xi1_x, xi1_y, xi1_z, xi2_x, xi2_y, xi2_z, xi3_x, xi3_y, xi3_z, L1, L2, L3, ftle"
        << endl;
    myfile << "zone I=" << nx << ", J=" << ny << ", K=" << nz << ", DATAPACKING=POINT" << endl;

    for (int k(0); k < nz; k++) {
      for (int j(0); j < ny; j++) {
        for (int i(0); i < nx; i++) {
          myfile << x[i] << " " << y[j] << " " << z[k] << " ";
          myfile << v1_[i][j][k][0] << " " << v1_[i][j][k][1] << " " << v1_[i][j][k][2] << " ";
          myfile << v2_[i][j][k][0] << " " << v2_[i][j][k][1] << " " << v2_[i][j][k][2] << " ";
          myfile << v3_[i][j][k][0] << " " << v3_[i][j][k][1] << " " << v3_[i][j][k][2] << " ";
          myfile << eig1_[i][j][k] << " " << eig2_[i][j][k] << " " << eig3_[i][j][k] << " " << ftle_[i][j][k] << endl;
        }
      }
    }
    // Closing the file
    myfile.close();
  } else {
    cerr << "Can't open output file name: " << file << endl;
    exit(-1);
  }
};

void lyapunov::outputVTK(string file) {
  ofstream myfile(file);
  if (myfile.is_open()) {
    myfile << "# vtk DataFile Version 2.0" << endl;
    myfile << "Solution" << endl;
    myfile << "ASCII" << endl << endl;

    myfile << "DATASET RECTILINEAR_GRID" << endl;
    myfile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;

    myfile << "X_COORDINATES " << nx << " float" << endl;
    for (int i(0); i < nx; i++)
      myfile << x[i] << " ";

    myfile << endl << "Y_COORDINATES " << ny << " float" << endl;
    for (int i(0); i < ny; i++)
      myfile << y[i] << " ";

    myfile << endl << "Z_COORDINATES " << nz << " float" << endl;
    for (int i(0); i < nz; i++)
      myfile << z[i] << " ";

    myfile << endl << "POINT_DATA " << nx * ny * nz << endl;

    myfile << "VECTORS " << "Eigenvector1" << " float" << endl;
    for (int k(0); k < nz; k++)
      for (int j(0); j < ny; j++)
        for (int i(0); i < nx; i++)
          myfile << v1_[i][j][k][0] << " " << v1_[i][j][k][1] << " " << v1_[i][j][k][2] << endl;

    myfile << "VECTORS " << "Eigenvector2" << " float" << endl;
    for (int k(0); k < nz; k++)
      for (int j(0); j < ny; j++)
        for (int i(0); i < nx; i++)
          myfile << v2_[i][j][k][0] << " " << v2_[i][j][k][1] << " " << v2_[i][j][k][2] << endl;

    myfile << "VECTORS " << "Eigenvector3" << " float" << endl;
    for (int k(0); k < nz; k++)
      for (int j(0); j < ny; j++)
        for (int i(0); i < nx; i++)
          myfile << v3_[i][j][k][0] << " " << v3_[i][j][k][1] << " " << v3_[i][j][k][2] << endl;

    myfile << "SCALARS " << "Eigenvalues" << " float" << " 3" << endl;
    myfile << "LOOKUP_TABLE default" << endl;
    for (int k(0); k < nz; k++)
      for (int j(0); j < ny; j++)
        for (int i(0); i < nx; i++)
          myfile << eig1_[i][j][k] << " " << eig2_[i][j][k] << " " << eig3_[i][j][k] << endl;

    myfile << "SCALARS " << "FTLE" << " float" << " 1" << endl;
    myfile << "LOOKUP_TABLE default" << endl;
    for (int k(0); k < nz; k++)
      for (int j(0); j < ny; j++)
        for (int i(0); i < nx; i++)
          myfile << ftle_[i][j][k] << endl;

    myfile.close();
  } else {
    cerr << "Can't open output file name: " << file << "." << endl;
    exit(-1);
  }
}