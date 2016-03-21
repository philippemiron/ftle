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

lyapunov::lyapunov(shared_ptr<parameter> &objpara, shared_ptr<flowmap> &objphi) :
    t(objpara->t()),
    nx(objpara->nx()),
    ny(objpara->ny()),
    nz(objpara->nz()),
    phi_x(objphi->phi_x()),
    phi_y(objphi->phi_y()),
    phi_z(objphi->phi_z()),
    phi_dudx(objphi->phi_dudx()),
    phi_dudy(objphi->phi_dudy()),
    phi_dudz(objphi->phi_dudz()),
    phi_dvdx(objphi->phi_dvdx()),
    phi_dvdy(objphi->phi_dvdy()),
    phi_dvdz(objphi->phi_dvdz()),
    phi_dwdx(objphi->phi_dwdx()),
    phi_dwdy(objphi->phi_dwdy()),
    phi_dwdz(objphi->phi_dwdz()) {
  cout << "Calculating FTLE" << endl;

  vecResize(ftle_, nx, ny, nz);
  vecResize(eig1_, nx, ny, nz);
  vecResize(eig2_, nx, ny, nz);
  vecResize(eig3_, nx, ny, nz);
  vecResize(v1_, nx, ny, nz, 3);
  vecResize(v2_, nx, ny, nz, 3);
  vecResize(v3_, nx, ny, nz, 3);

  // Lapack parameters
  int n(3), lda(n), lwork(1 + 6*n + 2*n*n), liwork(3 + 5*n), info(0);
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
        double D1 = phi_dudx[i][j][k]*phi_dudx[i][j][k] + phi_dvdx[i][j][k]*phi_dvdx[i][j][k]
            + phi_dwdx[i][j][k]*phi_dwdx[i][j][k];
        double D2 = phi_dudx[i][j][k]*phi_dudy[i][j][k] + phi_dvdx[i][j][k]*phi_dvdy[i][j][k]
            + phi_dwdx[i][j][k]*phi_dwdy[i][j][k];
        double D3 = phi_dudx[i][j][k]*phi_dudz[i][j][k] + phi_dvdx[i][j][k]*phi_dvdz[i][j][k]
            + phi_dwdx[i][j][k]*phi_dwdz[i][j][k];
        double D4 = phi_dudy[i][j][k]*phi_dudy[i][j][k] + phi_dvdy[i][j][k]*phi_dvdy[i][j][k]
            + phi_dwdy[i][j][k]*phi_dwdy[i][j][k];
        double D5 = phi_dudy[i][j][k]*phi_dudz[i][j][k] + phi_dvdy[i][j][k]*phi_dvdz[i][j][k]
            + phi_dwdy[i][j][k]*phi_dwdz[i][j][k];
        double D6 = phi_dudz[i][j][k]*phi_dudz[i][j][k] + phi_dvdz[i][j][k]*phi_dvdz[i][j][k]
            + phi_dwdz[i][j][k]*phi_dwdz[i][j][k];

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
          v3_[i][j][k][m] = a[2*n + m];
        }

        eig1_[i][j][k] = w[0];
        eig2_[i][j][k] = w[1];
        eig3_[i][j][k] = w[2];

        // Lyapunov exponent
        double eig_max = max(max(w[0], w[1]), w[2]);
        if (eig_max > numeric_limits<double>::epsilon())
          ftle_[i][j][k] = log(eig_max)/2.0/fabs(t);
        else
          ftle_[i][j][k] = 0;
      }
    }
  }
}