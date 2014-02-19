#include "lyapunov.h"
extern "C" {
    void dsyevd_(char& jobz, char& uplo, int& n, double* a, int& lda, double* w, double* work, int& lwork, int* iwork, int& liwork, int& info);
}

lyapunov::lyapunov(parameter* objpara, flowmap* objphi)
    :	T(objpara->get_T()),
        Nx(objpara->get_Nx()),
        Ny(objpara->get_Ny()),
        Nz(objpara->get_Nz()),
        phi_x(objphi->get_phi_x()),
        phi_y(objphi->get_phi_y()),
        phi_z(objphi->get_phi_z()),
        phi_dudx(objphi->get_phi_dudx()),
        phi_dudy(objphi->get_phi_dudy()),
        phi_dudz(objphi->get_phi_dudz()),
        phi_dvdx(objphi->get_phi_dvdx()),
        phi_dvdy(objphi->get_phi_dvdy()),
        phi_dvdz(objphi->get_phi_dvdz()),
        phi_dwdx(objphi->get_phi_dwdx()),
        phi_dwdy(objphi->get_phi_dwdy()),
        phi_dwdz(objphi->get_phi_dwdz())

{
cout << "Calculating FTLE" << endl;
// Valeur propre du tenseur

ftle = Construct3D(Nx, Ny, Nz);
eig1 = Construct3D(Nx, Ny, Nz);
eig2 = Construct3D(Nx, Ny, Nz);
eig3 = Construct3D(Nx, Ny, Nz);

V1 = Construct4D(Nx, Ny, Nz, 3);
V2 = Construct4D(Nx, Ny, Nz, 3);
V3 = Construct4D(Nx, Ny, Nz, 3);

// Lapack parameters
int n(2), lda(2), lwork(21), liwork(13), info(0); // iwork 5*n+3
char jobz = 'V';
char uplo = 'U';
lwork = 37; //2*n*n + 6*n +1
liwork = 18; // 5*n+3
            
// Calcul des coefficients du tenseur pour tous les noeuds du champ
for (int i(0); i<Nx; i++) {
    for (int j(0); j<Ny; j++) {
        for (int k(0); k<Nz; k++) {
            //  D1  D2  D3
            //  D2  D4  D5
            //  D3  D5  D6
            double D1 = phi_dudx[i][j][k]*phi_dudx[i][j][k] + phi_dvdx[i][j][k]*phi_dvdx[i][j][k] + phi_dwdx[i][j][k]*phi_dwdx[i][j][k];
            double D2 = phi_dudx[i][j][k]*phi_dudy[i][j][k] + phi_dvdx[i][j][k]*phi_dvdy[i][j][k] + phi_dwdx[i][j][k]*phi_dwdy[i][j][k];
            double D3 = phi_dudx[i][j][k]*phi_dudz[i][j][k] + phi_dvdx[i][j][k]*phi_dvdz[i][j][k] + phi_dwdx[i][j][k]*phi_dwdz[i][j][k];
            double D4 = phi_dudy[i][j][k]*phi_dudy[i][j][k] + phi_dvdy[i][j][k]*phi_dvdy[i][j][k] + phi_dwdy[i][j][k]*phi_dwdy[i][j][k];
            double D5 = phi_dudy[i][j][k]*phi_dudz[i][j][k] + phi_dvdy[i][j][k]*phi_dvdz[i][j][k] + phi_dwdy[i][j][k]*phi_dwdz[i][j][k];
            double D6 = phi_dudz[i][j][k]*phi_dudz[i][j][k] + phi_dvdz[i][j][k]*phi_dvdz[i][j][k] + phi_dwdz[i][j][k]*phi_dwdz[i][j][k];

            // Lapack parameter
            // and variable to store
            // output of function
	        int iwork[liwork];
	        double work[lwork];
	        double w[n]; // eigen value
            // Store collumn-wise!!!! this will
            // store the eigenvectors matrix
            double a[9] = {
                D1,  0.00, 0.00,
                D2,  D4,   0.00,
                D3,  D5,   D6
            };
            
            // Solve eigenproblem
			dsyevd_(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);
            if (info != 0)
                cout << "Error using Lapack dsyevd function" << endl;


            for (int m(0); m<n; m++) {
                V1[i][j][k][m] = a[m];
                V2[i][j][k][m] = a[n+m];
                V3[i][j][k][m] = a[2*n+m];
            }

            eig1[i][j][k] = w[0];
            eig2[i][j][k] = w[1];
            eig3[i][j][k] = w[2];
            
            // Lyapunov
            if (w[2] > std::numeric_limits<double>::epsilon())
                ftle[i][j][k]  = log(w[2])/2.0/fabs(T);
            else
                ftle[i][j][k] = 0;
            
        }
    }
}

}

lyapunov::~lyapunov(void)
{
Destruct3D(ftle);
Destruct3D(eig1);
Destruct3D(eig2);
Destruct3D(eig3);

Destruct4D(V1);
Destruct4D(V2);
Destruct4D(V3);
}

// Methode pour recuperer les donnees de la classe lyapunov
double*** lyapunov::get_Lyapunov(void)
{
   return ftle;
};

double*** lyapunov::get_Eig1(void)
{
   return eig1;
};

double*** lyapunov::get_Eig2(void)
{
   return eig2;
};

double*** lyapunov::get_Eig3(void)
{
   return eig3;
};

double**** lyapunov::get_V1(void)
{
   return V1;
};

double**** lyapunov::get_V2(void)
{
   return V2;
};

double**** lyapunov::get_V3(void)
{
   return V3;
};
