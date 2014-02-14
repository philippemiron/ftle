// Construction and destruction of 2D, 3D and 4D matrices.
#include "matrix.h"
/////////////////////
// Create a 2D matrix
double** Construct2D(int Nx, int Ny) {
	
double* ptr   = new double  [Nx*Ny];
double** x = new double* [Nx];
for(int i=0;i<Nx;i++,ptr+=Ny) 
	x[i]=ptr;

return x;
}

// Delete a 2D matrix
void Destruct2D(double** x) {
	delete [] x[0];
	delete [] x;
}

/////////////////////
// Create a 3D matrix
double*** Construct3D(int Nx, int Ny, int Nz) {

double* ptr = new double [Nx*Ny*Nz];
double** m = new double* [Nx*Ny];
double*** x = new double** [Nx];
for (int i=0;i<Nx;i++,m+=Ny) {
	x[i]=m;
	for(int j=0;j<Ny;j++,ptr+=Nz) {
		m[j]=ptr;
	}
}
return x;
}

// Delete a 3D matrix
void Destruct3D(double*** x) {
	delete [] x[0][0];
	delete [] x[0];
	delete [] x;
}

/////////////////////
// Create a 4D matrix
double**** Construct4D(int Nx, int Ny, int Nz, int N) {

double*  ptr = new double [Nx*Ny*Nz*N];
double**   m = new double* [Nx*Ny*Nz];
double***  n = new double** [Nx*Ny];
double**** x = new double*** [Nx];

for (int i=0;i<Nx;i++,n+=Ny) {
    x[i]=n;
    for (int j=0;j<Ny;j++,m+=Nz) {
        n[j]=m;
        for(int k=0;k<Nz;k++,ptr+=N) {
                m[k]=ptr;
        }
    }
}
return x;
}

// Delete a 4D matrix
void Destruct4D(double**** x) {
	delete [] x[0][0][0];
    delete [] x[0][0];
    delete [] x[0];
    delete [] x;
}

/////////////////////
// Create a 5D matrix
double***** Construct5D(int N1, int N2, int N3, int N4, int N5) {

double*  ptr = new double [N1*N2*N3*N4*N5];
double**   m = new double* [N1*N2*N3*N4];
double***  n = new double** [N1*N2*N3];
double**** o = new double*** [N1*N2];
double***** x = new double**** [N1];

for (int i=0;i<N1;i++,o+=N2) {
    x[i]=o;
    for (int j=0;j<N2;j++,n+=N3) {
        o[j]=n;
        for (int k=0;k<N3;k++,m+=N4) {
            n[k]=m;
            for(int l=0;l<N4;l++,ptr+=N5) {
                m[l]=ptr;
            }
        }
    }
}
return x;
}

// Delete a 5D matrix
void Destruct5D(double***** x) {
	delete [] x[0][0][0][0];
	delete [] x[0][0][0];
    delete [] x[0][0];
    delete [] x[0];
    delete [] x;
}
