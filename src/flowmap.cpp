#include "flowmap.h"

extern "C" {
    void dcopy_(int& n, double* x, int& incx, double* y, int& incy);
    void daxpy_(int& n, double& alpha, double* x, int& incx, double* y, int& incy);  
}

flowmap::flowmap(shared_ptr<parameter>& objpara):
  ode(objpara->mthode()),
  nx(objpara->nx()),
  ny(objpara->ny()),
  nz(objpara->nz()),
  xmin(objpara->xmin()),
  xmax(objpara->xmax()),
  ymin(objpara->ymin()),
  ymax(objpara->ymax()),
  zmin(objpara->zmin()),
  zmax(objpara->zmax()),
  function_u(objpara->fu()),
  function_v(objpara->fv()),
  function_w(objpara->fw()),
  function_dudx(objpara->fdudx()),
  function_dudy(objpara->fdudy()),
  function_dudz(objpara->fdudz()),
  function_dvdx(objpara->fdvdx()),
  function_dvdy(objpara->fdvdy()),
  function_dvdz(objpara->fdvdz()),
  function_dwdx(objpara->fdwdx()),
  function_dwdy(objpara->fdwdy()),
  function_dwdz(objpara->fdwdz()),
  t(objpara->t()),
  npt(objpara->npt()),
  t0(objpara->t0())
{

cout << "Flow map integrations (T=" << t << "s) for the " << nx*ny*nz << " particles." << endl;

// Choix de la méthode d'intégration numérique utilisée
switch(ode)
{
	case euler: ptr_ode=&flowmap::int_euler;
				cout << "	Using euler integration" << endl;
				break;
	case rk5:   ptr_ode=&flowmap::int_rk5;
				cout << "	Using RK5 integration" << endl;
				break;
	default:    ptr_ode=&flowmap::int_euler;
				cout << "	Using euler integration by default" << endl;
				break;
} 
cout << "	Using analytical equation" << endl;	
int number_of_threads = omp_get_max_threads();

// Setting up the expression evaluation
initialize_expression(number_of_threads);

// Formation des matrices de deplacements.
x_.resize(nx, 0.0);
y_.resize(ny, 0.0);
z_.resize(nz, 0.0);

for (int i(0); i<nx; i++)
    x_[i] = xmin + i*(xmax-xmin)/(nx-1);

for (int i(0); i<ny; i++)
    y_[i] = ymin + i*(ymax-ymin)/(ny-1);
    
for (int i(0); i<nz; i++)
    z_[i] = zmin + i*(zmax-zmin)/(nz-1);

phi_x_ = Construct3D(nx, ny, nz);
phi_y_ = Construct3D(nx, ny, nz);
phi_z_ = Construct3D(nx, ny, nz);

phi_dudx_ = Construct3D(nx, ny, nz);
phi_dudy_ = Construct3D(nx, ny, nz);
phi_dudz_ = Construct3D(nx, ny, nz);
phi_dvdx_ = Construct3D(nx, ny, nz);
phi_dvdy_ = Construct3D(nx, ny, nz);
phi_dvdz_ = Construct3D(nx, ny, nz);
phi_dwdx_ = Construct3D(nx, ny, nz);
phi_dwdy_ = Construct3D(nx, ny, nz);
phi_dwdz_ = Construct3D(nx, ny, nz);


// Pour déterminer combien de particules sortent du domaine
// lors du calcul des déplacements
                
// 1- Faire une boucle sur les noeuds du maillage
// 2- Trouver les coordonnees du noeud courant
// 3- Calculer le delacement Lagrangien du noeud  
//    a l'aide des vitesses mesurees par PIV
#pragma omp parallel for
for (int i=0; i<nx;i++) {
    for (int j=0; j<ny; j++) { 
        for (int k=0; k<nz; k++) {
            // deplacement du point courant
            vector<double> d(30,0);
            // Recupere les coordonnées du noeuds courants
            d[0]=x_[i];
            d[1]=y_[j];
            d[2]=z_[k];                

            // Les paramètres de l'intégration du déplacement sont:
            // Nous intégrons de t0 a t0+T     // T: Temps d'integration des particules
            // npt: Nombre de points (subdivisions) lors de l'integration
            // t0: Temps t0 de notre intégration
            // d: vecteurs de déplacements
            (this->*ptr_ode) (t, npt, t0, d);

            phi_x_[i][j][k] = d[0];
            phi_y_[i][j][k] = d[1];
            phi_z_[i][j][k] = d[2];

            phi_dudx_[i][j][k] = d[3];
            phi_dudy_[i][j][k] = d[4];
            phi_dudz_[i][j][k] = d[5];
            phi_dvdx_[i][j][k] = d[6];
            phi_dvdy_[i][j][k] = d[7];
            phi_dvdz_[i][j][k] = d[8];
            phi_dwdx_[i][j][k] = d[9];
            phi_dwdy_[i][j][k] = d[10];
            phi_dwdz_[i][j][k] = d[11];
        }
    }
}

}

void flowmap::int_euler (double T, int npt, double t0, vector<double>& d) {

int  nd = 12;
int incx = 1;
double u[nd];
double x0[nd];

// ==========
x0[0]=d[0];   // x_0
x0[1]=d[1];   // y_0
x0[2]=d[2];   // z_0
// CI Derivee premiere
x0[3]=1.0;      
x0[4]=0.0;      
x0[5]=0.0;      
x0[6]=0.0;      
x0[7]=1.0;		
x0[8]=0.0;
x0[9]=0.0;
x0[10]=0.0;
x0[11]=1.0;

// Nous intégrons de 0 a T
// Si T est négatif dt sera aussi négatif
double dt = t/(npt-1);

for (int k(0); k<npt; k++) {
  // velocity interpolation
  double t  = t0 + k*dt;
  interpolation(t,x0,u);

  // x0 += u*dt
	daxpy_(nd, dt, u, incx, x0, incx);
}
dcopy_(nd, x0, incx, &d[0], incx);
};

void flowmap::int_rk5 (double t, int npt, double t0, vector<double>& d)
{
int  nd = 12;
double g[nd];

double k1[nd];
double k2[nd];
double k3[nd];
double k4[nd];
double k5[nd];
double k6[nd];

double y1[nd];
double y2[nd];
double y3[nd];
double y4[nd];
double y5[nd];
double y6[nd];

// ==========
y1[0]=d[0];   // x_0
y1[1]=d[1];   // y_0
y1[2]=d[2];   // z_0

// CI Derivee premiere
y1[3]=1.0;      
y1[4]=0.0;      
y1[5]=0.0;      
y1[6]=0.0;      
y1[7]=1.0;		
y1[8]=0.0;
y1[9]=0.0;
y1[10]=0.0;
y1[11]=1.0;

// Si T est négatif, dt sera aussi négatif et nous intégrons de t0 à T qui se situe avant t0
double dt = t/(npt-1);

for(int k=0; k<npt-1; k++) {
	// Détermine les différents temps nécessaire pour l'intégration 
	double t  = t0 + k*dt;
	double t1 = t + dt*0.25;
	double t2 = t + dt*0.5 ;
	double t3 = t + dt*0.75;
	double t4 = t + dt*1.0;
	//cout << "k: " << k << " t: " << t << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << " t4: " << t4 << endl;
	// ==============================================
	interpolation(t,y1,g); 
	for(int i=0; i<nd; i++) {
		 k1[i] = g[i]*dt;
		 y2[i] = y1[i] + 0.25*k1[i];
	}
	// ==============================================
	interpolation(t1,y2,g); 
	for(int i=0; i<nd; i++) {
    k2[i] = g[i]*dt;
    y3[i] = y1[i] + (k1[i]+k2[i])/8.0;
	}
	// ==============================================
	interpolation(t1,y3,g); 
	for(int i=0; i<nd; i++) {
    k3[i] = g[i]*dt;
    y4[i] = y1[i] - k2[i]/2.0 + k3[i];
	}
	// ==============================================
	interpolation(t2,y4,g); 
	for(int i=0; i<nd; i++) {
    k4[i] = g[i]*dt;
    y5[i] = y1[i] + (3.0*k1[i] + 9.0*k4[i])/16.0;
	}
	// ==============================================
	interpolation(t3,y5,g); 
	for(int i=0; i<nd; i++) {
    k5[i] = g[i]*dt;
    y6[i] = y1[i] + (-3.0*k1[i] + 2.0*k2[i] + 12.0*k3[i] - 12.0*k4[i] + 8.0*k5[i])/7.0;
	}
	// ==============================================
	interpolation(t4,y6,g); 
	for(int i=0; i<nd; i++) {
    k6[i] = g[i]*dt;
    y6[i] = y1[i] + (7.0*k1[i] + 32.0*k3[i] + 12.0*k4[i] + 32.0*k5[i] + 7.0*k6[i])/90.0;
	}
	for(int i=0;i<nd;i++) y1[i]=y6[i];
}
for(int i=0;i<nd;i++) d[i]=y1[i];
}

void flowmap::interpolation (double t, double* y, double* g)
{
  //cout << "x:" << y[0] << " y: " << y[1] << endl;
  double f[12];

  // one variable and one expression
  // for every thread
  int tid = omp_get_thread_num();
  var_x[tid] = y[0];
  var_y[tid] = y[1];
  var_z[tid] = y[2];
  var_t[tid] = t;

  // Solution analytique des vitesses et des dérivées
  // u,v,w
  f[0] = exp_u[tid].value();
  f[1] = exp_v[tid].value();
  f[2] = exp_w[tid].value();

  // dudx dudy, dudz
  f[3] = exp_dudx[tid].value();
  f[4] = exp_dudy[tid].value();
  f[5] = exp_dudz[tid].value();

  // dvdx dvdy, dvdz
  f[6] = exp_dvdx[tid].value();
  f[7] = exp_dvdy[tid].value();
  f[8] = exp_dvdz[tid].value();

  // dwdx dwdy, dwdz
  f[9] =  exp_dwdx[tid].value();
  f[10] = exp_dwdy[tid].value();
  f[11] = exp_dwdz[tid].value();

  // Vecteur pour le calcul du deplacement
  g[0] = f[0];
  g[1] = f[1];
  g[2] = f[2];

  // Derivees premieres
  g[3] = f[3]*y[3] + f[4]*y[6] + f[5]*y[9];
  g[4] = f[3]*y[4] + f[4]*y[7] + f[5]*y[10];
  g[5] = f[3]*y[5] + f[4]*y[8] + f[5]*y[11];
  g[6] = f[6]*y[3] + f[7]*y[6] + f[8]*y[9];
  g[7] = f[6]*y[4] + f[7]*y[7] + f[8]*y[10];
  g[8] = f[6]*y[5] + f[7]*y[8] + f[8]*y[11];
  g[9] = f[9]*y[3] + f[10]*y[6]+ f[11]*y[9];
  g[10] = f[9]*y[4] + f[10]*y[7] + f[11]*y[10];
  g[11] = f[9]*y[5] + f[10]*y[8] + f[11]*y[11];
}

void flowmap::initialize_expression(int number_of_threads) {
    
  // one expression per threads
  exp_u = new exprtk::expression<double>[number_of_threads];
  exp_v = new exprtk::expression<double>[number_of_threads];
  exp_w = new exprtk::expression<double>[number_of_threads];
  exp_dudx = new exprtk::expression<double>[number_of_threads];
  exp_dudy = new exprtk::expression<double>[number_of_threads];
  exp_dudz = new exprtk::expression<double>[number_of_threads];
  exp_dvdx = new exprtk::expression<double>[number_of_threads];
  exp_dvdy = new exprtk::expression<double>[number_of_threads];
  exp_dvdz = new exprtk::expression<double>[number_of_threads];
  exp_dwdx = new exprtk::expression<double>[number_of_threads];
  exp_dwdy = new exprtk::expression<double>[number_of_threads];
  exp_dwdz = new exprtk::expression<double>[number_of_threads];

  // one variable per threads
  var_x.resize(number_of_threads, 0.0);
  var_y.resize(number_of_threads, 0.0);
  var_z.resize(number_of_threads, 0.0);
  var_t.resize(number_of_threads, 0.0);

  // initialization of all parameters
  for (int i(0); i<number_of_threads; i++) {
    exprtk::symbol_table<double> symbol_table;
    symbol_table.add_variable("x", var_x[i]);
    symbol_table.add_variable("y", var_y[i]);
    symbol_table.add_variable("z", var_z[i]);
    symbol_table.add_variable("t", var_t[i]);
    symbol_table.add_constants();

    exp_u[i].register_symbol_table(symbol_table);
    exp_v[i].register_symbol_table(symbol_table);    
    exp_w[i].register_symbol_table(symbol_table);
    exp_dudx[i].register_symbol_table(symbol_table);
    exp_dudy[i].register_symbol_table(symbol_table);
    exp_dudz[i].register_symbol_table(symbol_table);
    exp_dvdx[i].register_symbol_table(symbol_table);
    exp_dvdy[i].register_symbol_table(symbol_table);
    exp_dvdz[i].register_symbol_table(symbol_table);
    exp_dwdx[i].register_symbol_table(symbol_table);
    exp_dwdy[i].register_symbol_table(symbol_table);
    exp_dwdz[i].register_symbol_table(symbol_table);

    exprtk::parser<double> parser;

    parser.compile(function_u, exp_u[i]);
    parser.compile(function_v, exp_v[i]);
    parser.compile(function_w, exp_w[i]);
    parser.compile(function_dudx, exp_dudx[i]);
    parser.compile(function_dudy, exp_dudy[i]);
    parser.compile(function_dudz, exp_dudz[i]);
    parser.compile(function_dvdx, exp_dvdx[i]);
    parser.compile(function_dvdy, exp_dvdy[i]);
    parser.compile(function_dvdz, exp_dvdz[i]);
    parser.compile(function_dwdx, exp_dwdx[i]);
    parser.compile(function_dwdy, exp_dwdy[i]);
    parser.compile(function_dwdz, exp_dwdz[i]);
  }
};

flowmap::~flowmap(void) {
	Destruct3D(phi_x_);
	Destruct3D(phi_y_);
	Destruct3D(phi_z_);
	
	Destruct3D(phi_dudx_);
	Destruct3D(phi_dudy_);
	Destruct3D(phi_dudz_);
	Destruct3D(phi_dvdx_);
	Destruct3D(phi_dvdy_);
	Destruct3D(phi_dvdz_);
	Destruct3D(phi_dwdx_);
	Destruct3D(phi_dwdy_);
	Destruct3D(phi_dwdz_);
	
	delete [] exp_u;
	delete [] exp_v;
	delete [] exp_w;
	delete [] exp_dudx;
	delete [] exp_dudy;
	delete [] exp_dudz;
  delete [] exp_dvdx;
  delete [] exp_dvdy;
  delete [] exp_dvdz;
  delete [] exp_dwdx;
  delete [] exp_dwdy;
  delete [] exp_dwdz;
};
