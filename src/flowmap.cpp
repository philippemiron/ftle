#include "flowmap.h"
extern "C" {
    void dcopy_(int& n, double* x, int& incx, double* y, int& incy);
    void daxpy_(int& n, double& alpha, double* x, int& incx, double* y, int& incy);  
}

flowmap::flowmap(parameter* objpara)
	:     
    Nx(objpara->get_Nx()),
    Ny(objpara->get_Ny()),
    Nz(objpara->get_Nz()),
    xmin(objpara->getXmin()),
    xmax(objpara->getXmax()),
    ymin(objpara->getYmin()),
    ymax(objpara->getYmax()),
    zmin(objpara->getZmin()),
    zmax(objpara->getZmax()),
    T(objpara->get_T()),
    T0(objpara->get_T0()),
    Npt(objpara->get_Npt()),
    ode(objpara->get_MthOde()),
    function_u(objpara->getFunctionU()),
    function_v(objpara->getFunctionV()),
    function_w(objpara->getFunctionW()),
    function_dudx(objpara->getFunctionDudx()),
    function_dudy(objpara->getFunctionDudy()),
    function_dudz(objpara->getFunctionDudz()),
    function_dvdx(objpara->getFunctionDvdx()),
    function_dvdy(objpara->getFunctionDvdy()),
    function_dvdz(objpara->getFunctionDvdz()),
    function_dwdx(objpara->getFunctionDwdx()),
    function_dwdy(objpara->getFunctionDwdy()),
    function_dwdz(objpara->getFunctionDwdz())	        
{

cout << "Flow map integrations (T=" << T << "s) for the " << Nx*Ny*Nz << " particles." << endl;

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

// choix du right hand side choisi lors du calcul des déplacements
cout << "	Using analytical equation" << endl;	

// Setting up the expression evaluation

// Create a function with theses expressions
int number_of_threads = omp_get_max_threads();
initialize_expression(number_of_threads);

// Formation des matrices de deplacements.
x.resize(Nx, 0.0);
y.resize(Ny, 0.0);
z.resize(Nz, 0.0);

for (int i(0); i<Nx; i++)
    x[i] = xmin + i*(xmax-xmin)/(Nx-1);

for (int i(0); i<Ny; i++)
    y[i] = ymin + i*(ymax-ymin)/(Ny-1);
    
for (int i(0); i<Nz; i++)
    z[i] = zmin + i*(zmax-zmin)/(Nz-1);

phi_x = Construct3D(Nx, Ny, Nz);
phi_y = Construct3D(Nx, Ny, Nz);
phi_z = Construct3D(Nx, Ny, Nz);

phi_dudx = Construct3D(Nx, Ny, Nz);
phi_dudy = Construct3D(Nx, Ny, Nz);
phi_dudz = Construct3D(Nx, Ny, Nz);
phi_dvdx = Construct3D(Nx, Ny, Nz);
phi_dvdy = Construct3D(Nx, Ny, Nz);
phi_dvdz = Construct3D(Nx, Ny, Nz);
phi_dwdx = Construct3D(Nx, Ny, Nz);
phi_dwdy = Construct3D(Nx, Ny, Nz);
phi_dwdz = Construct3D(Nx, Ny, Nz);


// Pour déterminer combien de particules sortent du domaine
// lors du calcul des déplacements
                
// 1- Faire une boucle sur les noeuds du maillage
// 2- Trouver les coordonnees du noeud courant
// 3- Calculer le delacement Lagrangien du noeud  
//    a l'aide des vitesses mesurees par PIV
#pragma omp parallel for
for (int i=0; i<Nx;i++) {
    for (int j=0; j<Ny; j++) { 
        for (int k=0; k<Nz; k++) {
            // deplacement du point courant
            vector<double> d(30,0);
            // Recupere les coordonnées du noeuds courants
            d[0]=x[i];
            d[1]=y[j];
            d[2]=z[k];                

            // Les paramètres de l'intégration du déplacement sont:
            // Nous intégrons de T0 a T0+T     // T: Temps d'integration des particules
            // Npt: Nombre de points (subdivisions) lors de l'integration
            // T0: Temps t0 de notre intégration
            // d: vecteurs de déplacements
            (this->*ptr_ode) (T, Npt, T0, d);

            phi_x[i][j][k] = d[0];
            phi_y[i][j][k] = d[1];
            phi_z[i][j][k] = d[2];

            phi_dudx[i][j][k] = d[3];
            phi_dudy[i][j][k] = d[4];
            phi_dudz[i][j][k] = d[5];
            phi_dvdx[i][j][k] = d[6];
            phi_dvdy[i][j][k] = d[7];
            phi_dvdz[i][j][k] = d[8];
            phi_dwdx[i][j][k] = d[9];
            phi_dwdy[i][j][k] = d[10];
            phi_dwdz[i][j][k] = d[11];
        }
    }
}

}

flowmap::~flowmap(void) {
	Destruct3D(phi_x);
	Destruct3D(phi_y);
	Destruct3D(phi_z);
	
	Destruct3D(phi_dudx);
	Destruct3D(phi_dudy);
	Destruct3D(phi_dudz);
	Destruct3D(phi_dvdx);
	Destruct3D(phi_dvdy);
	Destruct3D(phi_dvdz);
	Destruct3D(phi_dwdx);
	Destruct3D(phi_dwdy);
	Destruct3D(phi_dwdz);
	
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
	
}

void flowmap::int_euler (double T, int Npt, double T0, vector<double>& D) {

int  nd = 12;
int incx = 1;
double u[nd];
double x0[nd];

// ==========
x0[0]=D[0];   // x_0
x0[1]=D[1];   // y_0
x0[2]=D[2];   // z_0
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
// Si T est négatif dT sera aussi négatif
double dt=T/(Npt-1);

for (int k(0); k<Npt; k++) {
    
    // velocity interpolation
    double t  = T0 + k*dt;
    
    interpolation(t,x0,u);
    //cout << "x=" << x0[0] << " y=" << x0[1] << endl;
    //cout << "u=" << g[0] << " v=" << g[1] << endl;

    // x0 += u*dt
	daxpy_(nd, dt, u, incx, x0, incx);
}
dcopy_(nd, x0, incx, &D[0], incx);
};

void flowmap::int_rk5 (double T, int Npt, double T0, vector<double>& D)
{
int  nd = 30;
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
y1[0]=D[0];   // x_0
y1[1]=D[1];   // y_0
y1[2]=D[2];   // z_0
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

// Si T est négatif, dT sera aussi négatif et nous intégrons de T0 à T qui se situe avant T0
double dT=T/(Npt-1);

for(int k=0;k<Npt-1;k++) {
	// Détermine les différents temps nécessaire pour l'intégration 
	double t  = T0 + k*dT;
	double t1 = t + dT*0.25;
	double t2 = t + dT*0.5 ;
	double t3 = t + dT*0.75;
	double t4 = t + dT*1.0;
	//cout << "k: " << k << " t: " << t << " t1: " << t1 << " t2: " << t2 << " t3: " << t3 << " t4: " << t4 << endl;
	// ==============================================
	interpolation(t,y1,g); 
	for(int i=0;i<nd;i++) {
		 k1[i]=g[i]*dT;
		 y2[i]=y1[i]+0.25*k1[i];
	}
	// ==============================================
	interpolation(t1,y2,g); 
	for(int i=0;i<nd;i++) {
		 k2[i]=g[i]*dT;
		 y3[i]=y1[i]+(k1[i]+k2[i])/8.0;
	}
	// ==============================================
	interpolation(t1,y3,g); 
	for(int i=0;i<nd;i++) {
		 k3[i]=g[i]*dT;
		 y4[i]=y1[i]-k2[i]/2.0+k3[i];
	}
	// ==============================================
	interpolation(t2,y4,g); 
	for(int i=0;i<nd;i++) {
		 k4[i]=g[i]*dT;
		 y5[i]=y1[i]+(3.0*k1[i]+9.0*k4[i])/16.0;
	}
	// ==============================================
	interpolation(t3,y5,g); 
	for(int i=0;i<nd;i++) {
		 k5[i]=g[i]*dT;
		 y6[i]=y1[i]+(-3.0*k1[i]+2.0*k2[i]+12.0*k3[i]-12.0*k4[i]+8.0*k5[i])/7.0;
	}
	// ==============================================
	interpolation(t4,y6,g); 
	for(int i=0;i<nd;i++) {
		 k6[i]=g[i]*dT;
		 y6[i]=y1[i]+(7.0*k1[i]+32.0*k3[i]+12.0*k4[i]+32.0*k5[i]+7.0*k6[i])/90.0;
	}
	for(int i=0;i<nd;i++) y1[i]=y6[i];
}
for(int i=0;i<nd;i++) D[i]=y1[i];
}

void flowmap::interpolation (double t, double* y, double* g)
{
    //cout << "x:" << y[0] << " y: " << y[1] << endl;
    double f[12];
    
    // TODO
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
    g[0]=f[0];
    g[1]=f[1];
    g[2]=f[2];

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

    //cout << "phi:" << f[0] << " " << f[1] << " " << f[2] << " " << f[3] << " " << f[4] << " " << f[5] << endl;
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
        symbol_table.add_variable("x",var_x[i]);
        symbol_table.add_variable("y",var_y[i]);
        symbol_table.add_variable("z",var_z[i]);
        symbol_table.add_variable("t",var_t[i]);
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
}

// Méthode pour récupérer les données de la classe flowmap
double*** flowmap::get_phi_x() {
	return phi_x;
}

double*** flowmap::get_phi_y() {
	return phi_y;
}

double*** flowmap::get_phi_z() {
	return phi_z;
}

double*** flowmap::get_phi_dudx() {
	return phi_dudx;
}

double*** flowmap::get_phi_dudy() {
	return phi_dudy;
}

double*** flowmap::get_phi_dudz() {
	return phi_dudz;
}

double*** flowmap::get_phi_dvdx() {
	return phi_dvdx;
}

double*** flowmap::get_phi_dvdy() {
	return phi_dvdy;
}

double*** flowmap::get_phi_dvdz() {
	return phi_dvdz;
}

double*** flowmap::get_phi_dwdx() {
	return phi_dwdx;
}

double*** flowmap::get_phi_dwdy() {
	return phi_dwdy;
}

double*** flowmap::get_phi_dwdz() {
	return phi_dwdz;
}

std::vector<double> flowmap::get_x() {
    return x;
}

std::vector<double> flowmap::get_y() {
    return y;
}

std::vector<double> flowmap::get_z() {
    return z;
}
