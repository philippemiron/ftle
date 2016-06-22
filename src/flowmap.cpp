#include "flowmap.h"

flowmap::flowmap(const shared_ptr<parameter> &objpara) :
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
    t0(objpara->t0()),
    tol(objpara->tolerance()),
    hmin(objpara->hmin()),
    hmax(objpara->hmax()) {
  // Setting up the expression evaluation
  initialize_expression(omp_get_max_threads());

  x_.resize(nx, 0.0);
  y_.resize(ny, 0.0);
  z_.resize(nz, 0.0);

  for (int i(0); i < nx; i++)
    x_[i] = xmin + i * (xmax - xmin) / (nx - 1);
  for (int i(0); i < ny; i++)
    y_[i] = ymin + i * (ymax - ymin) / (ny - 1);
  for (int i(0); i < nz; i++)
    z_[i] = zmin + i * (zmax - zmin) / (nz - 1);

  vecResize(phi_x_, nx, ny, nz);
  vecResize(phi_y_, nx, ny, nz);
  vecResize(phi_z_, nx, ny, nz);
  vecResize(phi_dudx_, nx, ny, nz);
  vecResize(phi_dudy_, nx, ny, nz);
  vecResize(phi_dudz_, nx, ny, nz);
  vecResize(phi_dvdx_, nx, ny, nz);
  vecResize(phi_dvdy_, nx, ny, nz);
  vecResize(phi_dvdz_, nx, ny, nz);
  vecResize(phi_dwdx_, nx, ny, nz);
  vecResize(phi_dwdy_, nx, ny, nz);
  vecResize(phi_dwdz_, nx, ny, nz);
};

void flowmap::calculate_trajectories() {
  // Calculate trajectory of every particle
  // at every node of the mesh
#pragma omp parallel for
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        vector<double> d(12, 0);
        // initial position
        d[0] = x_[i];
        d[1] = y_[j];
        d[2] = z_[k];

        // solve trajectories for each particle and retrieve final position
        rk45(t0, t0 + t, tol, hmin, hmax, 100000, d);

        // save values
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
};

void flowmap::rk45(double t0, double tend, double tol, double hmin, double hmax, size_t maxiter, vector<double> &d) {

  vector<double> x(12, 0.0);
  // x_0, y_0, z_0
  x[0] = d[0];
  x[1] = d[1];
  x[2] = d[2];
  // Derivative initial values
  x[3] = 1.0;
  x[4] = 0.0;
  x[5] = 0.0;
  x[6] = 0.0;
  x[7] = 1.0;
  x[8] = 0.0;
  x[9] = 0.0;
  x[10] = 0.0;
  x[11] = 1.0;
  double t(t0);
  double h(1.0e-02 * copysign(1.0, tend - t0)); // initial time step
  vector<double> K1, K2, K3, K4, K5, K6, K7;
  vector<double> x1(x.size()), x2(x.size()), x3(x.size()), x4(x.size()), x5(x.size()), x6(x.size());
  vector<double> error(x.size(), 0.0);
  double delta(1.0), err(1.0), error_max(1.0);

  for (size_t i(0); i < maxiter; i++) {
    K1 = velocity(t, x);
    for (size_t j(0); j < x.size(); j++)
      x1[j] = x[j] + h * (a21 * K1[j]);

    K2 = velocity(t + c2 * h, x1);
    for (size_t j(0); j < x.size(); j++)
      x2[j] = x[j] + h * (a31 * K1[j] + a32 * K2[j]);

    K3 = velocity(t + c3 * h, x2);
    for (size_t j(0); j < x.size(); j++)
      x3[j] = x[j] + h * (a41 * K1[j] + a42 * K2[j] + a43 * K3[j]);

    K4 = velocity(t + c4 * h, x3);
    for (size_t j(0); j < x.size(); j++)
      x4[j] = x[j] + h * (a51 * K1[j] + a52 * K2[j] + a53 * K3[j] + a54 * K4[j]);

    K5 = velocity(t + c5 * h, x4);
    for (size_t j(0); j < x.size(); j++)
      x5[j] = x[j] + h * (a61 * K1[j] + a62 * K2[j] + a63 * K3[j] + a64 * K4[j] + a65 * K5[j]);

    K6 = velocity(t + h, x5);
    for (size_t j(0); j < x.size(); j++)
      x6[j] = x[j] + h * (a71 * K1[j] + a73 * K3[j] + a74 * K4[j] + a75 * K5[j] + a76 * K6[j]);

    K7 = velocity(t + h, x6);
    for (size_t j(0); j < x.size(); j++)
      error[j] = h * (e1 * K1[j] + e3 * K3[j] + e4 * K4[j] + e5 * K5[j] + e6 * K6[j] + e7 * K7[j]);

    // error control
    err = norm(error);
    error_max = max(tol, tol * norm(x));
    if (err < error_max) {
      // accept time step
      t += h;
      x = x6;
    }

    // step variation
    if (err > 0) {
      delta = 0.9 * pow(error_max / err, 0.2);
      delta = min(max(delta, 0.1), 5.0);
    }
    h *= delta;
    h = min(fabs(h), hmax) * copysign(1.0, h);

    // adjust last step
    if (t == tend) {
      break;
    } else if (h > 0.0 and t + h > tend) {
      h = tend - t;
    } else if (h < 0.0 and t + h < tend) {
      h = tend - t;
    } else if (abs(h) < hmin) {
      break;
    }

    if (i == maxiter - 1)
      cerr << "reach max iterations" << endl;
  }
  d = x;
}

vector<double> flowmap::velocity(double t, vector<double> &y) {
  vector<double> f(12, 0.0);
  vector<double> g(12, 0.0);

  // one variable and one expression
  // for every thread
  int tid = omp_get_thread_num();
  var_x[tid] = y[0];
  var_y[tid] = y[1];
  var_z[tid] = y[2];
  var_t[tid] = t;

  // Analytical solution of velocity components and derivatives
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
  f[9] = exp_dwdx[tid].value();
  f[10] = exp_dwdy[tid].value();
  f[11] = exp_dwdz[tid].value();

  // vector velocity
  g[0] = f[0];
  g[1] = f[1];
  g[2] = f[2];

  // and first derivatives
  g[3] = f[3] * y[3] + f[4] * y[6] + f[5] * y[9];
  g[4] = f[3] * y[4] + f[4] * y[7] + f[5] * y[10];
  g[5] = f[3] * y[5] + f[4] * y[8] + f[5] * y[11];
  g[6] = f[6] * y[3] + f[7] * y[6] + f[8] * y[9];
  g[7] = f[6] * y[4] + f[7] * y[7] + f[8] * y[10];
  g[8] = f[6] * y[5] + f[7] * y[8] + f[8] * y[11];
  g[9] = f[9] * y[3] + f[10] * y[6] + f[11] * y[9];
  g[10] = f[9] * y[4] + f[10] * y[7] + f[11] * y[10];
  g[11] = f[9] * y[5] + f[10] * y[8] + f[11] * y[11];

  return g;
}

void flowmap::initialize_expression(int number_of_threads) {

  // one expression per threads to avoid race condition
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
  for (int i(0); i < number_of_threads; i++) {
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
  delete[] exp_u;
  delete[] exp_v;
  delete[] exp_w;
  delete[] exp_dudx;
  delete[] exp_dudy;
  delete[] exp_dudz;
  delete[] exp_dvdx;
  delete[] exp_dvdy;
  delete[] exp_dvdz;
  delete[] exp_dwdx;
  delete[] exp_dwdy;
  delete[] exp_dwdz;
};
