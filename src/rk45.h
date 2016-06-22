#ifndef LCS_RK45_H
#define LCS_RK45_H

// rk45 coefficients
const double a21 = 1.0 / 5.0;
const double a31 = 3.0 / 40.0;
const double a32 = 9.0 / 40.0;
const double a41 = 44.0 / 45.0;
const double a42 = -56.0 / 15.0;
const double a43 = 32.0 / 9.0;
const double a51 = 19372.0 / 6561.0;
const double a52 = -25360.0 / 2187.0;
const double a53 = 64448.0 / 6561.0;
const double a54 = -212.0 / 729.0;
const double a61 = 9017.0 / 3168.0;
const double a62 = -355.0 / 33.0;
const double a63 = 46732.0 / 5247.0;
const double a64 = 49.0 / 176.0;
const double a65 = -5103.0 / 18656.0;
const double a71 = 35.0 / 384.0;
//const double a72 = 0.0
const double a73 = 500.0 / 1113.0;
const double a74 = 125.0 / 192.0;
const double a75 = -2187.0 / 6784.0;
const double a76 = 11.0 / 84.0;

//const double c1 = 0.0
const double c2 = 1.0 / 5.0;
const double c3 = 3.0 / 10.0;
const double c4 = 4.0 / 5.0;
const double c5 = 8.0 / 9.0;
//const double c6 = 1.0;
//const double c7 = 1.0;

// bi = a7i
//const double b1 = 35.0 / 384.0;
//const double b2 = 0.0;
//const double b3 = 500.0 / 1113.0;
//const double b4 = 125.0 / 192.0;
//const double b5 = -2187.0 / 6784.0;
//const double b6 = 11.0 / 84.0;
//const double b7 = 0.0

//const double b1p = 5179.0 / 57600.0;
//const double b2p = 0.0
//const double b3p = 7571.0 / 16695.0;
//const double b4p = 393.0 / 640.0;
//const double b5p = -92097.0 / 339200.0;
//const double b6p = 187.0 / 2100.0;
//const double b7p = 1.0 / 40.0;

// ei = bi - bip
const double e1=71.0/57600.0;
const double e3=-71.0/16695.0;
const double e4=71.0/1920.0;
const double e5=-17253.0/339200.0;
const double e6=22.0/525.0;
const double e7=-1.0/40.0;

#endif //LCS_RK45_H