#ifndef _ftle_
#define _ftle_

// LES .h DU SYSTEME
#include <stdio.h>
#include <string>
#include <string.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <tecio.h>
#include <mkl.h>
#include <limits>
#include "exprtk.hpp"

// LES MACROS LES PLUS UTILES
#ifndef signe
#define signe(a)   ( ((a)>=0.0) ? (1.0) : (-1.0) )
#endif

#ifndef max
#define max(a,b)   ( ((a)>(b)) ? (a) : (b) )
#endif

#ifndef min
#define min(a,b)   ( ((a)>(b)) ? (b) : (a) )
#endif

// LES VARIABLES LOGIQUES
#define VRAI       (1)
#define FAUX       (0)

#define INT32      (1)

#endif

