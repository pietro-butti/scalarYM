#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <complex>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include "MersenneTwister.h"

time_t rawtime;
struct tm * ptm;
struct tm * timeinfo;

MTRand::uint32 oneSeed = 4357UL;
MTRand::uint32 bigSeed[4] = { 0x130, 0x359, 0x345, 0x587 };


/* Random seed: */
MTRand xx;

/* Fixed seed: */
// MTRand xx( bigSeed, 4 );


#define pi 3.1415926535898
#define twopi 6.28318530717958648
#define maximum_filename_length 2000
#define precision_digits 12

#ifndef __complex_double_structure__
#define __complex_double_structure__

struct complex_double {double re; double im;};

extern"C" {

  int zgeev_(const char *, const char *, int *, double *, int *, complex_double *, complex_double (*)[1], int *,complex_double (*)[1], int *, complex_double*, int *, complex_double *, int *);

}

#endif
