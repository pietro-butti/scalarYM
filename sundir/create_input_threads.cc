#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <math.h>
#include "parameters.h"
#define dim 4
#define Ncol 3
using namespace std;

const int maxthread=48;
double beta;
int nt, nx, ny;
#if dim>3
int nz;
#endif
double eta, delta_eta;
#if Ncol==3
double gnu, delta_gnu;
#endif
const int mode=0;//0=thermalization+measurements; 1=only measurements
int thermalization_time=200;
int reunitarization_period=100;
int number_of_measurements;
const int updates_between_measurements=1;
const int type_of_start=0;//0=hot, 1=cold, 2=read;
const int nsteps=200;
const char Jarzynski_direction='d';
int rowindex;
ofstream file;
char *filename;

int main() {
    
  filename=new char [2000];
  for (int thread=0; thread<maxthread; thread++) {
    rowindex=(thread/3);
    switch (rowindex) {
      case  (0) : { beta=8.7522; nx=10;  break; }
      case  (1) : { beta=8.8997; nx=12;  break; }
      case  (2) : { beta=9.035;  nx=14;  break; }
      case  (3) : { beta=9.1544; nx=16;  break; }
      
      case  (4) : { beta=8.1555; nx=10;  break; }
      case  (5) : { beta=8.3124; nx=12;  break; }
      case  (6) : { beta=8.4442; nx=14;  break; }
      case  (7) : { beta=8.5598; nx=16;  break; }
      
      case  (8) : { beta=7.5687; nx=10;  break; }
      case  (9) : { beta=7.717;  nx=12;  break; }
      case (10) : { beta=7.8521; nx=14;  break; }
      case (11) : { beta=7.9741; nx=16;  break; }
      case (12) : { beta=8.165;  nx=20; break; }
      
      case (13) : { beta=6.9671; nx=10;  break; }
      case (14) : { beta=7.1214; nx=12;  break; }
      case (15) : { beta=7.2549; nx=14;  break; }
      case (16) : { beta=7.3632; nx=16;  break; }
      case (17) : { beta=7.5525; nx=20; break; }
      
      case (18) : { beta=6.5512; nx=9;  break; }
      case (19) : { beta=6.786;  nx=12;  break; }
      case (20) : { beta=6.9748; nx=15; break; }
      case (21) : { beta=7.119;  nx=18; break; }
      
      default   : { break; }
    }
    switch (thread%3) {
      case  (0) : { delta_eta=0.001; break; }
      case  (1) : { delta_eta=0.0005; break; }
      case  (2) : { delta_eta=0.0002; break; }
      default   : { break; }
    }
    switch (nx) {
      case  (5)  : { number_of_measurements=400; break; }
      case  (6)  : { number_of_measurements=200; break; }
      case  (7)  : { number_of_measurements=100; break; }
      case  (8)  : { number_of_measurements=50; break; }
      case  (10) : { number_of_measurements=20; break; }
      case  (12) : { number_of_measurements=10; break; }
      default   : { break; }
    }
    nt=nx+1;
    ny=nx;
#if dim>3
    nz=nx;
#endif
    sprintf(filename,"other_SU3_input_thread_%d.dat",thread);
    file.open(filename);
    file << setprecision(12);
    file << fixed;
    file << " " << beta << endl;
    file << " " << nt << endl;
    file << " " << nx << endl;
    file << " " << ny << endl;
#if dim>3
    file << " " << nz << endl;
#endif
    file << " " << mode << endl; 
    file << " " << thermalization_time << endl; 
    file << " " << reunitarization_period << endl; 
    file << " " << number_of_measurements << endl; 
    file << " " << updates_between_measurements << endl; 
    file << " " << type_of_start << endl;
    file << " " << nsteps << endl;
    file << " " << Jarzynski_direction << endl;
    file << " 0.0" << endl; // This is eta
    file << " " << delta_eta << endl; 
#if Ncol==3
    file << " 0.0" << endl; // This is gnu
    file << " 0.0" << endl; // This is delta_gnu
#endif
    file.close();
  }
  delete [] filename;
  
}
