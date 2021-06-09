#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <math.h>
// #include "parameters.h"
#define dim 4
#define Ncol 2
using namespace std;

const int maxthread=68;
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
int number_of_measurements=25;
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
    if (thread<75) {
      rowindex=(thread/5);
      nx=(rowindex%5)+5;
      if (nx==9) {
          nx=10;
      }
      switch (rowindex) {
        case  (0) : { beta=3.4564; break; }
        case  (1) : { beta=3.5408; break; }
        case  (2) : { beta=3.6045; break; }
        case  (3) : { beta=3.6566; break; }
        case  (4) : { beta=3.7425; break; }
        
        case  (5) : { beta=3.1898; break; }
        case  (6) : { beta=3.2751; break; }
        case  (7) : { beta=3.3428; break; }
        case  (8) : { beta=3.4009; break; }
        case  (9) : { beta=3.5000; break; }
        
        case (10) : { beta=2.9568; break; }
        case (11) : { beta=3.0379; break; }
        case (12) : { beta=3.0961; break; }
        case (13) : { beta=3.1564; break; }
        case (14) : { beta=3.2433; break; }
        default   : { break; }
      }
      delta_eta=0.001*(1+(thread%5));
    }
    else {
      nx=10;
      switch (thread) {
        case (75) : { beta=3.7425; delta_eta=0.001; break; }
        case (76) : { beta=3.7425; delta_eta=0.002; break; }
        case (77) : { beta=3.7425; delta_eta=0.003; break; }
        case (78) : { beta=3.7425; delta_eta=0.004; break; }        
        case (79) : { beta=3.5000; delta_eta=0.001; break; }
        case (80) : { beta=3.5000; delta_eta=0.002; break; }
        case (81) : { beta=3.5000; delta_eta=0.003; break; }
        case (82) : { beta=3.5000; delta_eta=0.004; break; }        
        case (83) : { beta=3.2433; delta_eta=0.001; break; }
        case (84) : { beta=3.2433; delta_eta=0.002; break; }
        case (85) : { beta=3.2433; delta_eta=0.003; break; }
        default   : { break; }
      }
    };
    nt=nx+1;
    ny=nx;
#if dim>3
    nz=nx;
#endif
    sprintf(filename,"knl_input_thread_%d.dat",thread);
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
    file << " 0.785398163397448309615660845819875721049" << endl; // This is eta
    file << " " << delta_eta << endl; 
#if Ncol==3
    file << " " << gnu << endl; 
    file << " " << delta_gnu << endl; 
#endif
    file.close();
  }
  delete [] filename;
  
}
