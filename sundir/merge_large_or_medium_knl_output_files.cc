#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdio>
#include <math.h>
// #include "parameters.h"
#define dim 4
#define Ncol 3
using namespace std;

const int maximum_filename_length=2000;
const bool direct_transformation=true;
const int maxthread=68;
double beta;
int nt, nx, ny;
#if dim>3
int nz;
#endif
#if Ncol==2
double eta=0.785398163397448309615660845819875721049292349843776455243;
#endif
#if Ncol==3
double eta=0.;
#endif
double delta_eta;
#if Ncol==3
double gnu, delta_gnu;
#endif
const int mode=1;//0=thermalization+measurements; 1=only measurements
#if Ncol==2
int thermalization_time=0;//300;
int number_of_measurements;//=400;
#endif
#if Ncol==3
int thermalization_time=200;
int number_of_measurements;//=200;
#endif
int reunitarization_period=20;
const int updates_between_measurements=1;
const int type_of_start=2;//0=hot, 1=cold, 2=read;
const int nsteps=200;
const char Jarzynski_direction='d';
int rowindex;
ofstream input_file;
ofstream output_file;
char *filename, *measfilename;
char *rundetails = new char[maximum_filename_length];

int main() {
    
  filename=new char [2000];
  measfilename=new char [2000];
  for (int rowindex=0; rowindex<
#if Ncol==2
       20;
#endif
#if Ncol==3
       22;
#endif
       rowindex++) {
    
#if Ncol==2
    switch (rowindex%5) {
      case  (0) : { number_of_measurements=400; break; }
      case  (1) : { number_of_measurements=200; break; }
      case  (2) : { number_of_measurements=120; break; }
      case  (3) : { number_of_measurements=50; break; }        
      case  (4) : { number_of_measurements=25; break; }
      default   : { break; }        
    }
#endif
      
// #if Ncol==2
//     sprintf(filename,"%d_medium_knljob.sh",rowindex);
// #endif
// #if Ncol==3
//     sprintf(filename,"%d_large_knljob.sh",rowindex);
// #endif
//     file.open(filename);
//     file << "#!/bin/bash" << endl;
//     file << "#PBS -N " << rowindex;
// #if Ncol==2
//     file << "_medium";
// #endif
// #if Ncol==3
//     file << "_large";
// #endif
//     file << "_knl" << endl;
//     file << "#PBS -l select=1:ncpus=68:mpiprocs=68" << endl;
//     file << "#PBS -A INF17_sft_1" << endl;
//     file << "#PBS -l walltime=23:59:59" << endl;
//     file << "#PBS -M panero@to.infn.it" << endl;
//     file << "#PBS -j oe" << endl;
//     file << "#PBS -m ae" << endl;
//     file << "module purge" << endl;
//     file << "module load env-knl" << endl;
//     file << "module load autoload" << endl;
//     file << "module load lapack/3.6.1--intel--pe-xe-2017--binary" << endl;
//     file << "module load blas/3.6.0--intel--pe-xe-2017--binary" << endl;
//     file << "module load gnu/6.1.0" << endl;
//     file << "module load openmpi/1.10.3-threadmultiple--gnu--6.1.0" << endl;
//     file << "cd ~/knltestdir/" << endl;
// #if Ncol==2
//     file << "mpirun ./medium_knlrun.x " << rowindex << "_medium_knl_input" << endl;
// #endif
// #if Ncol==3
//     file << "mpirun ./large_knlrun.x " << rowindex << "_large_knl_input" << endl;
// #endif
//     file.close();
    cout << "   cat ";
    for (int thread=0; thread<maxthread; thread++) {
      switch (rowindex) {
#if Ncol==2
        case  (0) : { beta=3.4564; nx=10; break; }
        case  (1) : { beta=3.5408; nx=12; break; }
        case  (2) : { beta=3.6045; nx=14; break; }
        case  (3) : { beta=3.6566; nx=16; break; }        
        case  (4) : { beta=3.7425; nx=20; break; }
        
        case  (5) : { beta=3.1898; nx=10; break; }
        case  (6) : { beta=3.2751; nx=12; break; }
        case  (7) : { beta=3.3428; nx=14; break; }        
        case  (8) : { beta=3.4009; nx=16; break; }
        case  (9) : { beta=3.5000; nx=20; break; }
        
        case (10) : { beta=2.9568; nx=10; break; }
        case (11) : { beta=3.0379; nx=12; break; }
        case (12) : { beta=3.0961; nx=14; break; }
        case (13) : { beta=3.1564; nx=16; break; }
        case (14) : { beta=3.2433; nx=20; break; }
        
        case (15) : { beta=2.7124; nx=10; break; }
        case (16) : { beta=2.7938; nx=12; break; }
        case (17) : { beta=2.8598; nx=14; break; }
        case (18) : { beta=2.9115; nx=16; break; }
        case (19) : { beta=3.0071; nx=20; break; }
#endif

#if Ncol==3
        case  (0) : { beta=8.7522; nx=10; break; }
        case  (1) : { beta=8.8997; nx=12; break; }
        case  (2) : { beta=9.0350; nx=14; break; }
        case  (3) : { beta=9.1544; nx=16; break; }
        
        case  (4) : { beta=8.1555; nx=10; break; }        
        case  (5) : { beta=8.3124; nx=12; break; }
        case  (6) : { beta=8.4442; nx=14; break; }
        case  (7) : { beta=8.5598; nx=16; break; }
        
        case  (8) : { beta=7.5687; nx=10; break; }
        case  (9) : { beta=7.7170; nx=12; break; }        
        case (10) : { beta=7.8521; nx=14; break; }
        case (11) : { beta=7.9741; nx=16; break; }
        case (12) : { beta=8.1650; nx=20; break; }
        
        case (13) : { beta=6.9671; nx=10; break; }
        case (14) : { beta=7.1214; nx=12; break; }        
        case (15) : { beta=7.2549; nx=14; break; }
        case (16) : { beta=7.3632; nx=16; break; }
        case (17) : { beta=7.5525; nx=20; break; }
        
        case (18) : { beta=6.5512; nx=9;  break; }
        case (19) : { beta=6.7860; nx=12; break; }        
        case (20) : { beta=6.9748; nx=15; break; }
        case (21) : { beta=7.1190; nx=18; break; }
#endif
        default   : { break; }
      }
      delta_eta=0.0002;
      ny=nx;
#if dim>3
      nz=nx;
#endif
      nt=nx+1;
      sprintf(rundetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_", Ncol, nt, nx, ny);
#if dim==4
      sprintf(rundetails+strlen(rundetails),"nz_%d_", nz);
#endif
      sprintf(rundetails+strlen(rundetails),"beta_%12.10lf_", beta);
      sprintf(rundetails+strlen(rundetails),"nsteps_%d_%s_eta_%12.10lf_deltaeta_%12.10lf_", nsteps, (direct_transformation)?"direct":"reverse", eta, delta_eta);
#if Ncol==3
      sprintf(rundetails+strlen(rundetails),"gnu_%12.10lf_deltagnu_%12.10lf_", gnu, delta_gnu);
#endif
      sprintf(rundetails+strlen(rundetails),"thread_%d", thread);
  
      sprintf(measfilename,"datadir/Jarzynski_SF_%s.dat",rundetails);
      
      cout << " " << measfilename;
//       measfile.open(measfilename,fstream::in);

//       measfile.close();
    }
    cout << " > ";
    sprintf(rundetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_", Ncol, nt, nx, ny);
#if dim==4
    sprintf(rundetails+strlen(rundetails),"nz_%d_", nz);
#endif
    sprintf(rundetails+strlen(rundetails),"beta_%12.10lf_", beta);
    sprintf(rundetails+strlen(rundetails),"nsteps_%d_%s_eta_%12.10lf_deltaeta_%12.10lf_", nsteps, (direct_transformation)?"direct":"reverse", eta, delta_eta);
#if Ncol==3
    sprintf(rundetails+strlen(rundetails),"gnu_%12.10lf_deltagnu_%12.10lf", gnu, delta_gnu);
#endif
    sprintf(measfilename,"datadir/Jarzynski_SF_%s.dat",rundetails);
      
    cout << measfilename << endl;
  }
  delete [] filename;
  
}
