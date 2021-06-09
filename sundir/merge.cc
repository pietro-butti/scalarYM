#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <ctime>
#include <complex>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#define precision_digits 12


using namespace std;

char alldetails[2000];
char rundetails[2000];
char filename[2000];

double beta;
int Ncol;
int dim;
int nt;
int nspacelike;
int max_thread;
int max_lato;
int how_many_confs;
int intero;
double reale;

int nx;
int ny;
int nz;
ifstream readfile;
ofstream writefile;

int main(int argc, char **argv) {
  
  if (argc!=9) {
    cout << "Usage: ./merge.x <beta> <Ncol> <dim> <nt> <nspacelike> <max_thread> <max_lato> <how_many_confs>" << endl; 
    exit(0);
  };

  beta=(double) atof(argv[1]);
  Ncol=atoi(argv[2]);
  dim=atoi(argv[3]);
  nt=atoi(argv[4]);
  nspacelike=atoi(argv[5]);
  max_thread=atoi(argv[6]);
  max_lato=atoi(argv[7]);
  how_many_confs=atoi(argv[8]);
  
  if ((dim!=3)&&(dim!=4)) {
    cout << "The number of spacetime dimensions must be either 3 or 4" << endl;
  }
  nx=nspacelike;
  ny=nspacelike;
  nz=nspacelike;

  if (dim==3) {
    sprintf(alldetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf", Ncol, nt, nx, ny, beta);
  }
  else if (dim==4) {
    sprintf(alldetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf", Ncol, nt, nx, ny, nz, beta);
  };
  
  sprintf(filename,"datadir/smeared_Wilson_%s.dat",alldetails);
  writefile.open(filename);
  writefile << setprecision(precision_digits);
  writefile << fixed;
  for (int thread_number=0;thread_number<max_thread;thread_number++) {

    sprintf(rundetails,"%s_thread_%d", alldetails, thread_number);
    sprintf(filename,"datadir/smeared_Wilson_%s.dat",rundetails);
    readfile.open(filename);
    for (int config=0;config<how_many_confs;config++)
    for (int r=1;r<=max_lato;r++)
    for (int L=1;L<=max_lato;L++) {
      readfile >> intero;
      writefile << intero+how_many_confs*thread_number << " ";
      readfile >> intero;
      writefile << intero << " ";
      readfile >> intero;
      writefile << intero << " ";
      readfile >> reale;
      writefile << reale << " ";
      readfile >> reale;
      writefile << reale << endl;
    }
    readfile.close();
  
  }
  writefile.close();
  return 0;

}
