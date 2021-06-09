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
#define how_many_irreps 1

using namespace std;

char alldetails[2000];
char rundetails[2000];
char filename[2000];

double beta, deltabeta;
int Ncol=4;
const int dim=4;
const int nt=4;
const int nspacelike=12;
int max_thread;
int max_lato;
int how_many_confs;
int intero;
double reale;
double Pol_modulus;

int nx;
int ny;
int nz;
ifstream readfile;
ofstream writefile;
ofstream writemodulifile;


int main(int argc, char **argv) {

  deltabeta=0.15;
  how_many_confs=4000;
  nx=nspacelike;
  ny=nspacelike;
  nz=nspacelike;

  int thread_number;

  for (int beta_index=0;beta_index<1;beta_index++) {
    beta=7.42+deltabeta*beta_index;
    sprintf(alldetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf", Ncol, nt, nx, ny, nz, beta);

    sprintf(filename,"datadir/average_Polyakov_%s.dat",alldetails);
    writefile.open(filename);
    writefile << setprecision(precision_digits);
    writefile << fixed;
    sprintf(filename,"datadir/average_Pol_modulus_%s.dat",alldetails);
    writemodulifile.open(filename);
    writemodulifile << setprecision(precision_digits);
    writemodulifile << fixed;
    for (int thread_index=0;thread_index<1;thread_index++) {
      thread_number=thread_index+25*beta_index;
      sprintf(rundetails,"%s_thread_%d", alldetails, thread_number);
      sprintf(filename,"datadir/average_Polyakov_%s.dat",rundetails);
      readfile.open(filename);
cout << "Reading " << filename << endl;
      for (int config=0;config<how_many_confs;config++) {
        readfile >> intero;
	intero+=how_many_confs*thread_index;
        writefile << intero << " ";
        writemodulifile << intero << " ";
        for (int irrep_label=0; irrep_label<how_many_irreps; irrep_label++) {
          readfile >> reale;
          writefile << reale << " ";
	  Pol_modulus=reale*reale;
          readfile >> reale;
          writefile << reale << "  ";
	  Pol_modulus+=reale*reale;
          writemodulifile << sqrt(Pol_modulus) << " ";
        }
        writefile << endl;
        writemodulifile << endl;
      }
      readfile.close();

    }
    writefile.close();
    writemodulifile.close();
  }
  return 0;

}
