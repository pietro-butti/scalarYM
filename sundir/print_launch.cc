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
#include "variables.h"

using namespace std;

ofstream outputfile;
char *filename=new char[filename_length];

int main() {

  for (int beta_index=0; beta_index<how_many_betas; beta_index++)
  for (int dist=mindist; dist<=maxdist; dist++) {
    sprintf(filename,"%sdata_beta%d_dir/dist_%d_dir/launch.cmd",basicpath,betas[beta_index],dist);
    outputfile.open(filename);
/////////////////////////////////////////////////////////// 
// For vuori:
    outputfile << "#!/bin/csh" << endl;
///////////////////////////////////////////////////////////
// For taito:
//     outputfile << "#!/bin/bash" << endl;
    outputfile << "#SBATCH -J b_"<< betas[beta_index] << "_dist_" << dist << endl;
    outputfile << "#SBATCH -e err_%j" << endl;
    outputfile << "#SBATCH -o output_%j" << endl;
    outputfile << "#SBATCH --mail-type=END" << endl;
    outputfile << "#SBATCH --mail-user=marco.panero@helsinki.fi " << endl;
    outputfile << "#SBATCH --mem-per-cpu=" << mem_per_cpu << endl;
    outputfile << "#SBATCH -t "<< time_limit << endl;
    outputfile << "##the number of processes (number of cores)" << endl;
    outputfile << "#SBATCH -n 1" << endl;
    outputfile << "#SBATCH -p " << queue_name << endl;
    outputfile << "#" << endl;
/////////////////////////////////////////////////////////// 
// For vuori:
   outputfile << "module swap PrgEnv-pgi PrgEnv-gnu" << endl;
   outputfile << "srun /wrk/panero/bin/su3_ahiggs_multilevel" << endl;
///////////////////////////////////////////////////////////
// For taito:
//    outputfile << "module load PrgEnv-gnu" << endl;
//    outputfile << "module load openmpi/1.6-gcc" << endl;
//    outputfile << "module load openblas/0.2.6" << endl;
//    outputfile << "srun su3_ahiggs_multilevel" << endl;
    outputfile.close();
    printf("cd %sdata_beta%d_dir/dist_%d_dir/ ; sbatch launch.cmd\n",basicpath,betas[beta_index],dist);
  }
  
  delete filename;
  return 0;

}
