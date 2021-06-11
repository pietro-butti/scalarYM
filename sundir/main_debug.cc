
#include "headers.h"
#include "parameters.h"

int counter=0;

#include "more_headers.h"
#include "included_routines.h"


int main(int argc, char **argv) {


  time(&rawtime);
  ptm=gmtime(&rawtime);
  timeinfo=localtime(&rawtime);

  cout << setprecision(precision_digits);
  cout << fixed;

    // PARAMETERS  -------------------------
    beta = 6.0;
    nt = 10;
    nx = 10;
    ny = 10;
    mode = 0;
    thermalization_time = 500;               
    reunitarization_period = 5;            
    how_many_measurements = 100;
    updates_between_measurements = 1;
    type_of_start = 1;

    // derived variables
    nsites=nt*nx*ny;
    spatial_volume=nx*ny;
    nlinks=dim*nsites;
    ufielddimension=nlinks*Ncolsquare;
    ufield = new dc[ufielddimension];
    neighbor_plus = new int[nlinks];
    neighbor_minus = new int[nlinks];




    // =========================== PRODUCE INTERFACE       ============================ //
// ================================================================================ //
    
    output_stream << "# ======================================================================\n";
    output_stream << "# Local time and date at run start: " << asctime(timeinfo);
    output_stream << "# ----------------------------------------------------------------------\n";
    output_stream << "# Lattice setup:                 "          << endl;
    output_stream << "#                         Ncol = " << Ncol  << endl;
    output_stream << "#                          dim = " << dim   << endl; 
    output_stream << "#                                "          << endl;
    output_stream << "#                         beta = " << beta  << endl;
    output_stream << "#                           nt = " << nt    << endl;
    output_stream << "#                           nx = " << nx << endl;
    output_stream << "#                           ny = " << ny << endl;
    output_stream << "# ----------------------------------------------------------------------\n";
    output_stream << "#                         mode = " << mode << endl;
    output_stream << "#          thermalization time = " << thermalization_time << endl;
    output_stream << "#       reunitarization period = " << reunitarization_period << endl;
    output_stream << "#       number of measurements = " << how_many_measurements << endl;
    output_stream << "# updates between measurements = " << updates_between_measurements << endl;
    output_stream << "#                   start type = " << type_of_start << " ";
        output_stream << "(random start)";

    output_stream << endl;
    output_stream << "# ======================================================================\n";
// ================================================================================ //
// ================================================================================ //

  initialize();
  

  ofstream out;
  out.open("prova.dat");

  int Nsweep = 1000;
  int counter = 0;
  for(int tt=0; tt<Nsweep; tt++) {    
    if (Metropolis_sweep_gauge(.01)==true) counter++;
    out << plaquette() << endl;
    cout << plaquette() << endl;
  }

  out.close();

  cout << float(counter)/float(Nsweep) << endl;





    deallocate_arrays();





}