
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
    nt = 3;
    nx = 3;
    ny = 3;
    mode = 0;
    thermalization_time = 500;               
    reunitarization_period = 5;            
    how_many_measurements = 100;
    updates_between_measurements = 1;
    type_of_start = 0;

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
  /*
  for(int mu=0;  mu<dim ; mu++)
  for(int t=0 ; t<nt; t++ )
  for(int x=0 ; x<nx; x++ )
  for(int y=0 ; y<ny; y++ ){
  int site = xx.randInt(nsites);

  dc * test = get_ufield(site, 2);
  for(int i=0; i<4; i++ ) cout << test[i] << endl; //ufield[i + Ncolsquare*get_link(site, 2)] << endl; 
      cout << t << " "<< x <<" "<< y << " " << mu << " "<< get_site(t,x,y) <<  " " << get_link(get_site(t,x,y), mu) << endl;
    }
  */

  //get_plaq_index(5, 2,1);
  cout << get_wilson_action() << endl;
  //cout << plaquette() << endl;



 /*
  ofstream out;
  out.open("prova.dat");

  int Nsweep = 100000;
  int counter = 0;
  for(int tt=0; tt<Nsweep; tt++) {    
    if (Metropolis_sweep_gauge(0.01)==true) counter++;
    out << plaquette() << endl;
    cout << plaquette() << endl;
  }

  out.close();

  cout << float(counter)/float(Nsweep) << endl;





    deallocate_arrays();
*/




}