/* PROGRAM sun.cc */
/**/
/* Simulation of Yang-Mills theories */
/* with gauge group SU(N) (for 2<=N<=10), */
/* with or without double-trace deformation, */
/* in D=3 and D=4 spacetime dimensions. */
/* This program can measure: */
/* - average, timelike and spacelike plaquette in each configuration */
/* - average value of the Polyakov loop in each configuration, */
/*   for the first 12 independent irreducible representations */
/* - zero-transverse-momentum Polyakov loops */
/* - Wilson loops from unsmeared links (either averaging over all */
/*   planes, or over spacelike ones only) */
/* - Wilson loops from smeared links (averaging over timelike planes) */
/* - correlators of Polyakov loops, using the multi-level algorithm */
/* - Schroedinger functional, using Jarzinski's theorem */
/* Convention for directions: 0=t, 1=x, 2=y(, 3=z) */

#include "headers.h"
#include "parameters.h"

int counter=0;

#include "more_headers.h"
#include "included_routines.h"


int main(int argc, char **argv) {



// ============================ READING ARGUMENTS ================================= //
// ================================================================================ //
  if (argc>1) {
    sprintf(input_file_string,"%s",argv[1]);
  }

#ifdef __parallel__
  if (argc!=2) {
    cout << "Usage for parallel mode: ./run.x <input_file_string>" << endl; 
    exit(0);
  };
  MPI::Init(argc,argv);
  p=MPI::COMM_WORLD.Get_size();
  id=MPI::COMM_WORLD.Get_rank();
  thread_number=id+parallel_thread_offset;
  for (int seed_index=0;seed_index<4;seed_index++) {
    bigSeed[seed_index]+=xx.randInt(104729)+(seed_index+1)*thread_number;
  }
  xx.seed( bigSeed, 4 );
#else

  if (argc!=3) {
    cout << "Usage for scalar mode: ./scalar <input_file_string> <thread_number>" << endl; 
    exit(0);
  };
  thread_number=atoi(argv[2]);

#endif
// ================================================================================ //
// ================================================================================ //

  time(&rawtime);
  ptm=gmtime(&rawtime);
  timeinfo=localtime(&rawtime);

  cout << setprecision(precision_digits);
  cout << fixed;














// ============================= READ PARAMETERS and allocate (some memory) =================================== //
// ================================================================================ //
    char *inputfilename = new char[maximum_filename_length];

    sprintf(inputfilename,"%s_thread_%d.dat",input_file_string,thread_number);

    cout << inputfilename << endl;

    ifstream inputfile;
    inputfile.open(inputfilename);
    if (!inputfile.is_open()) {
      cout << "Error! " << inputfilename << " does not exists" << endl;
      cout << "Execution aborted\n";
    }
    inputfile >> beta;
    inputfile >> nt;
    inputfile >> nx;
    inputfile >> ny;
    nsites=nt*nx*ny;
    spatial_volume=nx*ny;
  #if dim>3
    inputfile >> nz;
    nsites*=nz;
    spatial_volume*=nz;
  #endif
    nlinks=dim*nsites;
    neighbor_plus = new int[nlinks];
    neighbor_minus = new int[nlinks];



    ufielddimension=nlinks*Ncolsquare;
    ufield = new dc[ufielddimension];

    
    inputfile >> mode;
    inputfile >> thermalization_time;               
    inputfile >> reunitarization_period;            
    inputfile >> how_many_measurements;
    inputfile >> updates_between_measurements;
    inputfile >> type_of_start;

    inputfile.close();
// ================================================================================ //
// ================================================================================ //






















// =========================== PRODUCE INTERFACE       ============================ //
// ================================================================================ //
    sprintf(rundetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_", Ncol, nt, nx, ny);
  #if dim==4
    sprintf(rundetails+strlen(rundetails),"nz_%d_", nz);
  #endif
    sprintf(rundetails+strlen(rundetails),"beta_%12.10lf_", beta);

    sprintf(rundetails+strlen(rundetails),"thread_%d", thread_number);
    
    sprintf(measfilename,"datadir/measurements_already_done_%s.dat",rundetails);
    measfile.open(measfilename,fstream::in);
  #if output_stream_index != 0
    sprintf(plaquettefilename,"datadir/plaquette_%s.dat", rundetails);
    plaquettefile.open(plaquettefilename
  #ifdef __wanna_append__
    , ios_base::app
  #endif
    );
    plaquettefile << setprecision(precision_digits);
  #endif
    
    output_stream << "# ======================================================================\n";
    output_stream << "# Local time and date at run start: " << asctime(timeinfo);
    output_stream << "# ----------------------------------------------------------------------\n";
    output_stream << "# Lattice setup:                 "          << endl;
    output_stream << "#                         Ncol = " << Ncol  << endl;
    output_stream << "#                          dim = " << dim   << endl; 
    output_stream << "#                                "          << endl;
    output_stream << "# Input parameters read from the " << inputfilename << " file:\n";
    output_stream << "#                         beta = " << beta  << endl;
    output_stream << "#                           nt = " << nt    << endl;
    output_stream << "#                           nx = " << nx << endl;
    output_stream << "#                           ny = " << ny << endl;
  #if dim>3
    output_stream << "#                           nz = " << nz << endl;
  #endif
    output_stream << "# ----------------------------------------------------------------------\n";
    output_stream << "#                         mode = " << mode << endl;
    output_stream << "#          thermalization time = " << thermalization_time << endl;
    output_stream << "#       reunitarization period = " << reunitarization_period << endl;
    output_stream << "#       number of measurements = " << how_many_measurements << endl;
    output_stream << "# updates between measurements = " << updates_between_measurements << endl;
    output_stream << "#                   start type = " << type_of_start << " ";
    switch (type_of_start) {
      case 0: {
        output_stream << "(random start)";
        break;
      }
      case 1: {
        output_stream << "(fully ordered configuration)";
        break;
      }
      case 2: {
        output_stream << "(reading saved configuration)";
        break;
      }
      default: {
        cout << "Error! Illegal start type";
        exit(0);
        break;
      }
    }
    output_stream << endl;
    output_stream << "# ======================================================================\n";
// ================================================================================ //
// ================================================================================ //


  /*
  if (measfile.is_open()) {
    measfile >> measurements_already_done;
    output_stream << "# " << measurements_already_done << " measurements already done in this thread" << endl;
  }
  else {
    output_stream << "# No previous measurements in this thread" << endl;
  }
  measfile.close();

  output_stream << "###############################################" << endl;
  */

	
  //simulation(mode);



/*
#ifdef __parallel__
  MPI::Finalize();
#endif
*/



  simulation();








  deallocate_arrays();


  return 0;

}
