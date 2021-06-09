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
  
#ifdef __wanna_trace_deformation__
#ifdef __wanna_multilevel__
  cout << "Error: The multilevel routines are not generically compatible with the trace deformation" << endl;
  exit(0);
#endif
#endif

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
    cout << "Usage for scalar mode: ./run.x <input_file_string> <thread_number>" << endl; 
    exit(0);
  };
  thread_number=atoi(argv[2]);

#endif

  time(&rawtime);
  ptm=gmtime(&rawtime);
  timeinfo=localtime(&rawtime);

  cout << setprecision(precision_digits);
  cout << fixed;

  char *inputfilename = new char[maximum_filename_length];

  sprintf(inputfilename,"%s_thread_%d.dat",input_file_string,thread_number);

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
  ufielddimension=nlinks*Ncolsquare;
  neighbor_plus = new int[nlinks];
  neighbor_minus = new int[nlinks];
#if defined __wanna_nonzero_mom_correlators__ || defined __wanna_multilevel__
  Pol_correlator_normalization_factor=1./(spatial_volume*Ncolsquare);
#endif
  
#ifdef __wanna_multilevel__
  how_many_slabs=nt/slab_size;
  locked_link= new bool[nlinks];
  atb_size=Ncolsquare*Ncolsquare*spatial_volume*distances;
  blocks_in_this_slab=new dc[atb_size];
  spatial_volume_times_distances=spatial_volume*distances;
  averaged_timelike_blocks=new dc[atb_size];
#endif

  ufield = new dc[ufielddimension];

#ifdef __wanna_smearing__
  smeared_ufield = new dc[ufielddimension];
  tempfield = new dc[ufielddimension];
#endif
  
  inputfile >> mode;
  inputfile >> thermalization_time;               
  inputfile >> reunitarization_period;            
  inputfile >> how_many_measurements;
  inputfile >> updates_between_measurements;
  inputfile >> type_of_start;
  
#ifdef __wanna_Jarzynski_SF__
  centralsite=nsites/2;
  centralsite_minus_one=centralsite-1;
  stored_ufield = new dc[ufielddimension];
  locked_link = new bool[nlinks];
  
  first_nt_equals_1_site_index=nx*ny;
  starting_site_index=(nt-2)*nx*ny;
  excluded_final_site_index=(nt-1)*nx*ny;
#if dim==4
  first_nt_equals_1_site_index*=nz;
  starting_site_index*=nz;
  excluded_final_site_index*=nz;
#endif
  first_nt_equals_2_site_index=2*first_nt_equals_1_site_index;

  division_by_L=1./nx;
  inputfile >> nsteps;
  nsteps_minus_one=nsteps-1;
  one_over_nsteps=1./nsteps;
  inputfile >> Jarzynski_direction;
  switch (Jarzynski_direction) {
    case 'd':
    case 'D': {
      direct_transformation=true;
      break;
    }
    case 'r':
    case 'R': {
      direct_transformation=false;
      break;
    }
    default: {
      cout << "Error while reading the input file: the direction must be 'd', 'D', 'r' or 'R'" << endl;
      exit(0);
      break;
    }
  }
  inputfile >> eta;
  inputfile >> deltaeta;
#if Ncol==3
  inputfile >> gnu;
  inputfile >> deltagnu;
#endif
#endif
#ifdef __wanna_improvement__
#ifdef __wanna_tadpole__
  inputfile >> averageplaquette;
#endif
#endif
#ifdef __wanna_trace_deformation__ 
  inputfile >> read_alpha;
#endif


#ifdef __wanna_Jarzynski_SF__
  Wilson_prefactor=-beta/Ncol;
#endif
  
#ifdef __wanna_trace_deformation__
  Wilson_prefactor=-beta/Ncol;
  for (int i=0; i<dim; i++) {
    is_compactified[i]=false;
  }
  is_compactified[1]=true;
  for (int i=0; i<how_many_trace_deformations; i++) {
    alpha_coeff[i]=0.;
  }
  alpha_coeff[0]=read_alpha;
#endif

  inputfile.close();

  sprintf(rundetails,"Ncol_%d_nt_%d_nx_%d_ny_%d_", Ncol, nt, nx, ny);
#if dim==4
  sprintf(rundetails+strlen(rundetails),"nz_%d_", nz);
#endif
  sprintf(rundetails+strlen(rundetails),"beta_%12.10lf_", beta);
#ifdef __wanna_Jarzynski_SF__
  sprintf(rundetails+strlen(rundetails),"nsteps_%d_%s_eta_%12.10lf_deltaeta_%12.10lf_", nsteps, (direct_transformation)?"direct":"reverse", eta, deltaeta);
#if Ncol==3
  sprintf(rundetails+strlen(rundetails),"gnu_%12.10lf_deltagnu_%12.10lf_", gnu, deltagnu);
#endif
#endif
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
  
  output_stream << "# Local time and date at run start: " << asctime(timeinfo);
  output_stream << "# Lattice setup:" << endl;
  output_stream << "# Ncol = " << Ncol << endl;
  output_stream << "# dim = " << dim << endl; 
  output_stream << "#" << endl;
  output_stream << "# Input parameters read from the " << inputfilename << " file:\n";
  output_stream << "# beta = " << beta << endl;
  output_stream << "# nt = " << nt 
#ifdef __wanna_Jarzynski_SF__
                << " sites (i.e. " << nt-1 << " links)"
#endif
                << endl;
  output_stream << "# nx = " << nx << endl;
  output_stream << "# ny = " << ny << endl;
#if dim>3
  output_stream << "# nz = " << nz << endl;
#endif
  output_stream << "# mode = " << mode << endl;
  output_stream << "# thermalization time = " << thermalization_time << endl;
  output_stream << "# reunitarization period = " << reunitarization_period << endl;
  output_stream << "# number of measurements = " << how_many_measurements << endl;
  output_stream << "# updates between measurements = " << updates_between_measurements << endl;
  output_stream << "# start type = " << type_of_start << " ";
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
#ifdef __wanna_Jarzynski_SF__
  output_stream << "# nsteps = " << nsteps << endl;
  output_stream << "# Jarzynski_direction = ";
  if (direct_transformation) {
    output_stream << "direct" << endl;
  }
  else {
    output_stream << "reverse" << endl;
  }
  output_stream << "# eta = " << eta << endl;
  output_stream << "# deltaeta = " << deltaeta << endl;
#if Ncol==3
  output_stream << "# gnu = " << gnu << endl;
  output_stream << "# deltagnu = " << deltagnu << endl;
#endif
#endif
#ifdef __wanna_improvement__
#ifdef __wanna_tadpole__
  output_stream << "# average plaquette from input file = " << averageplaquette << endl;
#endif
#endif
#ifdef __wanna_trace_deformation__ 
  output_stream << "# deformation parameter alpha = " << read_alpha << endl; 
#endif

#ifdef __wanna_Jarzynski_SF__
#if Ncol==2
  
  if (direct_transformation) {
    phi0[0]=-eta;
    newphi0[0]=phi0[0]-deltaeta;
    phi1[0]=eta-pi;
    newphi1[0]=phi1[0]+deltaeta;
  }
  else {
    newphi0[0]=-eta;
    phi0[0]=newphi0[0]-deltaeta;
    newphi1[0]=eta-pi;
    phi1[0]=newphi1[0]+deltaeta;
  }
  
  phi0[1]=-phi0[0];
  newphi0[1]=-newphi0[0];
  phi1[1]=-phi1[0];
  newphi1[1]=-newphi1[0];
#endif
#if Ncol==3
  
  if (direct_transformation) {
    phi0[0]=eta-(pi/3.);
    phi0[1]=eta*(gnu-0.5);
    phi0[2]=-eta*(gnu+0.5)+(pi/3.);
    newphi0[0]=eta+deltaeta-(pi/3.);
    newphi0[1]=(eta+deltaeta)*(gnu-0.5);
    newphi0[2]=-(eta+deltaeta)*(gnu+0.5)+(pi/3.);
    phi1[0]=-eta-pi;
    phi1[1]=(pi/3.)+eta*(gnu+0.5);
    phi1[2]=(2.*pi/3.)-eta*(gnu-0.5);
    newphi1[0]=-eta-deltaeta-pi;
    newphi1[1]=(pi/3.)+(eta+deltaeta)*(gnu+deltagnu+0.5);
    newphi1[2]=(2.*pi/3.)-(eta+deltaeta)*(gnu+deltagnu-0.5);
  }
  else {
    newphi0[0]=eta-(pi/3.);
    newphi0[1]=eta*(gnu-0.5);
    newphi0[2]=-eta*(gnu+0.5)+(pi/3.);
    phi0[0]=eta+deltaeta-(pi/3.);
    phi0[1]=(eta+deltaeta)*(gnu-0.5);
    phi0[2]=-(eta+deltaeta)*(gnu+0.5)+(pi/3.);
    newphi1[0]=-eta-pi;
    newphi1[1]=(pi/3.)+eta*(gnu+0.5);
    newphi1[2]=(2.*pi/3.)-eta*(gnu-0.5);
    phi1[0]=-(eta+deltaeta)-pi;
    phi1[1]=(pi/3.)+(eta+deltaeta)*(gnu+deltagnu+0.5);
    phi1[2]=(2.*pi/3.)-(eta+deltaeta)*(gnu+deltagnu-0.5);
  }
#endif
#endif
  
  
  if (measfile.is_open()) {
    measfile >> measurements_already_done;
    output_stream << "# " << measurements_already_done << " measurements already done in this thread" << endl;
  }
  else {
    output_stream << "# No previous measurements in this thread" << endl;
  }
  measfile.close();

  output_stream << "###############################################" << endl;

#ifdef __wanna_improvement__
#ifdef __wanna_tadpole__
  rectanglefactor/=sqrt(averageplaquette);
#endif
#endif

#ifdef __wanna_zero_mom_loops_and_correlators__
  zero_mom_loops = new dc[how_many_measurements*nx];
  zero_mom_corr = new dc[how_many_measurements*nx];
#endif
  
  simulation(mode);

#ifdef __parallel__
  MPI::Finalize();
#endif

  deallocate_arrays();
  
  return 0;

}
