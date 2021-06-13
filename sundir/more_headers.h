#ifndef __more_headers__
#define __more_headers__

#ifdef __parallel__
#include <mpi.h>
int id, p;
#endif

using namespace std;

typedef complex<double> dc;

int max_Wilson_loop_size;
const int min_Wilson_loop_size=1;

char *input_file_string = new char[maximum_filename_length];
char *inputfilename = new char[maximum_filename_length];
char *rundetails = new char[maximum_filename_length];
char *measfilename = new char[maximum_filename_length];
char *plaquettefilename = new char[maximum_filename_length];

unsigned long tempus;
unsigned long measurements_already_done=0;
unsigned long measurements_done_until_now=0;

int *neighbor_plus;
int *neighbor_minus;

const int Ncol_plus_one=Ncol+1;

#ifdef __wanna_Wilson_loops__
ofstream Wilson_loop_file;
#ifdef __wanna_smearing__
ofstream smeared_Wilson_loop_file;
#endif
#endif


#ifdef __wanna_zero_mom_loops_and_correlators__
ofstream average_Polyakov_file;
ofstream zero_mom_Polyakov_correlators_file;
#ifdef __wanna_nonzero_mom_correlators__
ofstream nonzero_mom_Polyakov_correlators_file;
#endif
#endif


#ifdef __wanna_multilevel__
bool *locked_link;
int atb_size;
dc *blocks_in_this_slab;
const int Ncolfourth=Ncolsquare*Ncolsquare;
int spatial_volume_times_distances;
dc *averaged_timelike_blocks=new dc[atb_size];
dc *largetemp1=new dc[Ncolfourth];
dc *largetemp2=new dc[Ncolfourth];
dc *largetemp3=new dc[Ncolfourth];
char multilevel_correlators_filename[maximum_filename_length];
ofstream multilevel_correlators_file;
#endif

double hb_tried=0.;
double hb_accepted=0.;
double or_tried=0.;
double or_accepted=0.;

dc *ufield;
dc *v = new dc[Ncolsquare];

#ifdef __wanna_smearing__
dc *smeared_ufield;
dc *tempfield;
#endif

#ifdef __wanna_Jarzynski_SF__
dc *stored_ufield;
dc *C0 = new dc[Ncolsquare];
dc *C1 = new dc[Ncolsquare];
dc *Ctemp = new dc[Ncolsquare];
bool *locked_link;
ofstream Jarzynski_SF_file;
unsigned long int n_MC_time;
int nsteps_minus_one;
double one_over_nsteps;
int first_nt_equals_1_site_index;
int first_nt_equals_2_site_index;
int centralsite, centralsite_minus_one;
int starting_site_index;
int excluded_final_site_index;
int site_index_times_Ncolsquare;
double Wilson_prefactor;
double *phi0 = new double[Ncol];
double *newphi0 = new double[Ncol];
double *phi1 = new double[Ncol];
double *newphi1 = new double[Ncol];
double *phitemp = new double[Ncol];
int step_index;
double division_by_L;
#endif

int site;
dc *staple = new dc[Ncolsquare];
dc *ak = new dc[Ncolsquare];
dc *bigesse = new dc[Ncolsquare];
dc *oldu = new dc[Ncolsquare];
dc *tempor = new dc[Ncolsquare];
dc *tempor2 = new dc[Ncolsquare];
dc *newblock = new dc[Ncolsquare];
dc *rkappa = new dc[4];

dc *u1 = new dc[Ncolsquare];
dc *u2 = new dc[Ncolsquare];
dc *u3 = new dc[Ncolsquare];
dc *u4 = new dc[Ncolsquare];
dc *product = new dc[Ncolsquare];


#ifdef __wanna_improvement__
dc *rectangularstaple = new dc[Ncolsquare];
dc *u5 = new dc[Ncolsquare];
dc *u6 = new dc[Ncolsquare];
dc *u7 = new dc[Ncolsquare];
dc *u8 = new dc[Ncolsquare];
dc *u9 = new dc[Ncolsquare];
dc *u10 = new dc[Ncolsquare];
dc *u11 = new dc[Ncolsquare];
dc *u12 = new dc[Ncolsquare];
dc *product2 = new dc[Ncolsquare];
dc *product3 = new dc[Ncolsquare];
#ifdef __wanna_tadpole__
double averageplaquette;
#endif
#endif

#ifdef __wanna_trace_deformation__
int max_t, max_x, max_y;
#if dim>3
int max_z;
#endif
double read_alpha;
dc *open_loop = new dc[Ncolsquare];
dc *closed_loop = new dc[Ncolsquare];
dc *aux_matrix = new dc[Ncolsquare];
dc *SUN_matrix = new dc[Ncolsquare];
double oldaction_minus_newaction;
double twiceamplitude=0.2;
dc temp_contribution=dc(0.,0.);
bool *is_compactified = new bool[dim];
int compactified_direction_size;
double Wilson_prefactor;
double double_trace_prefactor;
const int how_many_trace_deformations = (int) floor(0.5*Ncol);
double *alpha_coeff = new double[how_many_trace_deformations];
int metro_tried, metro_accepted;
const int max_tried = 10000;
double metro_acceptance;
double amplitude_resizing_factor=0.8;
#endif

dc *zero_mom_loops;
dc *zero_mom_corr;

int thread_number;

int nsites, nlinks, ufielddimension, spatial_volume;
fstream measfile;
ofstream plaquettefile;

// Parameters to be read from the parameter file:

double beta;
int nt;
int nx;
int ny;
#if dim>3
int nz;
#endif
int mode;
long int thermalization_time;
int reunitarization_period;
long int how_many_measurements;
int updates_between_measurements;
int type_of_start;
#ifdef __wanna_Jarzynski_SF__
int nsteps; // The interval is divided in nsteps steps
char Jarzynski_direction;
double eta;
double deltaeta;
#if Ncol==3
double gnu;
double deltagnu;
#endif
#endif


/////////////////////////////////////////////////////////////


// nuovi
dc pauli[12] = { 
  dc(0.,0.),
  dc(1.,0.),
  dc(1.,0.),
  dc(0.,0.),
  dc(0.,0.),
  dc(0.,-1.),
  dc(0.,1.),
  dc(0.,0.),
  dc(1.,0.),
  dc(0.,0.),
  dc(0.,0.),
  dc(-1.,0.)
};
#endif