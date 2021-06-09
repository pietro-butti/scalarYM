#define Ncol 2
const int Ncolsquare=Ncol*Ncol;

#define dim 3

#define how_many_hb 1
#define how_many_or 3

#define tol (1.5e-31)

#define how_many_irreps 1 // 12

// #define __debugging_mode__
// #define __parallel__
#define __wanna_save_last_conf__
// #define __wanna_improvement__
// #define __wanna_tadpole__
// #define __wanna_multilevel__
// #define __wanna_zero_mom_loops_and_correlators__
// #define __wanna_nonzero_mom_correlators__
// #define __wanna_Wilson_loops__
// #define __wanna_smearing__
// #define __wanna_trace_deformation__
#define __wanna_append__
// #define __wanna_Jarzynski_SF__

#define output_stream_index 0
// output_stream_index=0 directs plaquette output to screen; otherwise to plaquette file
#if output_stream_index == 0
#define output_stream cout
#else
#define output_stream plaquettefile
#endif

#ifdef __parallel__
const int parallel_thread_offset=0;
#endif

#ifdef __wanna_smearing__
#define first_smearing_direction 1
const double smearing_alpha=0.3;
const double link_smearing_coeff=1.-smearing_alpha;
const double staples_smearing_coeff=smearing_alpha/(2.*(dim-1-first_smearing_direction));
const int smearing_iterations=5;
const double convergence_tolerance=1.e-10;
long int max_number_of_projections=50;
#endif

#define first_mu 0
// 0: averages Wilson loops over all planes (for T=0 lattices)
// 1: averages spatial Wilson loops only (for lattices at finite temperature)

const int how_many_bins = 10;
int binsize;

#ifdef __wanna_improvement__
double plaquettefactor=(5./3.);
double rectanglefactor=-1./12.;
#endif

#if defined __wanna_nonzero_mom_correlators__ || defined __wanna_multilevel__
const int distances=6;
const int delta_x[distances]={0,1,2,3,4,5};//{0,1,2,3,4,5,6,7,8,9,10,11,12};
const int delta_y[distances]={0,0,0,0,0,0};//{0,0,0,0,0,0,0,0,0,0,0,0,0};
int spatial_volume;
#if dim==4
const int delta_z[distances]={0,0,0,0,0,0};//{0,0,0,0,0,0,0,0,0,0,0,0,0};
#endif
double Pol_correlator_normalization_factor;
#endif

#ifdef __wanna_multilevel__
const int slab_size=3;
const int how_many_slab_updates=10;
const double inverse_of_how_many_slab_updates=1./how_many_slab_updates;
// int how_many_slab;
#endif

#ifdef __wanna_Jarzynski_SF__
bool direct_transformation=true;
// const int how_many_Jarzynski_trajectories=10;
// const int how_many_Jarzynski_updates=10;
const int Jarzynski_thermalization=3;
#endif

#ifdef __wanna_trace_deformation__
#define how_many_metro 1
#endif

#define quantibins how_many_bins
