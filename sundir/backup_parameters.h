#define Ncol 4
const int Ncolsquare=Ncol*Ncol;
#define dim 4

#define nt 4
#define nspacelike 12

#define nx nspacelike
#define ny nspacelike

#if dim>3
#define nz nspacelike
const int nsites=nt*nx*ny*nz;
#else
const int nsites=nt*nx*ny;
#endif

const int nlinks=dim*nsites;
const int ufielddimension=nlinks*Ncolsquare;

#define how_many_hb 1
#define how_many_or 3

#define tol 1.5e-31

#define how_many_irreps 1 // 12

// #define __debugging_mode__

// #define __parallel__

// #define __wanna_save_last_conf__

#define __wanna_improvement__
// #define __wanna_tadpole__
// #define __wanna_multilevel__
#define __wanna_zero_mom_loops_and_correlators__
// #define __wanna_Wilson_loops__
// #define __wanna_smearing__
#define __wanna_q_hat__
// #define __wanna_trace_deformation__

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



#ifdef __wanna_multilevel__
const int maxdist=nx/2;
const int distances=maxdist+1;
const int slab_size=4;
const int how_many_slabs=nt/slab_size;
const int how_many_slab_updates=10;
const double inverse_of_how_many_slab_updates=1./how_many_slab_updates;
#if dim==3
const int spatial_volume=nx*ny;
#endif
#if dim==4
const int spatial_volume=nx*ny*nz;
#endif
#endif

#ifdef __wanna_trace_deformation__
#define how_many_metro 1
#endif


#define quantibins how_many_bins
