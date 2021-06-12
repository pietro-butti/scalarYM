// #include "print_matrix.cc"
// #include "print_diagonal_matrix.cc"
// #include "print_2x2_matrix.cc"

#if Ncol==2
  #include "multiplication_routines_for_SU2.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU2((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU2((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU2((C),(A),(B))
  #include "norm_su2.cc"
#endif

#if Ncol==3
  #include "multiplication_routines_for_SU3.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU3((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU3((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU3((C),(A),(B))
  #include "norm_su3.cc"
#endif

#if Ncol==4
  #include "multiplication_routines_for_SU4.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU4((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU4((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU4((C),(A),(B))
  #include "norm_su4.cc"
#endif

#if Ncol==5
  #include "multiplication_routines_for_SU5.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU5((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU5((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU5((C),(A),(B))
  #include "norm_su5.cc"
#endif

#if Ncol==6
  #include "multiplication_routines_for_SU6.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU6((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU6((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU6((C),(A),(B))
  #include "norm_su6.cc"
#endif

#if Ncol==7
  #include "multiplication_routines_for_SU7.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU7((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU7((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU7((C),(A),(B))
  #include "norm_su7.cc"
#endif

#if Ncol==8
  #include "multiplication_routines_for_SU8.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU8((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU8((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU8((C),(A),(B))
  #include "norm_su8.cc"
#endif

#if Ncol==9
  #include "multiplication_routines_for_SU9.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU9((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU9((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU9((C),(A),(B))
  #include "norm_su9.cc"
#endif

#if Ncol==10
  #include "multiplication_routines_for_SU10.cc"
  #define mult_C_equals_AB(C,A,B) mult_C_equals_AB_for_SU10((C),(A),(B))
  #define mult_C_equals_ABdagger(C,A,B) mult_C_equals_ABdagger_for_SU10((C),(A),(B))
  #define mult_C_equals_AdaggerB(C,A,B) mult_C_equals_AdaggerB_for_SU10((C),(A),(B))
  #include "norm_su10.cc"
#endif

#include "read_last_conf.cc"
#include "initialize.cc"

// Optional:
#include "check_unitarity.cc"
#include "locked_gauge_transform.cc"
#include "gauge_transform.cc"
////////////////

#include "pstaple.cc"
#include "nstaple.cc"

#ifdef __wanna_improvement__
  #include "prectangularstaple.cc"
  #include "nrectangularstaple.cc"
#endif

#include "plaquette.cc"
#include "space_plaq.cc"
#include "time_plaq.cc"


#include "heat_bath_for_one_link.cc"
#include "heat_bath.cc"
#include "overrelaxation_for_one_link.cc"
#include "overrelaxation.cc"


#ifdef __wanna_multilevel__
  #include "measure_correlators_from_averaged_timelike_blocks.cc"
  #include "slab_average.cc"
#endif

#ifdef __wanna_Jarzynski_SF__
  #include "update.cc"
  #include "update_first_boundary.cc"
  #include "update_second_boundary.cc"
  #include "reunitarize.cc"
  #include "compute_action.cc"
  #include "compute_action_from_first_boundary.cc"
  #include "compute_action_from_second_boundary.cc"
  #include "measure_Jarzynski_SF.cc"
#endif

#include "measure_Polyakov_loops.cc"

#ifdef __wanna_smearing__
  #include "smear.cc"
#endif

#ifdef __wanna_Wilson_loops__
  #include "measure_Wilson_loops.cc"
#endif

#ifdef __wanna_trace_deformation__
  #include "compute_action.cc"
  #include "open_Polyakov_loop.cc"
  #include "generate_SUN_matrix_near_identity.cc"
  #include "metropolis.cc"
#endif

#ifndef __wanna_Jarzynski_SF__
  #include "update.cc"
  #include "reunitarize.cc"
#endif


#include "scalarym_metropolis.cc"  //nuovo
#include "scalarym_getlattice.cc"  //nuovo

#include "time_to_stop.cc"
#include "save_last_conf.cc"
#include "terminate_run.cc"
#include "thermalize.cc"
#include "measurements.cc"
#include "simulation.cc"
#include "deallocate_arrays.cc"


