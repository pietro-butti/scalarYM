void deallocate_arrays() {
  
  delete [] neighbor_plus;
  delete [] neighbor_minus;
  delete [] ufield;
  delete [] v;
  delete [] input_file_string;
  delete [] inputfilename;
  delete [] rundetails;
  delete [] measfilename;
  delete [] plaquettefilename;
  delete [] staple;
  delete [] ak;
  delete [] bigesse;
  delete [] oldu;
  delete [] tempor;
  delete [] tempor2;
  delete [] newblock;
  delete [] rkappa;
  delete [] u1;
  delete [] u2;
  delete [] u3;
  delete [] u4;
  delete [] product;


#ifdef __wanna_improvement__
  delete [] rectangularstaple;
  delete [] u5;
  delete [] u6;
  delete [] u7;
  delete [] u8;
  delete [] u9;
  delete [] u10;
  delete [] u11;
  delete [] u12;
  delete [] product2;
  delete [] product3;
#endif

  
#ifdef __wanna_trace_deformation__
  delete [] open_loop;
  delete [] closed_loop;
  delete [] aux_matrix;
  delete [] SUN_matrix;
  delete [] is_compactified;
  delete [] alpha_coeff;
#endif
  
  
#ifdef __wanna_Wilson_loops__
  delete [] Wilson_filename;
#ifdef __wanna_smearing__
  delete [] smeared_Wilson_filename
#endif
#endif

#ifdef __wanna_multilevel__
  delete [] locked_link;
  delete [] blocks_in_this_slab;
  delete [] averaged_timelike_blocks;
  delete [] largetemp1;
  delete [] largetemp2;
  delete [] largetemp3;
#endif
  
#ifdef __wanna_smearing__
  delete [] smeared_ufield;
  delete [] tempfield;
#endif
  
#ifdef __wanna_Jarzynski_SF__
  delete [] stored_ufield;
  delete [] locked_link;
  delete [] C0;
  delete [] C1;
  delete [] Ctemp;
  delete [] phi0;
  delete [] newphi0;
  delete [] phi1;
  delete [] newphi1;
  delete [] phitemp;
#endif

#ifdef __wanna_zero_mom_loops_and_correlators__
  delete [] zero_mom_loops;
  delete [] zero_mom_corr;
#endif

}