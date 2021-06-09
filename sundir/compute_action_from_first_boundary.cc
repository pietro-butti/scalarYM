double compute_action_from_first_boundary() {

  int site_index, link_index_times_Ncolsquare, i;
//   int dir;
  double action=0.;

#ifdef __wanna_improvement_
  cout << "Error! This routine works only for the Wilson action!" << endl;
  exit(0);
#endif  

// Assuming that the Wilson (i.e., unimproved) action is used, as is the case here, 
// there are two types of contributions to the action from the first boundary:
//   
// - timelike plaquettes between t=0 and t=1, with weight w=1, and
// - spatial plaquettes on the first boundary (t=0), with weight w=1/2.
//   
// Given that, throughout this work, spatially constant boundary conditions are
// always assumed, the latter contribution is identically vanishing.
// As a consequence, the action is obtained from the contribution of timelike
// plaquettes between t=0 and t=1:

  for (site_index=0; site_index<first_nt_equals_1_site_index; site_index++) {
    // This cycle runs over all temporal links in the timeslice between t=0 and t=1:
    link_index_times_Ncolsquare=site_index*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      oldu[i]=ufield[link_index_times_Ncolsquare+i];
    }
    pstaple(staple,site_index,0,1,0);
    mult_C_equals_ABdagger(tempor,oldu,staple);
    for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
      action+=real(tempor[i]);
    }
  }
  action*=Wilson_prefactor;

  
  
// Old version:
/*
  // This cycle runs over all timelike links at t=0
  
  for (site_index=0; site_index<first_nt_equals_1_site_index; site_index++) {
    link_index_times_Ncolsquare=site_index*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      oldu[i]=ufield[link_index_times_Ncolsquare+i];
    }
    action+=compute_action(site_index,0,oldu);
  }
  
  // This cycle runs over all spatial links at t=1
  for (site_index=first_nt_equals_1_site_index; site_index<first_nt_equals_2_site_index; site_index++)
  for (dir=1; dir<dim; dir++) {
    link_index_times_Ncolsquare=(dir*nsites+site_index)*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      oldu[i]=ufield[link_index_times_Ncolsquare+i];
    }
    action+=compute_action(site_index,dir,oldu);
  }
 */ 
  return action;

}
