double compute_action_from_second_boundary() {

  int site_index, link_index_times_Ncolsquare, i;
//   int dir;
  double action=0.;

#ifdef __wanna_improvement_
  cout << "Error! This routine works only for the Wilson action!" << endl;
  exit(0);
#endif  

// Assuming that the Wilson (i.e., unimproved) action is used, as is the case here, 
// there are two types of contributions to the action from the second boundary:
//   
// - timelike plaquettes between t=nt-2 and t=nt-1, with weight w=1, and
// - spatial plaquettes on the second boundary (t=nt-1), with weight w=1/2.
//   
// Given that, throughout this work, spatially constant boundary conditions are
// always assumed, the latter contribution is identically vanishing.
// As a consequence, the action is obtained from the contribution of timelike
// plaquettes between t=nt-2 and t=nt-1:

  for (site_index=starting_site_index; site_index<excluded_final_site_index; site_index++) {
    // This cycle runs over all temporal links in the timeslice between t=nt-2 and t=nt-1:
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
  for (site_index=starting_site_index; site_index<excluded_final_site_index; site_index++)
  for (dir=0; dir<dim; dir++) {
    // This cycle runs over all links in the time-slice next to the second boundary:
    link_index_times_Ncolsquare=(dir*nsites+site_index)*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      oldu[i]=ufield[link_index_times_Ncolsquare+i];
    }
    action+=compute_action(site_index,dir,oldu);
// cout << action << endl;
  }
// cout << endl;
*/
  return action;

}
