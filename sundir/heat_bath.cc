void heat_bath(bool locked_mode) {

  int site, dir;

#ifndef __wanna_Jarzynski_SF__
  for (site=0;site<nsites;site++)
  for (dir=0;dir<dim;dir++) {
    heat_bath_for_one_link(locked_mode, site, dir);
  }
#else
// When the Schroedinger functional is studied with the Jarzynski method,
// it is best to start updating the links starting from those closest to the boundaries,
// then moving towards the bulk of the lattice, so that the physical information
// about the new value of the boundary links can propagate throughout the lattice
// in one sweep. Moreover, we skip updating the links that are locked, i.e. the 
// spatial links at t=0 and all links at t=nt-1:
  for (site=0;site<first_nt_equals_1_site_index;site++) {
    heat_bath_for_one_link(locked_mode, site, 0);
  }
  for (site=first_nt_equals_1_site_index;site<centralsite;site++)
  for (dir=0;dir<dim;dir++) {
    heat_bath_for_one_link(locked_mode, site, dir);
  }
  for (site=excluded_final_site_index-1;site>centralsite_minus_one;site--)
  for (dir=0;dir<dim;dir++) {
    heat_bath_for_one_link(locked_mode, site, dir);
  }
#endif

}
