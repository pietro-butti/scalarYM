double plaquette() {

  int site, mu, i;
  dc pltot=dc(0.,0.);

// When computing the Jarzynski SF, we exclude the (vanishing) contribution from
// spacelike plaquettes in the t=0 timeslice, as well as from all
// plaquettes in the t=Nt-1 timeslice. Note that the contribution from
// space-time plaquettes between the t=0 and t=1 timeslice is included, instead.
  
  for (site=
#ifdef __wanna_Jarzynski_SF__
            first_nt_equals_1_site_index
#else
            0
#endif
            ;site<
#ifdef __wanna_Jarzynski_SF__
            excluded_final_site_index
#else
            nsites
#endif
            ;site++)
  for (mu=0;mu<dim;mu++) {
    for (i=0;i<Ncolsquare;i++) {
      u4[i] = ufield[(mu*nsites+site)*Ncolsquare+i];
    }
    pstaple( staple, site, mu, 0, 0);
    mult_C_equals_ABdagger(tempor, u4, staple);
    for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
      pltot+=tempor[i];
    }
  }

#ifdef __wanna_Jarzynski_SF__
  for (site=0;site<first_nt_equals_1_site_index;site++) {
// Here we are setting mu=0 (space-time plaquettes), hence we can
// redefine this variable to save some multiplications
    mu=site*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      u4[i] = ufield[mu+i];
    }
    pstaple( staple, site, 0, 0, 0);
    mult_C_equals_ABdagger(tempor, u4, staple);
    for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
      pltot+=tempor[i];
    }
  }
#endif

#ifdef __wanna_Jarzynski_SF__
  return real(pltot)/(Ncol*(dim-1.)*((nx*ny
#if dim>3    
                                           *nz
#endif
                                           )*(dim*(nt-2)+1.)));
#else
  return real(pltot)/(Ncol*dim*(dim-1.)*nsites);
#endif

}
