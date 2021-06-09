void pstaple(dc *stot, int site, int mu, int first_direction, int use_smeared) {

// Computes the staples from (site) to (site+a mu) in the positive 
// directions (starting from first_direction) and stores their sum into stot.
    
  int jmu, knu, nu, i, aux_ind1, aux_ind2, aux_ind3;
  
  jmu=neighbor_plus[mu*nsites+site];

  for (i=0;i<Ncolsquare;i++) {
    stot[i]=dc(0.,0.);
  }

  for (nu=first_direction;nu<dim;nu++)
  if (nu!=mu) {
    knu=neighbor_plus[nu*nsites+site];
    aux_ind1=(nsites*nu+site)*Ncolsquare;
    aux_ind2=(nsites*mu+knu)*Ncolsquare;
    aux_ind3=(nsites*nu+jmu)*Ncolsquare;
    for (i=0;i<Ncolsquare;i++) {
      if (use_smeared==0) {
        u1[i]=ufield[aux_ind1+i];
        u2[i]=ufield[aux_ind2+i];
        u3[i]=ufield[aux_ind3+i];
      }
#ifdef __wanna_smearing__
      else {
        u1[i]=smeared_ufield[aux_ind1+i];
        u2[i]=smeared_ufield[aux_ind2+i];
        u3[i]=smeared_ufield[aux_ind3+i];
      }
#endif
    }
    mult_C_equals_AB(product,u1,u2);
    mult_C_equals_ABdagger(u1,product,u3);
    for (i=0;i<Ncolsquare;i++) {
      stot[i]+=u1[i];
    }
  }

}
