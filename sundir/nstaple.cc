void nstaple(dc *stot, int site, int mu, int first_direction, int use_smeared) {

// Computes the staples from (site) to (site+a mu) in the negative 
// directions (starting from first_direction) and adds them to stot.
// 
// Important: As opposed to pstaple, here the result is ADDED to stot!
    
  int jmu, knu, nu, i, aux_ind1, aux_ind2, aux_ind3;


  for (nu=first_direction;nu<dim;nu++)
  if (nu!=mu) {
    knu=neighbor_minus[nu*nsites+site];
    jmu=neighbor_plus[mu*nsites+knu];
    aux_ind1=(nu*nsites+knu)*Ncolsquare;
    aux_ind2=(mu*nsites+knu)*Ncolsquare;
    aux_ind3=(nu*nsites+jmu)*Ncolsquare;
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
    mult_C_equals_AdaggerB(product,u1,u2);
    mult_C_equals_AB(u2,product,u3);
    for (i=0;i<Ncolsquare;i++) {
      stot[i]+=u2[i];
    }
  }

}
