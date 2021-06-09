void prectangularstaple(dc *stot, int site, int mu, int first_direction, int use_smeared) {

/* Given the link from site to site+mu, */
/* this routine calculates the sum of the */
/* 2x1 rectangular staples in the positive nu directions */
/* relevant in the update of U_mu(site), namely: */
/*                _         */
/*  __     __    | |        */
/* | _| + |_ | + | |        */
/*                          */
/* and stores it into stot. */
  
  int jmu, jtwomu, jminusmu; 
  int nu, knu, kminusmuplusnu, knumu, ktwonu;
  int i, j;
  int aux_ind1, aux_ind2, aux_ind3, aux_ind4, aux_ind5, aux_ind6;
  int aux_ind7, aux_ind8, aux_ind9, aux_ind10, aux_ind11, aux_ind12;

  jmu=neighbor_plus[mu*nsites+site];
  jtwomu=neighbor_plus[mu*nsites+jmu];
  jminusmu=neighbor_minus[mu*nsites+site];

  for (i=0;i<Ncolsquare;i++) {
    stot[i]=dc(0.,0.);
  }

  for (nu=first_direction;nu<dim;nu++)
  if (nu!=mu) {
    knu=neighbor_plus[nu*nsites+site];
    knumu=neighbor_plus[mu*nsites+knu];
    kminusmuplusnu=neighbor_minus[mu*nsites+knu];
    ktwonu=neighbor_plus[nu*nsites+knu];
    
    aux_ind1 =(nu*nsites+site)*Ncolsquare;
    aux_ind2 =(mu*nsites+knu)*Ncolsquare;
    aux_ind3 =(nu*nsites+jmu)*Ncolsquare);
    aux_ind4 =(mu*nsites+knumu)*Ncolsquare;
    aux_ind5 =(nu*nsites+jtwomu)*Ncolsquare;
    aux_ind6 =(mu*nsites+jmu)*Ncolsquare;
    aux_ind7 =(mu*nsites+jminusmu)*Ncolsquare;
    aux_ind8 =(nu*nsites+jminusmu)*Ncolsquare;
    aux_ind9 =(mu*nsites+kminusmuplusnu)*Ncolsquare;
    aux_ind10=(nu*nsites+knu)*Ncolsquare;
    aux_ind11=(mu*nsites+ktwonu)*Ncolsquare;
    aux_ind12=(nu*nsites+knumu)*Ncolsquare;
    
    for (i=0;i<Ncolsquare;i++) {
      if (use_smeared==0) {
        u1[i] =ufield[aux_ind1+i];
        u2[i] =ufield[aux_ind2+i];
        u3[i] =ufield[aux_ind3+i];
        u4[i] =ufield[aux_ind4+i];
        u5[i] =ufield[aux_ind5+i];
        u6[i] =ufield[aux_ind6+i];
        u7[i] =ufield[aux_ind7+i];
        u8[i] =ufield[aux_ind8+i];
        u9[i] =ufield[aux_ind9+i];
        u10[i]=ufield[aux_ind10+i];
        u11[i]=ufield[aux_ind11+i];
        u12[i]=ufield[aux_ind12+i];
      }
#ifdef __wanna_smearing__
      else {
        u1[i] =smeared_ufield[aux_ind1+i];
        u2[i] =smeared_ufield[aux_ind2+i];
        u3[i] =smeared_ufield[aux_ind3+i];
        u4[i] =smeared_ufield[aux_ind4+i];
        u5[i] =smeared_ufield[aux_ind5+i];
        u6[i] =smeared_ufield[aux_ind6+i];
        u7[i] =smeared_ufield[aux_ind7+i];
        u8[i] =smeared_ufield[aux_ind8+i];
        u9[i] =smeared_ufield[aux_ind9+i];
        u10[i]=smeared_ufield[aux_ind10+i];
        u11[i]=smeared_ufield[aux_ind11+i];
        u12[i]=smeared_ufield[aux_ind12+i];
      }
#endif
    }
    
    mult_C_equals_AB(tempor,u1,u2);
    mult_C_equals_AB(tempor2,tempor,u4);
    mult_C_equals_ABdagger(tempor,tempor2,u5);
    mult_C_equals_ABdagger(product,tempor,u6);

    mult_C_equals_AdaggerB(tempor,u7,u8);
    mult_C_equals_AB(tempor2,tempor,u9);
    mult_C_equals_AB(tempor,tempor2,u2);
    mult_C_equals_ABdagger(product2,tempor,u3);
    
    mult_C_equals_AB(tempor,u1,u10);
    mult_C_equals_AB(tempor2,tempor,u11);
    mult_C_equals_ABdagger(tempor,tempor2,u12);
    mult_C_equals_ABdagger(product3,tempor,u3);
    
    for (i=0;i<Ncolsquare;i++) {
      stot[i]+=product[i]+product2[i]+product3[i];
    }
    
  }

}
