void nrectangularstaple(dc *stot, int site, int mu, int first_direction, int use_smeared) {

  int jmu, jminusmu;
  int nu, knu, knumu, knutwomu, knuminusmu, ktwonu, ktwonumu;
  int i, j;

/* Given the link from site to site+mu, */
/* this routine calculates the sum of the */
/* 2x1 rectangular staples in the negative nu directions */
/* relevant in the update of U_mu(site), namely: */
/*   _     _              */
/* |__| + |__| + | |      */
/*               | |      */
/*                -       */
/*                        */
/* and adds them to stot. */
 
// Important: As opposed to prectangularstaple, here the result is ADDED to stot!


  jmu=neighbor_plus[mu*nsites+site];
  jminusmu=neighbor_minus[mu*nsites+site];
  for (nu=first_direction;nu<dim;nu++)
  if (nu!=mu) {
    knu=neighbor_minus[nu*nsites+site];
    knumu=neighbor_plus[mu*nsites+knu];
    knutwomu=neighbor_plus[mu*nsites+knumu];
    knuminusmu=neighbor_minus[mu*nsites+knu];
    ktwonu=neighbor_minus[nu*nsites+knu];
    ktwonumu=neighbor_plus[mu*nsites+ktwonu];
    
    for (i=0;i<Ncol;i++)
    for (j=0;j<Ncol;j++) {
      if (use_smeared==0) {
        u1dag[i*Ncol+j]=conj(ufield[(nu*nsites+knu)*Ncolsquare+j*Ncol+i]);
        u2[i*Ncol+j]=ufield[(mu*nsites+knu)*Ncolsquare+i*Ncol+j];
        u3[i*Ncol+j]=ufield[(nu*nsites+knumu)*Ncolsquare+i*Ncol+j];
        u4[i*Ncol+j]=ufield[(mu*nsites+knumu)*Ncolsquare+i*Ncol+j];
        u5[i*Ncol+j]=ufield[(nu*nsites+knutwomu)*Ncolsquare+i*Ncol+j];
        u6dag[i*Ncol+j]=conj(ufield[(mu*nsites+jmu)*Ncolsquare+j*Ncol+i]);
        u7dag[i*Ncol+j]=conj(ufield[(nu*nsites+ktwonu)*Ncolsquare+j*Ncol+i]);
        u8[i*Ncol+j]=ufield[(mu*nsites+ktwonu)*Ncolsquare+i*Ncol+j];
        u9[i*Ncol+j]=ufield[(nu*nsites+ktwonumu)*Ncolsquare+i*Ncol+j];
        u10dag[i*Ncol+j]=conj(ufield[(mu*nsites+jminusmu)*Ncolsquare+j*Ncol+i]);
        u11dag[i*Ncol+j]=conj(ufield[(nu*nsites+knuminusmu)*Ncolsquare+j*Ncol+i]);
        u12[i*Ncol+j]=ufield[(mu*nsites+knuminusmu)*Ncolsquare+i*Ncol+j];
      }
#ifdef __wanna_smearing__
      else {
        u1dag[i*Ncol+j]=conj(smeared_ufield[(nu*nsites+knu)*Ncolsquare+j*Ncol+i]);
        u2[i*Ncol+j]=smeared_ufield[(mu*nsites+knu)*Ncolsquare+i*Ncol+j];
        u3[i*Ncol+j]=smeared_ufield[(nu*nsites+knumu)*Ncolsquare+i*Ncol+j];
        u4[i*Ncol+j]=smeared_ufield[(mu*nsites+knumu)*Ncolsquare+i*Ncol+j];
        u5[i*Ncol+j]=smeared_ufield[(nu*nsites+knutwomu)*Ncolsquare+i*Ncol+j];
        u6dag[i*Ncol+j]=conj(smeared_ufield[(mu*nsites+jmu)*Ncolsquare+j*Ncol+i]);
        u7dag[i*Ncol+j]=conj(smeared_ufield[(nu*nsites+ktwonu)*Ncolsquare+j*Ncol+i]);
        u8[i*Ncol+j]=smeared_ufield[(mu*nsites+ktwonu)*Ncolsquare+i*Ncol+j];
        u9[i*Ncol+j]=smeared_ufield[(nu*nsites+ktwonumu)*Ncolsquare+i*Ncol+j];
        u10dag[i*Ncol+j]=conj(smeared_ufield[(mu*nsites+jminusmu)*Ncolsquare+j*Ncol+i]);
        u11dag[i*Ncol+j]=conj(smeared_ufield[(nu*nsites+knuminusmu)*Ncolsquare+j*Ncol+i]);
        u12[i*Ncol+j]=smeared_ufield[(mu*nsites+knuminusmu)*Ncolsquare+i*Ncol+j];
      }
#endif
    }
    
    mult_C_equals_AB(tempor,u1dag,u2);
    mult_C_equals_AB(tempor2,tempor,u4);
    mult_C_equals_AB(tempor,tempor2,u5);
    mult_C_equals_AB(product,tempor,u6dag);
    
    mult_C_equals_AB(tempor,u10dag,u11dag);
    mult_C_equals_AB(tempor2,tempor,u12);
    mult_C_equals_AB(tempor,tempor2,u2);
    mult_C_equals_AB(product2,tempor,u3);
    
    mult_C_equals_AB(tempor,u1dag,u7dag);
    mult_C_equals_AB(tempor2,tempor,u8);
    mult_C_equals_AB(tempor,tempor2,u9);
    mult_C_equals_AB(product3,tempor,u3);
    
    for (i=0;i<Ncolsquare;i++) {
      stot[i]+=product[i] + product2[i] + product3[i];
    }
  }

}
