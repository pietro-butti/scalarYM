double time_plaq() {

  int mu, site, i, primo, secondo;
  dc pltot=dc(0.,0.);

  for (site=0;site<nsites;site++)
  for (mu=1;mu<dim;mu++) {

    primo=neighbor_plus[mu*nsites+site];
    secondo=neighbor_plus[site];

    for (i=0;i<Ncolsquare;i++) {
      u1[i] = ufield[(mu*nsites+site)*Ncolsquare+i];
      u2[i] = ufield[(primo)*Ncolsquare+i];
      u3[i] = ufield[(site)*Ncolsquare+i];
      u4[i] = ufield[(mu*nsites+secondo)*Ncolsquare+i];
    }

    mult_C_equals_AB(tempor,u1,u2);
    mult_C_equals_AB(tempor2,u3,u4);

    for (i=0;i<Ncolsquare;i++) {
      pltot+=tempor[i]*conj(tempor2[i]);
    }
  }

  return real(pltot/(Ncol*(dim-1.)*nsites));

}
