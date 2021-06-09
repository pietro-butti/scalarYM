double space_plaq() {

  int mu, nu, site, i, primo, secondo;
  dc pltot = dc(0.,0.);

  for (site=0;site<nsites;site++)
  for (mu=1;mu<dim-1;mu++)
  for (nu=mu+1;nu<dim;nu++) {

    primo=neighbor_plus[mu*nsites+site];
    secondo=neighbor_plus[nu*nsites+site];

    for (i=0;i<Ncolsquare;i++) {
      u1[i] = ufield[(mu*nsites+site)*Ncolsquare+i];
      u2[i] = ufield[(nu*nsites+primo)*Ncolsquare+i];
      u3[i] = ufield[(nu*nsites+site)*Ncolsquare+i];
      u4[i] = ufield[(mu*nsites+secondo)*Ncolsquare+i];
    }

    mult_C_equals_AB(tempor,u1,u2);
    mult_C_equals_AB(tempor2,u3,u4);

    for (i=0;i<Ncolsquare;i++) {
      pltot+=tempor[i]*conj(tempor2[i]);
    }
  }

  return real(2.*pltot/(Ncol*(dim-1.)*(dim-2.)*nsites));

}
