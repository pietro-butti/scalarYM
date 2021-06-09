void check_unitarity() {

  const double tolerance=1.e-10;
  int site, dir, i;

  for (site=0;site<nsites;site++)
  for (dir=0;dir<dim;dir++) {
    for (i=0;i<Ncolsquare;i++) {
      u1[i]=ufield[(dir*nsites+site)*Ncolsquare+i];
      u2[i]=u1[i];
    }
    mult_C_equals_ABdagger(tempor,u1,u2);
    mult_C_equals_AdaggerB(tempor2,u1,u2);
    for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
      tempor[i]-=dc(1.,0.);
      tempor2[i]-=dc(1.,0.);
    }
    for (i=0;i<Ncolsquare;i++) {
      if ( (norm(tempor[i])>tolerance) ||
	   (norm(tempor2[i])>tolerance) 
         ) {
	printf("ERRORE! %d %d\n",site,dir);
	exit(0);
      };
    }
  }
}
