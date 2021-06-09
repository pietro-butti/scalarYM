void update_second_boundary(unsigned long int step_index) {

  int site, dir, i, ind;

  for (i=0;i<Ncolsquare;i++) {
    Ctemp[i]=dc(0.,0.);
  }

  for (i=0; i<Ncol; i++) {
    phitemp[i]=phi1[i]+step_index*one_over_nsteps*(newphi1[i]-phi1[i]);
    Ctemp[i*(Ncol_plus_one)]=exp(dc(0.,phitemp[i]));
  }
  
// cout << "(" << real(Ctemp[0]) << "," << imag(Ctemp[0]) << ")  ";
// cout << "(" << real(Ctemp[1]) << "," << imag(Ctemp[1]) << ")" << endl;
// cout << "(" << real(Ctemp[2]) << "," << imag(Ctemp[2]) << ")  "; 
// cout << "(" << real(Ctemp[3]) << "," << imag(Ctemp[3]) << ")" << endl;
// exit(0);

    
  for (site=excluded_final_site_index;site<nsites;site++)
  for (dir=0;dir<dim;dir++) {
    ind=(dir*nsites+site)*Ncolsquare;
// Set all links in the t=nt-1 slice to Ctemp:
    for (i=0;i<Ncolsquare;i++) {
      ufield[ind+i]=Ctemp[i];
    }
  }

}
