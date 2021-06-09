void open_Polyakov_loop(dc *open_loop, int dir, int initial_site) {

// Computes the ordered product of the link matrices in the direction \hat{dir},
// from the site x+a\hat{dir} to site x, and stores it into the 
// complex NxN matrix open_loop 
  
  int i, link_counter;
  int number_of_factor_links;
  if (dir==0) { 
    number_of_factor_links=nt-1;
  }
  if (dir==1) { 
    number_of_factor_links=nx-1;
  }
  if (dir==2) { 
    number_of_factor_links=ny-1;
  }
#if dim>3
  if (dir==3) { 
    number_of_factor_links=nz-1;
  }
#endif

  int site=neighbor_plus[dir*nsites+initial_site];
  
  if (number_of_factor_links==1) {
    for (i=0;i<Ncolsquare;i++) {
      open_loop[i]=ufield[(dir*nsites+site)*Ncolsquare+i];
    }
  }
  else {
  
    int next_site=neighbor_plus[dir*nsites+site];
  
    for (i=0;i<Ncolsquare;i++) {
      u1[i] = ufield[(dir*nsites+site)*Ncolsquare+i];
      u2[i] = ufield[(dir*nsites+next_site)*Ncolsquare+i];
    }
    mult_C_equals_AB(tempor, u1, u2);
  
    for (link_counter=2;link_counter<number_of_factor_links;link_counter++) {
      next_site=neighbor_plus[dir*nsites+next_site];
      for (i=0;i<Ncolsquare;i++) {
        u1[i] = tempor[i];
        u2[i] = ufield[(dir*nsites+next_site)*Ncolsquare+i];
      }
      mult_C_equals_AB(tempor, u1, u2);
    }
  
// cout << "initial site = " << initial_site << "; final site = " << neighbor_plus[dir*nsites+next_site] << endl;
 
    for (i=0;i<Ncolsquare;i++) {
      open_loop[i]=tempor[i];
    }
  }
}
