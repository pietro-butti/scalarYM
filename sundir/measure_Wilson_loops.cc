void measure_Wilson_loops(int configuration, std::ofstream& output_file, int smearing_mode) {

  int i, site, mu, nu, L_mu, L_nu;
  int length, next, check_final_site;
  dc average;
  
  
  for (L_mu=min_Wilson_loop_size; L_mu<=max_Wilson_loop_size; L_mu++)
  for (L_nu=min_Wilson_loop_size; L_nu<=max_Wilson_loop_size; L_nu++) {
    average = dc(0.,0.);
    for (site=0; site<nsites; site++)
    for (mu=(smearing_mode==0)?first_mu:0; mu<((smearing_mode==0)?(dim-1):1); mu++)
    for (nu=mu+1; nu<dim; nu++) {  
      next=site;
      for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
        if (smearing_mode==0) {
          tempor[i]=ufield[(mu*nsites+site)*Ncolsquare+i];
	}
	else  {
          tempor[i]=smeared_ufield[(mu*nsites+site)*Ncolsquare+i];
	}
#else
        tempor[i]=ufield[(mu*nsites+site)*Ncolsquare+i];
#endif
      }
      next=neighbor_plus[mu*nsites+next];
      
      for (length=1; length<L_mu; length++) {
        for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
          if (smearing_mode==0) {
            u1[i]=ufield[(mu*nsites+next)*Ncolsquare+i];
	  }
	  else {
            u1[i]=smeared_ufield[(mu*nsites+next)*Ncolsquare+i];
	  }
#else
          u1[i]=ufield[(mu*nsites+next)*Ncolsquare+i];
#endif
        }
        mult_C_equals_AB(product,tempor,u1);
        for (i=0; i<Ncolsquare; i++) {
          tempor[i]=product[i];
        }
        next=neighbor_plus[mu*nsites+next];
      }
      
      for (length=0; length<L_nu; length++) {
        for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
          if (smearing_mode==0) {
            u1[i]=ufield[(nu*nsites+next)*Ncolsquare+i];
	  }
	  else {
            u1[i]=smeared_ufield[(nu*nsites+next)*Ncolsquare+i];
	  }
#else
          u1[i]=ufield[(nu*nsites+next)*Ncolsquare+i];
#endif
        }
        mult_C_equals_AB(product,tempor,u1);
        for (i=0; i<Ncolsquare; i++) {
          tempor[i]=product[i];
        }
        next=neighbor_plus[nu*nsites+next];
      }
     
      for (i=0; i<Ncolsquare; i++) {
        u2[i]=tempor[i];
      }
      
      check_final_site=next;
    


/* Now, instead, we go first along nu, then mu: */
  
      next=site;
      for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
        if (smearing_mode==0) {
          tempor[i]=ufield[(nu*nsites+site)*Ncolsquare+i];
	}
	else {
          tempor[i]=smeared_ufield[(nu*nsites+site)*Ncolsquare+i];
	}
#else
        tempor[i]=ufield[(nu*nsites+site)*Ncolsquare+i];
#endif
      }
      next=neighbor_plus[nu*nsites+next];

      for (length=1; length<L_nu; length++) {
        for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
          if (smearing_mode==0) {
            u1[i]=ufield[(nu*nsites+next)*Ncolsquare+i];
	  }
	  else {
            u1[i]=smeared_ufield[(nu*nsites+next)*Ncolsquare+i];
	  }
#else
          u1[i]=ufield[(nu*nsites+next)*Ncolsquare+i];
#endif
        }
        mult_C_equals_AB(product,tempor,u1);
        for (i=0; i<Ncolsquare; i++) {
          tempor[i]=product[i];
        }
        next=neighbor_plus[nu*nsites+next];
      }
      
      for (length=0; length<L_mu; length++) {
        for (i=0; i<Ncolsquare; i++) {
#ifdef __wanna_smearing__
          if (smearing_mode==0) {
            u1[i]=ufield[(mu*nsites+next)*Ncolsquare+i];
	  }
	  else {
            u1[i]=smeared_ufield[(mu*nsites+next)*Ncolsquare+i];
	  }
#else
          u1[i]=ufield[(mu*nsites+next)*Ncolsquare+i];
#endif
        }
        mult_C_equals_AB(product,tempor,u1);
        for (i=0; i<Ncolsquare; i++) {
          tempor[i]=product[i];
        }
        next=neighbor_plus[mu*nsites+next];
      } 
  
if (next!=check_final_site) {printf("ERROR! The paths do not end at the same site: first path end = %d  second path end = %d, %s\n",next,check_final_site,(smearing_mode==0)?"not smeared":"smeared"); exit(0); /*do {;} while (1);*/}
  
      for (i=0; i<Ncolsquare; i++) {
        average+=u2[i]*conj(tempor[i]);
      }
    }
#ifdef __wanna_smearing__
    average/=(smearing_mode==0)?(Ncol*nsites*0.5*(dim-first_mu)*(dim-first_mu-1)):(Ncol*nsites*(dim-1));
#else
    average/=(Ncol*nsites*0.5*(dim-first_mu)*(dim-first_mu-1));
#endif
//     output_file << configuration << " " << L_mu << " " << L_nu << " " << real(average) << " " << imag(average) << endl;
    output_file << configuration+measurements_already_done << " " << L_mu << " " << L_nu << " " << real(average) << " " << imag(average) << endl;
  }

}
