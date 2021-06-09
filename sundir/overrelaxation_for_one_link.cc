void overrelaxation_for_one_link(bool locked_mode, int site, int dir) {

  double det_modulus;

  int first_entry, second_entry, i;
  
#ifdef __wanna_multilevel__
  if ((locked_link[dir*nsites+site]==false) || (locked_mode==false)) 
#endif
#ifdef __wanna_Jarzynski_SF__
  if (locked_link[dir*nsites+site]==false) 
#endif
  {
    pstaple(staple,site,dir,0,0);
    nstaple(staple,site,dir,0,0);
#ifdef __wanna_improvement__
    prectangularstaple(rectangularstaple,site,dir,0,0);
    nrectangularstaple(rectangularstaple,site,dir,0,0);
    for (i=0;i<Ncolsquare;i++) {
      staple[i]=
          plaquettefactor*staple[i]
        + rectanglefactor*rectangularstaple[i];
    }    
#endif
    for (first_entry=0;first_entry<Ncol-1;first_entry++)
    for (second_entry=first_entry+1;second_entry<Ncol;second_entry++) {
      for (i=0;i<Ncolsquare;i++) {
        ak[i]=dc(0.,0.);
        oldu[i] = ufield[(dir*nsites+site)*Ncolsquare+i];
      }
      for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
        ak[i]=dc(1.,0.);
      }

      mult_C_equals_ABdagger(bigesse,oldu,staple);
      
      rkappa[0] = 0.5*dc(
        real( bigesse[first_entry*Ncol+first_entry] + bigesse[second_entry*Ncol+second_entry] ),
        imag( bigesse[first_entry*Ncol+first_entry] - bigesse[second_entry*Ncol+second_entry] ) );
      rkappa[1] = 0.5*dc(
        real( bigesse[first_entry*Ncol+second_entry] - bigesse[second_entry*Ncol+first_entry] ),
        imag( bigesse[first_entry*Ncol+second_entry] + bigesse[second_entry*Ncol+first_entry] ) );
      rkappa[2] = -conj( rkappa[1] );
      rkappa[3] = conj( rkappa[0] );

      det_modulus= norm(rkappa[0]) + norm(rkappa[1]) ;
      det_modulus=1./det_modulus;

      ak[first_entry*Ncol+first_entry] = rkappa[3]*rkappa[3]+rkappa[1]*rkappa[2];
      ak[first_entry*Ncol+second_entry] = -rkappa[3]*rkappa[1]-rkappa[1]*rkappa[0];
      ak[second_entry*Ncol+first_entry] = -conj(ak[first_entry*Ncol+second_entry]);
      ak[second_entry*Ncol+second_entry] = conj(ak[first_entry*Ncol+first_entry]);

      ak[first_entry*Ncol+first_entry]*=det_modulus;
      ak[first_entry*Ncol+second_entry]*=det_modulus;
      ak[second_entry*Ncol+first_entry]*=det_modulus;
      ak[second_entry*Ncol+second_entry]*=det_modulus;

      mult_C_equals_AB(newblock,ak,oldu);

      for (i=0;i<Ncolsquare;i++) {
        ufield[(dir*nsites+site)*Ncolsquare+i]=newblock[i];
      }
    }
  }

}
