void heat_bath_for_one_link(bool locked_mode, int site, int dir) {
  
  const double prefactor = (2.*beta)/Ncol;

  double raddet, alpha, cosq, biga, deltabar, phi, ics, azero, auno, adue, atre;

  int first_entry, second_entry, i, j;
  
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
      for (i=0;i<Ncol;i++) {
        for (j=0;j<Ncol;j++) {
          ak[i*Ncol+j]=dc(0.,0.);
          bigesse[i*Ncol+j] = conj(staple[j*Ncol+i]);
          oldu[i*Ncol+j] = ufield[(dir*nsites+site)*Ncolsquare+i*Ncol+j];
        }
        ak[i*(Ncol_plus_one)]=dc(1.,0.);
      }

      mult_C_equals_AB(tempor,oldu,bigesse);

      for (i=0;i<Ncolsquare;i++) {
        bigesse[i] = tempor[i];
      }

      rkappa[0] = 0.5*dc(
        real( bigesse[first_entry*Ncol+first_entry] + bigesse[second_entry*Ncol+second_entry] ),
        imag( bigesse[first_entry*Ncol+first_entry] - bigesse[second_entry*Ncol+second_entry] ) );
      rkappa[1] = 0.5*dc(
        real( bigesse[first_entry*Ncol+second_entry] - bigesse[second_entry*Ncol+first_entry] ),
        imag( bigesse[first_entry*Ncol+second_entry] + bigesse[second_entry*Ncol+first_entry] ) );
      rkappa[2] = -conj( rkappa[1] );
      rkappa[3] = conj( rkappa[0] );

      raddet=sqrt( norm(rkappa[0]) + norm(rkappa[1]) );

      alpha=1./(prefactor*raddet);

      do {
        cosq=cos(twopi*xx.randExc());
        cosq*=cosq;
        biga=-log(xx.randDblExc())*alpha*cosq;
        deltabar=-log(xx.randDblExc())*alpha +biga;
        ics=xx.randDblExc();
        ics*=ics;
        hb_tried+=1.;
      } while ( ics > (1. - 0.5*deltabar) );

      hb_accepted+=1.;
      azero=1.-deltabar;

      phi=twopi*xx.randExc();
      ics=1.-2.*xx.randDblExc();
      atre=ics;
      ics=sqrt( 1. - ics*ics );
      auno=ics*cos(phi);
      adue=ics*sin(phi);

      raddet=1./raddet;
      ics=sqrt( 1. - azero*azero )*raddet;
      azero*=raddet;

      auno*=ics;
      adue*=ics;
      atre*=ics;

      ak[first_entry*Ncol+first_entry] = dc(azero,atre)*rkappa[3] 
        - dc(adue,auno)*rkappa[2];
      ak[first_entry*Ncol+second_entry] = -dc(azero,atre)*rkappa[1]
        + dc(adue,auno)*rkappa[0];
      ak[second_entry*Ncol+first_entry] = dc(-adue,auno)*rkappa[3]
        + dc(-azero,atre)*rkappa[2];
      ak[second_entry*Ncol+second_entry] = dc(adue,-auno)*rkappa[1]
        + dc(azero,-atre)*rkappa[0];

      mult_C_equals_AB(newblock,ak,oldu);

      for (i=0;i<Ncolsquare;i++) {
        ufield[(dir*nsites+site)*Ncolsquare+i]=newblock[i];
      }
    }
  }

}
