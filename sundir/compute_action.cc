inline double compute_action(int site, int dir, dc *link
#ifdef __wanna_trace_deformation__
, dc *open_loop
#endif    
) {
        
  int i;
  dc temp_contribution;
  double tempaction;
  
#ifdef __wanna_trace_deformation__
#if Ncol>3
  int deformation_index;
#endif
#endif
  
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
  mult_C_equals_ABdagger(tempor,link,staple);
  tempaction=real(tempor[0]);
  for (i=Ncol_plus_one;i<Ncolsquare;i+=Ncol_plus_one) {
    tempaction+=real(tempor[i]);
  }
  tempaction*=Wilson_prefactor;
#ifdef __wanna_trace_deformation__
  if (is_compactified[dir]==true) { 
    mult_C_equals_AB(closed_loop,link,open_loop);

    temp_contribution=closed_loop[0];
    for (i=Ncol_plus_one;i<Ncolsquare;i+=Ncol_plus_one) {
      temp_contribution+=closed_loop[i];
    }
    tempaction+=double_trace_prefactor*alpha_coeff[0]*norm(temp_contribution);

#if Ncol>3
    for (i=0;i<Ncolsquare;i++) {
      tempor[i]=closed_loop[i];
    }
    for (deformation_index=1; deformation_index<how_many_trace_deformations; deformation_index++) {
      mult_C_equals_AB(tempor2,tempor,closed_loop);
      temp_contribution=tempor2[0];
      for (i=Ncol_plus_one;i<Ncolsquare;i+=Ncol_plus_one) {
        temp_contribution+=tempor2[i];
      }
      tempaction+=double_trace_prefactor*alpha_coeff[deformation_index]*norm(temp_contribution);
 
      if (deformation_index+1<how_many_trace_deformations)
        for (i=0;i<Ncolsquare;i++) {
          tempor[i]=tempor2[i];
        }
      }
#endif
  }
#endif

  return tempaction;

}
