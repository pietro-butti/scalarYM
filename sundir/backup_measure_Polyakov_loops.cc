dc measure_Polyakov_loops(int configuration, dc **zero_mom_loops_now, dc **zero_mom_corr_now
#ifdef __wanna_nonzero_mom_correlators__
, dc *nonzero_mom_corr_now
#endif
) {

  int i, t, x, y, steps_along_x, spatial_site;
#ifdef __wanna_nonzero_mom_correlators__
  spatial_coordinate;
#endif
  int transverse_volume=ny;
#if dim==4
  int z;
  transverse_volume=ny*nz;
#endif

  dc average = dc(0.,0.);
  dc trace = dc(0.,0.);
  dc temporary_loop = dc(0.,0.);
  dc trace_of_one_loop= dc(0.,0.);
  dc f= dc(0.,0.);
  dc fbar= dc(0.,0.);
  dc temporary_trace[how_many_irreps];
  dc aux;
#ifdef __wanna_nonzero_mom_correlators__
  dc *loop_through_x=new dc [spatial_volume];
  dc nonzero_mom_Polyakov_correlator;
  int dist;
#endif
  
#if Ncol>3
  dc fundrep_eigenvalues[Ncol];
  dc fundrep_eigenvalue_squares[Ncol];
  dc fundrep_eigenvalue_cubes[Ncol];
  double fundrep_phases[Ncol];
  for (i=0;i<Ncol;i++) {
    fundrep_eigenvalues[i]=dc(0.,0.);
    fundrep_phases[i]=0.;
  }
  struct complex_double b[Ncol], DUMMY[1][1], WORK[2*Ncol];
  
  double AT[2*Ncol*Ncol];
  int ok, c1, c2, c3;
  char c4;
 
#endif

  for (x=0;x<nx;x++) {
    for (y=0;y<ny;y++)
#if dim==4
    for (z=0;z<nz;z++)
#endif
    {    
  
      spatial_site=y+ny*x;
#if dim==4
      spatial_site*=nz;
      spatial_site+=z;
#endif
#ifdef __wanna_nonzero_mom_correlators__
      spatial_coordinate=spatial_site;
#endif
      for (i=0;i<Ncolsquare;i++) {
        tempor[i]=ufield[spatial_site*Ncolsquare+i];
#ifdef __debugging_mode__
tempor[i]=dc(0.,0.);
#endif
      }
#ifdef __debugging_mode__
for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {tempor[i]=dc(1.,0.);}
#endif

      for (t=1;t<nt-1;t++) {
        spatial_site=y+ny*(x+nx*t);
#if dim==4
        spatial_site*=nz;
        spatial_site+=z;
#endif
        for (i=0;i<Ncolsquare;i++) {
          tempor2[i]=ufield[spatial_site*Ncolsquare+i];
#ifdef __debugging_mode__
tempor2[i]=dc(0.,0.);
#endif
        }
#ifdef __debugging_mode__
for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {tempor2[i]=dc(1.,0.);}
#endif
        mult_C_equals_AB(product,tempor,tempor2);
        for (i=0;i<Ncolsquare;i++) {
          tempor[i]=product[i];
        }
      }

      spatial_site=y+ny*(x+nx*(nt-1));
#if dim==4
      spatial_site*=nz;
      spatial_site+=z;
#endif
      for (i=0;i<Ncolsquare;i++) {
        tempor2[i]=ufield[spatial_site*Ncolsquare+i];
#ifdef __debugging_mode__
tempor2[i]=dc(0.,0.);
#endif
      }
#ifdef __debugging_mode__
for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {tempor2[i]=dc(1.,0.);}
#endif
      mult_C_equals_AB(product,tempor,tempor2);
      


#if Ncol>3
      for (i=0; i<Ncolsquare; i++) {
        AT[2*i]=real(product[i]);
        AT[2*i+1]=imag(product[i]);
      }
      c1=Ncol;
      c2=2*Ncol;
      c3=1;
      c4='N';

      zgeev_(&c4, &c4,&c1, AT, &c1, b, DUMMY, &c3, DUMMY, &c3, WORK, &c2, WORK, &ok);

// dc determinant=dc(1.,0.);

      if (ok==0) {
        for (i=0; i<Ncol; i++) {
          fundrep_eigenvalues[i]=dc(b[i].re, b[i].im);
          fundrep_eigenvalue_squares[i]=
            fundrep_eigenvalues[i]*fundrep_eigenvalues[i];
          fundrep_eigenvalue_cubes[i]=
            fundrep_eigenvalue_squares[i]*fundrep_eigenvalues[i];
          fundrep_phases[i]=imag(log(fundrep_eigenvalues[i]));
// cout << real(fundrep_eigenvalues[i]) << " " << imag(fundrep_eigenvalues[i]) << endl;
// cout << fundrep_phases[i] << endl;
// determinant*=fundrep_eigenvalues[i];
        }
// cout << "Determinant = " << determinante << endl;
// exit(0);
      }
      else { printf("Error in the computation of the Polyakov line eigenvalues\n"); exit(0); }
#endif

      trace_of_one_loop=product[0];
      for (i=Ncol_plus_one;i<Ncolsquare;i+=Ncol_plus_one) {
        trace_of_one_loop+=product[i];
      }
//       trace_of_one_loop/=Ncol;

// Fundamental representation:
      temporary_trace[0]=trace_of_one_loop;
#ifdef __wanna_nonzero_mom_correlators__
      loop_through_x[spatial_coordinate]=trace_of_one_loop;
#endif
      
      zero_mom_loops_now[0][x]+=trace_of_one_loop;
      trace+=trace_of_one_loop;
            
      f=trace_of_one_loop;
      fbar=conj(f);
      
#if Ncol==2
#if how_many_irreps>1
// Representation 3:
      temporary_trace[1]=f*fbar-1.;
      zero_mom_loops_now[1][x]+=temporary_trace[1];
     
#if how_many_irreps>2       
// Representation 4:
      temporary_trace[2]=f*(temporary_trace[1]-1.);
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 5:
      temporary_trace[3]=f*temporary_trace[2]-temporary_trace[1];
      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
// Representation 6:
      temporary_trace[4]=f*temporary_trace[3]-temporary_trace[2];
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
 #if how_many_irreps>5    
// Representation 7:
      temporary_trace[5]=f*temporary_trace[4]-temporary_trace[3];
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
// Representation 8:
      temporary_trace[6]=f*temporary_trace[5]-temporary_trace[4];
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
// Representation 9:
      temporary_trace[7]=f*temporary_trace[6]-temporary_trace[5];
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 10:
      temporary_trace[8]=f*temporary_trace[7]-temporary_trace[6];
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>8
// Representation 11:
      temporary_trace[9]=f*temporary_trace[8]-temporary_trace[7];
      zero_mom_loops_now[9][x]+=temporary_trace[9];
     
#if how_many_irreps>10 
// Representation 12:
      temporary_trace[10]=f*temporary_trace[9]-temporary_trace[8];
      zero_mom_loops_now[10][x]+=temporary_trace[10];
    
#if how_many_irreps>11 
// Representation 13:
      temporary_trace[11]=f*temporary_trace[10]-temporary_trace[9];
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      
      
#if Ncol==3
#if how_many_irreps>1
// Representation 6:
      temporary_trace[1]=f*f - fbar;
      zero_mom_loops_now[1][x]+=temporary_trace[1];
     
#if how_many_irreps>2       
// Representation 8:
      temporary_trace[2]=f*fbar - 1.;
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 10:
      temporary_trace[3]=f*temporary_trace[1]-temporary_trace[2];
      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
// Representation 15:
      temporary_trace[4]=fbar*temporary_trace[1]-f;
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
 #if how_many_irreps>5    
// Representation 15':
      temporary_trace[5]=f*temporary_trace[3]-temporary_trace[4];
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
/**********************************************************/
/* Now we already calculate the 24, to use it for the 21: */
      aux=fbar*temporary_trace[3]-temporary_trace[1];
/**********************************************************/

// Representation 21:
      temporary_trace[6]=f*temporary_trace[5]-aux;
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
// Representation 24:
      temporary_trace[7]=aux;
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 27:
      temporary_trace[8]=temporary_trace[1]*conj(temporary_trace[1])
         -temporary_trace[2]-1.;
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>8
/**********************************************************/
/* Now we already calculate the 35, to use it for the 28: */
      aux=temporary_trace[3]*temporary_trace[2]
         -temporary_trace[8]
         -temporary_trace[3]
         -temporary_trace[2];
/**********************************************************/
     
// Representation 28:
      temporary_trace[9]=f*temporary_trace[6]-aux;
      zero_mom_loops_now[9][x]+=temporary_trace[9];
     
#if how_many_irreps>10 
// Representation 35:
      temporary_trace[10]=aux;
      zero_mom_loops_now[10][x]+=temporary_trace[10];
    
#if how_many_irreps>11 
/**********************************************************/
/* Now we calculate the 48, to use it for the 36: */
      aux=fbar*temporary_trace[6]-temporary_trace[5];
/**********************************************************/
   
// Representation 36:
      temporary_trace[11]=f*temporary_trace[9]-aux;
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      
#if Ncol==4

#if how_many_irreps>1
// Representation 6:
      temporary_trace[1]=dc(
        2.*( cos(fundrep_phases[0]+fundrep_phases[1]) 
           + cos(fundrep_phases[0]+fundrep_phases[2]) 
           + cos(fundrep_phases[1]+fundrep_phases[2]) 
           ),
        0.);
      zero_mom_loops_now[1][x]+=temporary_trace[1];
      
#if how_many_irreps>2
// Representation 10:
      temporary_trace[2]=f*f-temporary_trace[1];
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 15:
      temporary_trace[3]=f*fbar-1.;		
      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
// Representation 20:
      temporary_trace[4]=f*temporary_trace[1]-fbar;
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
#if how_many_irreps>5
// Representation 20':
      temporary_trace[5]=dc(
        2.*( 1. + cos(fundrep_phases[0]-fundrep_phases[1]) 
                + cos(fundrep_phases[0]-fundrep_phases[2])  
                + cos(fundrep_phases[0]-fundrep_phases[3])
                + cos(fundrep_phases[1]-fundrep_phases[2]) 
                + cos(fundrep_phases[1]-fundrep_phases[3])
                + cos(fundrep_phases[2]-fundrep_phases[3])
                + cos( 2.*(fundrep_phases[0]+fundrep_phases[1]) )
                + cos( 2.*(fundrep_phases[0]+fundrep_phases[2]) )
                + cos( 2.*(fundrep_phases[1]+fundrep_phases[2]) ) ),
                0.);
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
// Representation 20'':
      temporary_trace[6]=temporary_trace[2]*f-temporary_trace[4];
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
/**********************************************************/
/* Now we already calculate the 45, to use it for the 35: */
      aux=temporary_trace[1]*temporary_trace[2]-temporary_trace[3];
/**********************************************************/
      
// Representation 35:
      temporary_trace[7]=f*temporary_trace[6]-aux;
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 36:
      temporary_trace[8]=f*temporary_trace[3]-f-conj(temporary_trace[4]);
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>9
// Representation 45:
      temporary_trace[9]=aux;
      zero_mom_loops_now[9][x]+=temporary_trace[9];
      
#if how_many_irreps>10
// Representation 50:
      temporary_trace[10]=dc( 
        2.*(
               cos(2.*fundrep_phases[0])
             + cos(2.*fundrep_phases[1])
             + cos(2.*fundrep_phases[2])
             + cos(2.*fundrep_phases[3])
             + 2.*(
                    cos(fundrep_phases[0]+fundrep_phases[1])
                  + cos(fundrep_phases[0]+fundrep_phases[2])
                  + cos(fundrep_phases[1]+fundrep_phases[2])
                  )
             + cos( 3.*( fundrep_phases[0]+fundrep_phases[1]) )
             + cos( 3.*( fundrep_phases[0]+fundrep_phases[2]) )
             + cos( 3.*( fundrep_phases[1]+fundrep_phases[2]) )
             
             + cos(2.*fundrep_phases[0]+fundrep_phases[1]-fundrep_phases[2])
             + cos(2.*fundrep_phases[1]+fundrep_phases[2]-fundrep_phases[0])
             + cos(2.*fundrep_phases[2]+fundrep_phases[0]-fundrep_phases[1])
             
             + cos(2.*fundrep_phases[0]-fundrep_phases[1]+fundrep_phases[2])
             + cos(2.*fundrep_phases[1]-fundrep_phases[2]+fundrep_phases[0])
             + cos(2.*fundrep_phases[2]-fundrep_phases[0]+fundrep_phases[1])
             
             + cos(3.*fundrep_phases[0]+2.*fundrep_phases[1]+fundrep_phases[2])
             + cos(3.*fundrep_phases[1]+2.*fundrep_phases[2]+fundrep_phases[0])
             + cos(3.*fundrep_phases[2]+2.*fundrep_phases[0]+fundrep_phases[1])
             
             + cos(3.*fundrep_phases[0]+fundrep_phases[1]+2.*fundrep_phases[2])
             + cos(3.*fundrep_phases[1]+fundrep_phases[2]+2.*fundrep_phases[0])
             + cos(3.*fundrep_phases[2]+fundrep_phases[0]+2.*fundrep_phases[1])
           ),
        0.);
      zero_mom_loops_now[10][x]+=temporary_trace[10];
      
#if how_many_irreps>11
/**********************************************************/
/* Now we already calculate the 60 and the 84, to use them for the 56; */
/* the 60 is: */
      aux=f*temporary_trace[5]-conj(temporary_trace[4]);

/* From this, we can obtain the 84: */
      aux= temporary_trace[2]*temporary_trace[4]
          -temporary_trace[8]-conj(temporary_trace[4])
          -aux;
/**********************************************************/
      
// Representation 56:
      temporary_trace[11]=f*temporary_trace[7]-aux;
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      
#if Ncol==5

#if how_many_irreps>1
/**********************************************************/
/* Now we already calculate the 15, to use it for the 10: */
      aux=fundrep_eigenvalue_squares[0]
        + fundrep_eigenvalue_squares[1]
        + fundrep_eigenvalue_squares[2]
        + fundrep_eigenvalue_squares[3]
        + fundrep_eigenvalue_squares[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4];
/**********************************************************/

// Representation 10:
      temporary_trace[1]=f*f-aux;
      zero_mom_loops_now[1][x]+=temporary_trace[1];
      
#if how_many_irreps>2
// Representation 15:
      temporary_trace[2]=f*f-temporary_trace[1];
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 24:
      temporary_trace[3]=f*fbar-1.;
      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
/**********************************************************/
/* Now we already calculate the 40, to use it for the 35: */
      aux=f*temporary_trace[1]-conj(temporary_trace[1]);
/**********************************************************/

// Representation 35:
      temporary_trace[4]=f*temporary_trace[2]-aux;
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
#if how_many_irreps>5
// Representation 40:
      temporary_trace[5]=aux;
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
// Representation 45:
      temporary_trace[6]=f*conj(temporary_trace[1])-fbar;
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
// Representation 50:
      temporary_trace[7]=temporary_trace[1]*temporary_trace[1]-temporary_trace[6]-fbar;
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 70:
      temporary_trace[8]=f*temporary_trace[3]-conj(temporary_trace[6])-f;
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>9
/**********************************************************/
/* Now we already calculate the 105, to use it for the 70': */
      aux=temporary_trace[1]*temporary_trace[2]-temporary_trace[6];
/**********************************************************/
      
// Representation 70':
      temporary_trace[9]=f*temporary_trace[4]-aux;
      zero_mom_loops_now[9][x]+=temporary_trace[9];
      
#if how_many_irreps>10
// Representation 75:
      temporary_trace[10]=conj(temporary_trace[1])*temporary_trace[1]
           -temporary_trace[3]-1.;
      zero_mom_loops_now[10][x]+=temporary_trace[10];
      
#if how_many_irreps>11
// Representation 105:
      temporary_trace[11]=aux;
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      
#if Ncol==6

#if how_many_irreps>1
// Representation 15:
      temporary_trace[1]= fundrep_eigenvalues[0]*fundrep_eigenvalues[1]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[5];
      zero_mom_loops_now[1][x]+=temporary_trace[1];
      
#if how_many_irreps>2
// Representation 20:
      temporary_trace[2]=2.*(
          cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[2])
            + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[3])
            + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[4])
            + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[3])
            + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[4])            
            + cos(fundrep_phases[0]+fundrep_phases[3]+fundrep_phases[4])
            + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[3])
            + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[4])
            + cos(fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[4])
            + cos(fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[4])
        );
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 21:
      temporary_trace[3]=f*f-temporary_trace[1];	
      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
// Representation 35:
      temporary_trace[4]=f*fbar-1.;
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
#if how_many_irreps>5
/**********************************************************/
/* Now we already calculate the 70, to use it for the 56: */
      aux=f*temporary_trace[1]-temporary_trace[2];
/**********************************************************/
      
// Representation 56:
      temporary_trace[5]=f*temporary_trace[3]-aux;
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
// Representation 70:
      temporary_trace[6]=aux;
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
// Representation 84:
      temporary_trace[7]=f*conj(temporary_trace[1])-fbar;
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 105 - the one with canonical label (1,0,1,0,0):
      temporary_trace[8]=f*temporary_trace[2]-conj(temporary_trace[1]);
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>9
// Representation 105' - the one with canonical label (0,2,0,0,0):
      temporary_trace[9]=temporary_trace[1]*temporary_trace[1]-temporary_trace[8]-conj(temporary_trace[1]);
      zero_mom_loops_now[9][x]+=temporary_trace[9];
      
#if how_many_irreps>10
// Representation 120:
      temporary_trace[10]=f*temporary_trace[4]-conj(temporary_trace[7])-f;
      zero_mom_loops_now[10][x]+=temporary_trace[10];
      
#if how_many_irreps>11
/**********************************************************/
/* Now we calculate the 210, to use it for the 126: */
      aux=f*temporary_trace[6]-temporary_trace[8]-temporary_trace[9];
/**********************************************************/
      
// Representation 126:
      temporary_trace[11]=f*temporary_trace[5]-aux;
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif



#if Ncol==7

#if how_many_irreps>1
/**********************************************************/
/* Now we calculate a quantity, which appears also later: */
      aux=fundrep_eigenvalues[0]*(
                 fundrep_eigenvalues[1]
                +fundrep_eigenvalues[2]
                +fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
         +fundrep_eigenvalues[1]*(
                 fundrep_eigenvalues[2]
                +fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
         +fundrep_eigenvalues[2]*(
                 fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
         +fundrep_eigenvalues[3]*(
                 fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
         +fundrep_eigenvalues[4]*(
                 fundrep_eigenvalues[5]);
/**********************************************************/

// Representation 21:
      temporary_trace[1]= aux 
                        + fundrep_eigenvalues[0]
                        + fundrep_eigenvalues[1]
                        + fundrep_eigenvalues[2]
                        + fundrep_eigenvalues[3]
                        + fundrep_eigenvalues[4]
                        + fundrep_eigenvalues[5];
			
      zero_mom_loops_now[1][x]+=temporary_trace[1];
			
#if how_many_irreps>2
// Representation 28:
      temporary_trace[2]=f*f-temporary_trace[1];
      zero_mom_loops_now[2][x]+=temporary_trace[2];

#if how_many_irreps>3
// Representation 35:
      temporary_trace[3]= aux               
         +fundrep_eigenvalues[0]*(
             fundrep_eigenvalues[1]*(
                 fundrep_eigenvalues[2]
                +fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[2]*(
                 fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[3]*(
                 fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[4]*(
                 fundrep_eigenvalues[5]))
                 
         +fundrep_eigenvalues[1]*(
             fundrep_eigenvalues[2]*(
                 fundrep_eigenvalues[3]
                +fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[3]*(
                 fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[4]*(
                 fundrep_eigenvalues[5]))
                 
         +fundrep_eigenvalues[2]*(
             fundrep_eigenvalues[3]*(
                 fundrep_eigenvalues[4]
                +fundrep_eigenvalues[5])
            +fundrep_eigenvalues[4]*(
                 fundrep_eigenvalues[5]))
                 
         +fundrep_eigenvalues[3]*(
             fundrep_eigenvalues[4]*(
                 fundrep_eigenvalues[5]));
	 
      zero_mom_loops_now[3][x]+=temporary_trace[3];

#if how_many_irreps>4
// Representation 48:
      temporary_trace[4]=f*fbar-1.;
      zero_mom_loops_now[4][x]+=temporary_trace[4];

#if how_many_irreps>5
/**********************************************************/
/* Now we calculate the 112, to use it for the 84: */
      aux=f*temporary_trace[1]-temporary_trace[3];
/**********************************************************/
    
// Representation 84:
      temporary_trace[5]=temporary_trace[2]*f-aux;
      zero_mom_loops_now[5][x]+=temporary_trace[5];

#if how_many_irreps>6
// Representation 112:
      temporary_trace[6]=aux;
      zero_mom_loops_now[6][x]+=temporary_trace[6];

#if how_many_irreps>7
// Representation 140:
      temporary_trace[7]=conj(temporary_trace[1])*f-fbar;
      zero_mom_loops_now[7][x]+=temporary_trace[7];

#if how_many_irreps>8
// Representation 189:
      temporary_trace[8]=temporary_trace[4]*f-conj(temporary_trace[7])-fbar;
      zero_mom_loops_now[8][x]+=temporary_trace[8];

#if how_many_irreps>9
// Representation 196:
      temporary_trace[9]= f*temporary_trace[6]
                         -temporary_trace[1]*temporary_trace[2];
      zero_mom_loops_now[9][x]+=temporary_trace[9];

#if how_many_irreps>10
// Representation 210:
      temporary_trace[10]=temporary_trace[3]*f-conj(temporary_trace[3]);
      zero_mom_loops_now[10][x]+=temporary_trace[10];

#if how_many_irreps>11
// Representation 210':
      temporary_trace[11]=f*temporary_trace[5]
                         -temporary_trace[1]*temporary_trace[2]
                         +temporary_trace[10];
      zero_mom_loops_now[11][x]+=temporary_trace[11];

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif


      
#if Ncol==8
      
#if how_many_irreps>1
/**********************************************************/
/* Now we already calculate the 36, to use it for the 28: */
      aux=  fundrep_eigenvalue_squares[0]
        + fundrep_eigenvalue_squares[1]
        + fundrep_eigenvalue_squares[2]
        + fundrep_eigenvalue_squares[3]
        + fundrep_eigenvalue_squares[4]
        + fundrep_eigenvalue_squares[5]
        + fundrep_eigenvalue_squares[6]
        + fundrep_eigenvalue_squares[7]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        + fundrep_eigenvalues[6]*fundrep_eigenvalues[7];
/**********************************************************/
	
// Representation 28:
      temporary_trace[1]=f*f-aux;
      zero_mom_loops_now[1][x]+=temporary_trace[1];
      
#if how_many_irreps>2
// Representation 36:
      temporary_trace[2]=aux;
      zero_mom_loops_now[2][x]+=temporary_trace[2];
      
#if how_many_irreps>3
// Representation 56:
      temporary_trace[3]=
          fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[2]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[1]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[2]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[3]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[4]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[0]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7]

        
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]*fundrep_eigenvalues[3]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[2]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[3]*fundrep_eigenvalues[7]
                              
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[4]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[1]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7]
        
        
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]*fundrep_eigenvalues[4]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[3]*fundrep_eigenvalues[7]
                              
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[4]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[2]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7]

        
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4]*fundrep_eigenvalues[5]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[4]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[3]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7]

        
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[5]*fundrep_eigenvalues[6]
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[5]*fundrep_eigenvalues[7]
        
        + fundrep_eigenvalues[4]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7]
        
        
        + fundrep_eigenvalues[5]*fundrep_eigenvalues[6]*fundrep_eigenvalues[7];

      zero_mom_loops_now[3][x]+=temporary_trace[3];
      
#if how_many_irreps>4
// Representation 63:
      temporary_trace[4]=f*fbar-1.;
      zero_mom_loops_now[4][x]+=temporary_trace[4];
      
#if how_many_irreps>5
// Representation 70:
      temporary_trace[5]=
        dc(2.*(
                  cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[3])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[4])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[6])
                
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[4])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[6])
                
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[4]+fundrep_phases[6])
                
                + cos(fundrep_phases[0]+fundrep_phases[1]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[4])
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[6])
                
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[4]+fundrep_phases[6])
               
                + cos(fundrep_phases[0]+fundrep_phases[2]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[0]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[0]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[6])
               
                + cos(fundrep_phases[0]+fundrep_phases[3]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[0]+fundrep_phases[4]+fundrep_phases[5]+fundrep_phases[6])
                
                
                
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[4])
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[5])
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[6])
                
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[4]+fundrep_phases[6])
               
                + cos(fundrep_phases[1]+fundrep_phases[2]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[6])
               
                + cos(fundrep_phases[1]+fundrep_phases[3]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[1]+fundrep_phases[4]+fundrep_phases[5]+fundrep_phases[6])
                
                
                
                + cos(fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[5])
                + cos(fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[6])
               
                + cos(fundrep_phases[2]+fundrep_phases[3]+fundrep_phases[5]+fundrep_phases[6])
                
                
                + cos(fundrep_phases[2]+fundrep_phases[4]+fundrep_phases[5]+fundrep_phases[6])
                
                
                
                + cos(fundrep_phases[3]+fundrep_phases[4]+fundrep_phases[5]+fundrep_phases[6])
              ),
           0.);
      zero_mom_loops_now[5][x]+=temporary_trace[5];
      
#if how_many_irreps>6
/**********************************************************/
/* Now we calculate the 168, to use it for the 120: */
      aux=f*temporary_trace[1]-temporary_trace[3];
/**********************************************************/
      
// Representation 120:
      temporary_trace[6]=f*temporary_trace[2]-aux;
      zero_mom_loops_now[6][x]+=temporary_trace[6];
      
#if how_many_irreps>7
// Representation 168:
      temporary_trace[7]=aux;
      zero_mom_loops_now[7][x]+=temporary_trace[7];
      
#if how_many_irreps>8
// Representation 216:
      temporary_trace[8]=f*conj(temporary_trace[1])-fbar;
      zero_mom_loops_now[8][x]+=temporary_trace[8];
      
#if how_many_irreps>9
// Representation 280:
      temporary_trace[9]=fbar*temporary_trace[2]-f;
      zero_mom_loops_now[9][x]+=temporary_trace[9];
      
#if how_many_irreps>10
// Representation 330:
      temporary_trace[10]=f*(temporary_trace[6]+temporary_trace[3])
          - temporary_trace[2]*temporary_trace[1] 
          - temporary_trace[5];
      zero_mom_loops_now[10][x]+=temporary_trace[10];
      
#if how_many_irreps>11
// Representation 336:
      temporary_trace[11]=temporary_trace[1]*temporary_trace[1] 
          -f*temporary_trace[3];
      zero_mom_loops_now[11][x]+=temporary_trace[11];
      
      
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif



    }
    
    zero_mom_loops_now[0][x]/=(transverse_volume*Ncol);
    
    
#if Ncol==2

#if how_many_irreps>1
// Representation 3:
    zero_mom_loops_now[1][x]/=(transverse_volume*3);
      
#if how_many_irreps>2
// Representation 4:
    zero_mom_loops_now[2][x]/=(transverse_volume*4);
      
#if how_many_irreps>3
// Representation 5:
    zero_mom_loops_now[3][x]/=(transverse_volume*5);
      
#if how_many_irreps>4
// Representation 6:
    zero_mom_loops_now[4][x]/=(transverse_volume*6);
      
#if how_many_irreps>5
// Representation 7:
    zero_mom_loops_now[5][x]/=(transverse_volume*7);
      
#if how_many_irreps>6
// Representation 8:
    zero_mom_loops_now[6][x]/=(transverse_volume*8);
      
#if how_many_irreps>7
// Representation 9:
    zero_mom_loops_now[7][x]/=(transverse_volume*9);
      
#if how_many_irreps>8
// Representation 10:
    zero_mom_loops_now[8][x]/=(transverse_volume*10);
      
#if how_many_irreps>9
// Representation 11:
    zero_mom_loops_now[9][x]/=(transverse_volume*11);
      
#if how_many_irreps>10
// Representation 12:
    zero_mom_loops_now[10][x]/=(transverse_volume*12);
      
#if how_many_irreps>11
// Representation 13:
    zero_mom_loops_now[11][x]/=(transverse_volume*13);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==3

#if how_many_irreps>1
// Representation 6:
    zero_mom_loops_now[1][x]/=(transverse_volume*6);
      
#if how_many_irreps>2
// Representation 8:
    zero_mom_loops_now[2][x]/=(transverse_volume*8);
      
#if how_many_irreps>3
// Representation 10:
    zero_mom_loops_now[3][x]/=(transverse_volume*10);
      
#if how_many_irreps>4
// Representation 15:
    zero_mom_loops_now[4][x]/=(transverse_volume*15);
      
#if how_many_irreps>5
// Representation 15':
    zero_mom_loops_now[5][x]/=(transverse_volume*15);
      
#if how_many_irreps>6
// Representation 21:
    zero_mom_loops_now[6][x]/=(transverse_volume*21);
      
#if how_many_irreps>7
// Representation 24:
    zero_mom_loops_now[7][x]/=(transverse_volume*24);
      
#if how_many_irreps>8
// Representation 27:
    zero_mom_loops_now[8][x]/=(transverse_volume*27);
      
#if how_many_irreps>9
// Representation 28:
    zero_mom_loops_now[9][x]/=(transverse_volume*28);
      
#if how_many_irreps>10
// Representation 35:
    zero_mom_loops_now[10][x]/=(transverse_volume*35);
      
#if how_many_irreps>11
// Representation 36:
    zero_mom_loops_now[11][x]/=(transverse_volume*36);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==4

#if how_many_irreps>1
// Representation 6:
    zero_mom_loops_now[1][x]/=(transverse_volume*6);
      
#if how_many_irreps>2
// Representation 10:
    zero_mom_loops_now[2][x]/=(transverse_volume*10);
      
#if how_many_irreps>3
// Representation 15:
    zero_mom_loops_now[3][x]/=(transverse_volume*15);
      
#if how_many_irreps>4
// Representation 20:
    zero_mom_loops_now[4][x]/=(transverse_volume*20);
      
#if how_many_irreps>5
// Representation 20':
    zero_mom_loops_now[5][x]/=(transverse_volume*20);
      
#if how_many_irreps>6
// Representation 20'':
    zero_mom_loops_now[6][x]/=(transverse_volume*20);
      
#if how_many_irreps>7
// Representation 35:
    zero_mom_loops_now[7][x]/=(transverse_volume*35);
      
#if how_many_irreps>8
// Representation 36:
    zero_mom_loops_now[8][x]/=(transverse_volume*36);
      
#if how_many_irreps>9
// Representation 45:
    zero_mom_loops_now[9][x]/=(transverse_volume*45);
      
#if how_many_irreps>10
// Representation 50:
    zero_mom_loops_now[10][x]/=(transverse_volume*50);
      
#if how_many_irreps>11
// Representation 56:
    zero_mom_loops_now[11][x]/=(transverse_volume*56);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==5

#if how_many_irreps>1
// Representation 10:
    zero_mom_loops_now[1][x]/=(transverse_volume*10);
      
#if how_many_irreps>2
// Representation 15:
    zero_mom_loops_now[2][x]/=(transverse_volume*15);
      
#if how_many_irreps>3
// Representation 24:
    zero_mom_loops_now[3][x]/=(transverse_volume*24);
      
#if how_many_irreps>4
// Representation 35:
    zero_mom_loops_now[4][x]/=(transverse_volume*35);
      
#if how_many_irreps>5
// Representation 40:
    zero_mom_loops_now[5][x]/=(transverse_volume*40);
      
#if how_many_irreps>6
// Representation 45:
    zero_mom_loops_now[6][x]/=(transverse_volume*45);
      
#if how_many_irreps>7
// Representation 50:
    zero_mom_loops_now[7][x]/=(transverse_volume*50);
      
#if how_many_irreps>8
// Representation 70:
    zero_mom_loops_now[8][x]/=(transverse_volume*70);
      
#if how_many_irreps>9
// Representation 70':
    zero_mom_loops_now[9][x]/=(transverse_volume*70);
      
#if how_many_irreps>10
// Representation 75:
    zero_mom_loops_now[10][x]/=(transverse_volume*75);
      
#if how_many_irreps>11
// Representation 105:
    zero_mom_loops_now[11][x]/=(transverse_volume*105);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==6

#if how_many_irreps>1
// Representation 15:
    zero_mom_loops_now[1][x]/=(transverse_volume*15);
      
#if how_many_irreps>2
// Representation 20:
    zero_mom_loops_now[2][x]/=(transverse_volume*20);
      
#if how_many_irreps>3
// Representation 21:
    zero_mom_loops_now[3][x]/=(transverse_volume*21);
      
#if how_many_irreps>4
// Representation 35:
    zero_mom_loops_now[4][x]/=(transverse_volume*35);
      
#if how_many_irreps>5
// Representation 56:
    zero_mom_loops_now[5][x]/=(transverse_volume*56);
      
#if how_many_irreps>6
// Representation 70:
    zero_mom_loops_now[6][x]/=(transverse_volume*70);
      
#if how_many_irreps>7
// Representation 84:
    zero_mom_loops_now[7][x]/=(transverse_volume*84);
      
#if how_many_irreps>8
// Representation 105:
    zero_mom_loops_now[8][x]/=(transverse_volume*105);
      
#if how_many_irreps>9
// Representation 105':
    zero_mom_loops_now[9][x]/=(transverse_volume*105);
      
#if how_many_irreps>10
// Representation 120:
    zero_mom_loops_now[10][x]/=(transverse_volume*120);
      
#if how_many_irreps>11
// Representation 126:
    zero_mom_loops_now[11][x]/=(transverse_volume*126);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==7

#if how_many_irreps>1
// Representation 21:
    zero_mom_loops_now[1][x]/=(transverse_volume*21);
      
#if how_many_irreps>2
// Representation 28:
    zero_mom_loops_now[2][x]/=(transverse_volume*28);
      
#if how_many_irreps>3
// Representation 35:
    zero_mom_loops_now[3][x]/=(transverse_volume*35);
      
#if how_many_irreps>4
// Representation 48:
    zero_mom_loops_now[4][x]/=(transverse_volume*48);
      
#if how_many_irreps>5
// Representation 84:
    zero_mom_loops_now[5][x]/=(transverse_volume*84);
      
#if how_many_irreps>6
// Representation 112:
    zero_mom_loops_now[6][x]/=(transverse_volume*112);

#if how_many_irreps>7     
// Representation 140:
    zero_mom_loops_now[7][x]/=(transverse_volume*140);
      
#if how_many_irreps>8
// Representation 189:
    zero_mom_loops_now[8][x]/=(transverse_volume*189);
      
#if how_many_irreps>9
// Representation 196:
    zero_mom_loops_now[9][x]/=(transverse_volume*196);
      
#if how_many_irreps>10
// Representation 210:
    zero_mom_loops_now[10][x]/=(transverse_volume*210);
      
#if how_many_irreps>11
// Representation 210':
    zero_mom_loops_now[11][x]/=(transverse_volume*210);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    
    
#if Ncol==8

#if how_many_irreps>1
// Representation 28:
    zero_mom_loops_now[1][x]/=(transverse_volume*28);
      
#if how_many_irreps>2
// Representation 36:
    zero_mom_loops_now[2][x]/=(transverse_volume*36);
      
#if how_many_irreps>3
// Representation 56:
    zero_mom_loops_now[3][x]/=(transverse_volume*56);
      
#if how_many_irreps>4
// Representation 63:
    zero_mom_loops_now[4][x]/=(transverse_volume*63);
      
#if how_many_irreps>5
// Representation 70:
    zero_mom_loops_now[5][x]/=(transverse_volume*70);
      
#if how_many_irreps>6
// Representation 120:
    zero_mom_loops_now[6][x]/=(transverse_volume*120);

#if how_many_irreps>7     
// Representation 168:
    zero_mom_loops_now[7][x]/=(transverse_volume*168);
      
#if how_many_irreps>8
// Representation 216:
    zero_mom_loops_now[8][x]/=(transverse_volume*216);
      
#if how_many_irreps>9
// Representation 280:
    zero_mom_loops_now[9][x]/=(transverse_volume*280);
      
#if how_many_irreps>10
// Representation 330:
    zero_mom_loops_now[10][x]/=(transverse_volume*330);
      
#if how_many_irreps>11
// Representation 336:
    zero_mom_loops_now[11][x]/=(transverse_volume*336);

#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif
      
  }

  trace/=(nx*transverse_volume);
  average=trace;

  for (steps_along_x=0;steps_along_x<nx;steps_along_x++) {
    for (x=0;x<nx;x++) {
      zero_mom_corr_now[0][steps_along_x]+=zero_mom_loops_now[0][x]*conj(zero_mom_loops_now[0][(x+steps_along_x)%nx]);
    }
    zero_mom_corr_now[0][steps_along_x]/=nx;
  }

  for (x=0;x<nx;x++) {
    temporary_loop+=zero_mom_loops_now[0][x];
  }
  temporary_loop/=nx;

#ifdef __wanna_nonzero_mom_correlators__
  for (dist=0;dist<distances;dist++) {
    nonzero_mom_Polyakov_correlator=dc(0.,0.);
    for (x=0;x<nx;x++) 
    for (y=0;y<ny;y++)
#if dim==4
    for (z=0;z<nz;z++)
#endif
    {    
      spatial_site=y+ny*x;
      spatial_coordinate=((y+delta_y[dist])%ny)+ny*((x+delta_x[dist])%nx);
#if dim==4
      spatial_site*=nz;
      spatial_site+=z;
      spatial_coordinate*=nz;
      spatial_coordinate+=((z+delta_z[dist])%nz);
#endif
      nonzero_mom_Polyakov_correlator+=loop_through_x[spatial_site]*conj(loop_through_x[spatial_coordinate]);
    
    }
    nonzero_mom_Polyakov_correlator*=Pol_correlator_normalization_factor;
    nonzero_mom_corr_now[dist]=nonzero_mom_Polyakov_correlator;
  }
  
  delete [] loop_through_x;

#endif

  return average;

}
