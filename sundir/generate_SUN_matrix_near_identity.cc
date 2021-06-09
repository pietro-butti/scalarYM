inline void generate_SUN_matrix_near_identity(dc *v) {
  
  int i;
  
  void (*norm_sun)(dc *v);
  
#if Ncol==2 
  norm_sun=norm_su2;
#endif

#if Ncol==3 
  norm_sun=norm_su3;
#endif

#if Ncol==4
  norm_sun=norm_su4;
#endif

#if Ncol==5
  norm_sun=norm_su5;
#endif

#if Ncol==6
  norm_sun=norm_su6;
#endif

#if Ncol==7
  norm_sun=norm_su7;
#endif

#if Ncol==8
  norm_sun=norm_su8;
#endif

#if Ncol==9
  norm_sun=norm_su9;
#endif

#if Ncol==10
  norm_sun=norm_su10;
#endif

  for (i=0;i<Ncolsquare;i++) {
    aux_matrix[i]=dc(twiceamplitude*(0.5-xx.rand()),twiceamplitude*(0.5-xx.rand()));
  }
  for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
    aux_matrix[i]+=dc(1.,0.);
  }
  norm_sun(aux_matrix);
  for (i=0;i<Ncolsquare;i++) {
    v[i]=aux_matrix[i];
  }
  
}

