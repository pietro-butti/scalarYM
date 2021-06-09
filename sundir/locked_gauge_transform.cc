void locked_gauge_transform(bool locked_mode) {

  int site, mu, i, neighbor;

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

  for (site=0;site<nsites;site++)
#ifdef __wanna_multilevel__
  if ((locked_mode==false) || (
#ifdef __wanna_improvement__
      ((((site/spatial_volume)+1)%slab_size)!=0) &&
      ((((site/spatial_volume)+slab_size-1)%slab_size)!=0) &&
#endif    
      (((site/spatial_volume)%slab_size)!=0) )
     )
#endif
  {

    for (i=0;i<Ncolsquare;i++) {
      v[i]=dc(xx.randNorm(0.,1.),xx.randNorm(0.,1.));
    }
    norm_sun(v);

    for (mu=0;mu<dim;mu++) {

      for (i=0;i<Ncolsquare;i++) {
        oldu[i] = ufield[(mu*nsites+site)*Ncolsquare+i];
      }

      mult_C_equals_AB(u1,v,oldu);

      for (i=0;i<Ncolsquare;i++) {
        ufield[(mu*nsites+site)*Ncolsquare+i] = u1[i];
      }


      neighbor=neighbor_minus[mu*nsites+site];

      for (i=0;i<Ncolsquare;i++) {
        oldu[i] = ufield[(mu*nsites+neighbor)*Ncolsquare+i];
      }

      mult_C_equals_ABdagger(u1,oldu,v);

      for (i=0;i<Ncolsquare;i++) {
        ufield[(mu*nsites+neighbor)*Ncolsquare+i] = u1[i];
      }
    }
  }

}
