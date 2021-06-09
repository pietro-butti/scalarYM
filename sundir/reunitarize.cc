void reunitarize() {

  int dir, site, i;
  void (*norm_sun)(dc *v);

#  if Ncol==2 
  norm_sun=norm_su2;
#endif

#  if Ncol==3 
  norm_sun=norm_su3;
#endif

#  if Ncol==4
  norm_sun=norm_su4;
#endif

#  if Ncol==5
  norm_sun=norm_su5;
#endif

#  if Ncol==6
  norm_sun=norm_su6;
#endif

#  if Ncol==7
  norm_sun=norm_su7;
#endif

#  if Ncol==8
  norm_sun=norm_su8;
#endif

#  if Ncol==9
  norm_sun=norm_su9;
#endif

#  if Ncol==10
  norm_sun=norm_su10;
#endif

  for (site=0;site<nsites;site++)
  for (dir=0;dir<dim;dir++) {
    for (i=0;i<Ncolsquare;i++) {
      v[i]=ufield[(dir*nsites+site)*Ncolsquare+i];
    };
    norm_sun(v);
    for (i=0;i<Ncolsquare;i++) {
      ufield[(dir*nsites+site)*Ncolsquare+i]=v[i];
    }
  }

}
