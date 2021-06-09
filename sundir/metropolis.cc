void metropolis(bool locked_mode) {

  int origin, site, dir, i, next_site;

  for (dir=0;dir<dim;dir++) {
#ifdef __wanna_trace_deformation__
    if (is_compactified[dir]==true) {
      switch (dir) {
        case 0:
	  compactified_direction_size=nt;
	  break;
        case 1:
	  compactified_direction_size=nx;
	  break;
        case 2:
	  compactified_direction_size=ny;
	  break;
#if dim>3
        case 3:
	  compactified_direction_size=nz;
	  break;
#endif
	default:
	  break;
      }
      double_trace_prefactor=compactified_direction_size*compactified_direction_size;
#if dim>3
      double_trace_prefactor*=compactified_direction_size;
#endif
    }        
#endif
    double_trace_prefactor=1./double_trace_prefactor;
    max_t=nt;
    max_x=nx;
    max_y=ny;
#if dim>3
    max_z=nz;
#endif
    switch (dir) {
      case 0:
        max_t=1;
        break;
      case 1:
        max_x=1;
        break;
      case 2:
        max_y=1;
        break;
#if dim>3
      case 3:
        max_z=1;
        break;
#endif
      default:
        break;
    }
  
    for (int t=0;t<max_t;t++)
    for (int x=0;x<max_x;x++)
    for (int y=0;y<max_y;y++)
#if dim>3
    for (int z=0;z<max_z;z++)
#endif
    {
      origin=y+ny*(x+nx*t);
#if dim>3
      origin*=nz;
      origin+=z;
#endif
      site=origin;
      open_Polyakov_loop(open_loop, dir, site);
      do {
#ifdef __wanna_multilevel__
        if ((locked_mode==false) ||
          (locked_link[dir*nsites+site]==false) 
           )
#endif
        {
	
          for (i=0;i<Ncolsquare;i++) {
	    oldu[i]=ufield[(dir*nsites+site)*Ncolsquare+i];
          }
        
          oldaction_minus_newaction=compute_action(site, dir, oldu
#ifdef __wanna_trace_deformation__
            , open_loop
#endif
          );
          generate_SUN_matrix_near_identity(SUN_matrix);
	  metro_tried++;
	  mult_C_equals_AB(newblock,oldu,SUN_matrix);
          oldaction_minus_newaction-=compute_action(site, dir, newblock
#ifdef __wanna_trace_deformation__
            , open_loop
#endif
          );
	
	  if ( (oldaction_minus_newaction>0.) || 
	       (xx.rand()<exp(oldaction_minus_newaction)) ) {
            for (i=0;i<Ncolsquare;i++) {
	      ufield[(dir*nsites+site)*Ncolsquare+i]=newblock[i];
	    }
	    metro_accepted++;
	  }
	
	  next_site=neighbor_plus[dir*nsites+site];
          for (i=0;i<Ncolsquare;i++) {
	    u1[i]=ufield[(dir*nsites+next_site)*Ncolsquare+i];
	    u2[i]=open_loop[i];
	    u3[i]=ufield[(dir*nsites+site)*Ncolsquare+i];
          }
          mult_C_equals_AdaggerB(u4,u1,u2);
          mult_C_equals_AB(open_loop,u4,u3);
        
        }     
#ifdef __wanna_multilevel__
        else {
	  next_site=neighbor_plus[dir*nsites+site];
	}
#endif      
        site=next_site;
      
      } while (site!=origin);
    }
  }
    
  if (metro_tried>=max_tried) {
    metro_acceptance=((double) metro_accepted)/metro_tried;
// cout << "# Metropolis acceptance = " << metro_acceptance;
// cout << "  twiceamplitude = " << twiceamplitude << endl;
    if (metro_acceptance<0.25) {
      twiceamplitude*=amplitude_resizing_factor;
    }
    if (metro_acceptance>0.5) {
      twiceamplitude/=amplitude_resizing_factor;
    }
    metro_accepted=0;
    metro_tried=0;
  }

}
