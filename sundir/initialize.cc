void initialize() {

  int t, x, y, site, dir, next_t, next_x, next_y, next_site, i;
#if dim==4 
  int z, next_z;
#endif

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

#ifdef __wanna_Wilson_loops__
  if (nx<=ny) {
    max_Wilson_loop_size=nx/2;
  }
  else {
    max_Wilson_loop_size=ny/2;
  }
#if dim==4  
  if ((nz/2)<max_Wilson_loop_size) {
    max_Wilson_loop_size=nz/2;
  }
#endif
  if ((first_mu==0) && ((nt/2)<max_Wilson_loop_size)) {
    max_Wilson_loop_size=nt/2;
  }
#endif

#ifdef __wanna_Jarzynski_SF__

//   if ((nx!=ny)
// #if dim==4 
//       || (nx!=nz)
// #endif
//   ) {
//     printf("In the present version of the code, the Jarzynski's theorem routines require the spatial sizes of the lattice to be equal\n");
//     exit(0);
//   }
  
#ifdef __wanna_multilevel__
  printf("In the present version of the code, it is not possible to use the multilevel in combination with Jarzynski's theorem\n");
  exit(0);
#endif
  
#ifdef __wanna_improvement__
  printf("In the present version of the code, the Jarzynski's theorem routines do not support the improved gauge action\n");
  exit(0);
#endif

#endif


// Conventions for the site labelling:
//
//   site =          z +
//                y*nz +
//             x*ny*nz +
//          t*nx*ny*nz
//  
// Note that, in particular, on a 10^4 lattice, the site of index
// txyz has coordinates (t,x,y,z).
  
  for (t=0;t<nt;t++)
  for (x=0;x<nx;x++)
  for (y=0;y<ny;y++)
#if dim==4  
  for (z=0;z<nz;z++)
#endif
  {
    site=y+ny*(x+nx*t);
#if dim==4
    site*=nz;
    site+=z;
#endif
    for (dir=0;dir<dim;dir++) {
#ifdef __wanna_Jarzynski_SF__
      if ((t==(nt-1)) || ((t==0) && (dir!=0))) {
// Links to be held fixed: spatial links at t=0
// and all links at t=nt-1:
        locked_link[dir*nsites+site]=true;
      }
      else {
        locked_link[dir*nsites+site]=false;
      };
#endif
#ifdef __wanna_multilevel__
      if ((t%slab_size==0) && (dir!=0)) {
        locked_link[dir*nsites+site]=true;
      }
      else {
        locked_link[dir*nsites+site]=false;
      };
#ifdef __wanna_improvement__
      if ((t+1)%slab_size==0) {
        locked_link[dir*nsites+site]=true;
      }
#endif
#endif
      next_t=t;
      next_x=x;
      next_y=y;
#if dim==4  
      next_z=z;
#endif
      switch (dir) {
        case 0: {
          next_t=(t+1)%nt;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z;
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_t=(t+nt-1)%nt;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4 
          next_site*=nz;
          next_site+=next_z; 
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        case 1: {
          next_x=(x+1)%nx;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4 
          next_site*=nz;
          next_site+=next_z; 
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_x=(x+nx-1)%nx;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z; 
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
        case 2: {
          next_y=(y+1)%ny;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z; 
#endif
          neighbor_plus[dir*nsites+site]=next_site;
          next_y=(y+ny-1)%ny;
          next_site=next_y+ny*(next_x+nx*next_t);
#if dim==4  
          next_site*=nz;
          next_site+=next_z; 
#endif
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
#if dim==4  
        case 3: {
          next_z=(z+1)%nz;
          next_site=next_z+nz*(next_y+ny*(next_x+nx*next_t));
          neighbor_plus[dir*nsites+site]=next_site;
          next_z=(z+nz-1)%nz;
          next_site=next_z+nz*(next_y+ny*(next_x+nx*next_t));
          neighbor_minus[dir*nsites+site]=next_site;
          break;
        }
#endif
        default: {
          break;
        }
      }

      if (type_of_start==0) {
        for (i=0;i<Ncolsquare;i++) {
          v[i]=dc(1.-2.*xx.rand(),1.-2.*xx.rand());
        }
        norm_sun(v);
        for (i=0;i<Ncolsquare;i++) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=v[i];
        }
      }

      if (type_of_start==1) {



        for (i=0;i<Ncolsquare;i++) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
        }
        for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
          ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
        }

        
      }
/*
//       if (type_of_start==3) {
// 	if ( ( (site!=0) && (site!=ny) ) || (dir!=0) ) {
//           for (i=0;i<Ncolsquare;i++) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
//           }
//           for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
//           }
//         }
//         else {
//           for (i=0;i<Ncolsquare;i++) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(0.,0.);
//           }
//           for (i=0;i<Ncolsquare;i+=Ncol_plus_one) {
//             ufield[(dir*nsites+site)*Ncolsquare+i]=dc(1.,0.);
//           }
//           ufield[(dir*nsites+site)*Ncolsquare+0]=dc(0.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+1]=dc(1.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+Ncol_plus_one]=dc(0.,0.);
//           ufield[(dir*nsites+site)*Ncolsquare+Ncol]=dc(-1.,0.);
// 	}
//       }
*/
    }
  }
  
#ifdef __wanna_Jarzynski_SF__

  int ind;
  for (i=0; i<Ncolsquare; i++) {
    C0[i]=dc(0.,0.);
    C1[i]=dc(0.,0.);
  }
  for (i=0; i<Ncol; i++) {
    phi0[i]*=division_by_L;
    newphi0[i]*=division_by_L;
    phi1[i]*=division_by_L;
    newphi1[i]*=division_by_L;
    C0[i*(Ncol_plus_one)]=exp(dc(0.,phi0[i]));
    C1[i*(Ncol_plus_one)]=exp(dc(0.,phi1[i]));
  }

  for (x=0;x<nx;x++)
  for (y=0;y<ny;y++)
#if dim==4  
  for (z=0;z<nz;z++)
#endif
  {
    site=(y+ny*(x+nx*(nt-1)));
#if dim==4
    site*=nz;
    site+=z;
#endif
    for (dir=0;dir<dim;dir++) {
        ind=(dir*nsites+site)*Ncolsquare;
// Set all links in t=nt-1 slice to C1:
        for (i=0;i<Ncolsquare;i++) {
          ufield[ind+i]=C1[i];
        }
    }
    site=y+ny*x;
#if dim==4
    site*=nz;
    site+=z;
#endif
    for (dir=1;dir<dim;dir++) {
        ind=(dir*nsites+site)*Ncolsquare;
// Set spatial links in t=0 slice to C0:
        for (i=0;i<Ncolsquare;i++) {
          ufield[ind+i]=C0[i];
        }
    }
  }

#endif

  if (type_of_start==2) {
    read_last_conf();
//     char inputfilename[maximum_filename_length];
//     FILE *inputfile;
//     size_t fread_result;
// #if dim==3
//     sprintf(inputfilename,"saved_conf_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf_thread_%d.dat", Ncol, nt, nx, ny, beta, thread_number);
// #endif
// #if dim==4
//     sprintf(inputfilename,"saved_conf_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf_thread_%d.dat", Ncol, nt, nx, ny, nz, beta, thread_number);
// #endif
//     inputfile=fopen(inputfilename,"r");
//     fread_result=fread( ufield, sizeof(ufield), 1, inputfile);
//     fclose(inputfile);
  }

}
