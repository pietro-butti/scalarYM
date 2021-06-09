void slab_average(int index, std::ofstream& multilevel_correlators_file) {

  int i, j, k, l, m, p, q, dir, update_ind;
  int slab, local_space, dist, local_t;
  int start_site, end_site;
  int t, x, y, site;
  int second_x, second_y, second_site;
  dc aux;
  
#if dim==3
  const int sites_in_one_slab=slab_size*nx*ny;
#endif

#if dim==4
  const int sites_in_one_slab=slab_size*nx*ny*nz;
  int z, second_z;
#endif

  for (slab=0;slab<how_many_slabs;slab++) {
    
    start_site=sites_in_one_slab*slab;
    end_site=sites_in_one_slab*(slab+1);

    for (i=0; i<atb_size; i++) {
      blocks_in_this_slab[i]=dc(0.,0.);
    }
    
    for (update_ind=0;update_ind<how_many_slab_updates;update_ind++) {
// locked_gauge_transform(true);
/**/
      for (unsigned short int hb_counter=0;hb_counter<how_many_hb;hb_counter++) {  
        for (site=start_site; site<end_site; site++)
        for (dir=0; dir<dim; dir++) {
          heat_bath_for_one_link(true, site, dir);
        }
        for (unsigned short int or_counter=0;or_counter<how_many_or;or_counter++)  
        for (site=start_site; site<end_site; site++)
        for (dir=0; dir<dim; dir++) {
          overrelaxation_for_one_link(true, site, dir);
        }
      }
/**/
      t=slab_size*slab;
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
        for (i=0;i<Ncolsquare;i++) {
          tempor[i]=ufield[site*Ncolsquare+i];
// tempor[i]=((i%Ncol_plus_one)==0)?dc(1.,0.): dc(0.,0.);
        }
        for (local_t=1;local_t<slab_size; local_t++) {
          for (i=0;i<Ncolsquare;i++) {
	    u2[i]=ufield[(site+spatial_volume*local_t)*Ncolsquare+i];
// u2[i]=((i%Ncol_plus_one)==0)?dc(1.,0.): dc(0.,0.);
          }
          mult_C_equals_AB(u1, tempor, u2);
          for (i=0;i<Ncolsquare;i++) {
	    tempor[i]=u1[i];
          }
        }
    
        for (dist=0; dist<distances; dist++) {

          second_x=(x+delta_x[dist])%nx;
          second_y=(y+delta_y[dist])%ny;
          second_site=second_y+ny*(second_x+nx*t);

#if dim==4
          second_z=(z+delta_z[dist])%nz;
          second_site*=nz;
          second_site+=second_z;
#endif

          for (i=0;i<Ncolsquare;i++) {
            tempor2[i]=ufield[second_site*Ncolsquare+i];
// tempor2[i]=((i%Ncol_plus_one)==0)?dc(1.,0.): dc(0.,0.);
          }
          for (local_t=1;local_t<slab_size; local_t++) {
            for (i=0;i<Ncolsquare;i++) {
              u4[i]=ufield[(second_site+spatial_volume*local_t)*Ncolsquare+i];
// u4[i]=((i%Ncol_plus_one)==0)?dc(1.,0.): dc(0.,0.);
            }
            mult_C_equals_AB(u3, tempor2, u4);
            for (i=0;i<Ncolsquare;i++) {
              tempor2[i]=u3[i];
            }
          }
        
          local_space=y+ny*x;

#if dim==4
          local_space*=nz;
          local_space+=z;
#endif
        
          for (i=0;i<Ncolsquare;i++)
          for (j=0;j<Ncolsquare;j++) {
            blocks_in_this_slab[
              ((i*Ncolsquare+j)*spatial_volume
                 +local_space)*distances
                 +dist]+=conj(u1[i])*u3[j];
          }
          
        } /* End of cycle over distances */

      } /* End of cycle over sites within a slab */

    } /* End of cycle over updates in this slab */
    
    if (slab==0) {
      for (i=0; i<atb_size; i++) {
        averaged_timelike_blocks[i]=blocks_in_this_slab[i]*inverse_of_how_many_slab_updates;
      }
    }
    else {
      for (m=0;m<spatial_volume_times_distances;m++) {
        for (i=0;i<Ncolfourth;i++) {
          largetemp1[i]=
            averaged_timelike_blocks[i*spatial_volume_times_distances+m];
          largetemp2[i]=
            blocks_in_this_slab[i*spatial_volume_times_distances+m]*inverse_of_how_many_slab_updates;
          largetemp3[i]=dc(0.,0.);
        }
        for (i=0;i<Ncol;i++)
        for (j=0;j<Ncol;j++)
        for (k=0;k<Ncol;k++)
        for (l=0;l<Ncol;l++)
        for (p=0;p<Ncol;p++)
        for (q=0;q<Ncol;q++) {
          largetemp3[(i*Ncol+j)*Ncolsquare+k*Ncol+l]+=
            largetemp1[(i*Ncol+p)*Ncolsquare+k*Ncol+q]*
            largetemp2[(p*Ncol+j)*Ncolsquare+q*Ncol+l];
	}
        for (i=0;i<Ncolfourth;i++) {
          averaged_timelike_blocks[i*spatial_volume_times_distances+m]=largetemp3[i];
	}
      }
    }
    
  } /* End of cycle over slabs */

  measure_correlators_from_averaged_timelike_blocks(index,  multilevel_correlators_file);
  
}
