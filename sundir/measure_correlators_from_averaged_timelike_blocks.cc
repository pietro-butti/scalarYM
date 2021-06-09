void measure_correlators_from_averaged_timelike_blocks(int index, std::ofstream& multilevel_correlators_file) {

  dc loop_correlator, contribution;
  int dist, local_space, i, j;

  multilevel_correlators_file << index << "   ";
    
  for (dist=0; dist<distances; dist++) {
    loop_correlator=dc(0.,0.);    
    for (local_space=0; local_space<spatial_volume; local_space++) 
    for (i=0; i<Ncolsquare; i+=Ncol_plus_one)
    for (j=0; j<Ncolsquare; j+=Ncol_plus_one) {
      loop_correlator+=averaged_timelike_blocks[
            ((i*Ncolsquare+j)*spatial_volume
               +local_space)*distances
               +dist];
    }
    loop_correlator*=Pol_correlator_normalization_factor;
    multilevel_correlators_file << real(loop_correlator) << " "<<imag(loop_correlator);
    if (dist<distances-1) {
      multilevel_correlators_file << "  ";
    }
    else {
      multilevel_correlators_file << endl;
    }
  }

}
