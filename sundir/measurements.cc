void measurements() {

/* Produces the files with: average Wilson loops, */
/* average Polyakov loops and their _disconnected_ correlators */

  unsigned long int tempus, total_tempus, measurement;
  
#ifdef __wanna_Wilson_loops__
  char *Wilson_filename = new char[maximum_filename_length];
#ifdef __wanna_smearing__
  char *smeared_Wilson_filename = new char[maximum_filename_length];
#endif
#endif


#ifdef __wanna_zero_mom_loops_and_correlators__
  int steps_along_x, irrep, x;
  char average_Polyakov_filename[maximum_filename_length];
  char zero_mom_Polyakov_correlators_filename[maximum_filename_length];
  dc average_Polyakov_loop;
  dc temp=dc(0.,0.);
  dc **zero_mom_loops_now = new dc* [how_many_irreps];
  dc **zero_mom_corr_now = new dc* [how_many_irreps];
  for (irrep=0;irrep<how_many_irreps;irrep++) {
    zero_mom_loops_now[irrep]= new dc[nx];
    zero_mom_corr_now[irrep]= new dc[nx];
    for (x=0;x<nx;x++) {
      zero_mom_loops_now[irrep][x]=dc(0.,0.);
      zero_mom_corr_now[irrep][x]=dc(0.,0.);
    }  
  }
#ifdef __wanna_nonzero_mom_correlators__
  char nonzero_mom_Polyakov_correlators_filename[maximum_filename_length];
  int dist;
  dc *nonzero_mom_corr_now = new dc [distances];
  for (dist=0;dist<distances;dist++) {
    nonzero_mom_corr_now[dist]=dc(0.,0.);
  }  
#endif
#endif

#ifdef __wanna_zero_mom_loops_and_correlators__
  sprintf(average_Polyakov_filename,"datadir/average_Polyakov_%s.dat", rundetails);
  sprintf(zero_mom_Polyakov_correlators_filename,"datadir/zero_mom_Polyakov_correlators_%s.dat", rundetails);
#ifdef __wanna_nonzero_mom_correlators__
  sprintf(nonzero_mom_Polyakov_correlators_filename,"datadir/nonzero_mom_Polyakov_correlators_%s.dat", rundetails);
#endif
#endif

#ifdef __wanna_Wilson_loops__
  sprintf(Wilson_filename,"datadir/Wilson_%s.dat", rundetails);
#ifdef __wanna_smearing__
  sprintf(smeared_Wilson_filename,"datadir/smeared_Wilson_%s.dat", rundetails);
#endif
#endif

#ifdef __wanna_multilevel__
  sprintf(multilevel_correlators_filename,"datadir/multilevel_correlators_%s.dat", rundetails);
#endif

#ifdef __wanna_Jarzynski_SF__
  char *Jarzynski_SF_filename = new char[maximum_filename_length];
  sprintf(Jarzynski_SF_filename,"datadir/Jarzynski_SF_%s.dat", rundetails);
#endif


#ifdef __wanna_zero_mom_loops_and_correlators__
  average_Polyakov_file.open(average_Polyakov_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  average_Polyakov_file << setprecision(precision_digits);
  average_Polyakov_file << scientific;

  zero_mom_Polyakov_correlators_file.open(zero_mom_Polyakov_correlators_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  zero_mom_Polyakov_correlators_file << setprecision(precision_digits);
  zero_mom_Polyakov_correlators_file << scientific;
  
#ifdef __wanna_nonzero_mom_correlators__
  nonzero_mom_Polyakov_correlators_file.open(nonzero_mom_Polyakov_correlators_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  nonzero_mom_Polyakov_correlators_file << setprecision(precision_digits);
  nonzero_mom_Polyakov_correlators_file << scientific;
#endif
#endif


#ifdef __wanna_Wilson_loops__
  Wilson_loop_file.open(Wilson_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  Wilson_loop_file << setprecision(precision_digits);
  Wilson_loop_file << scientific;
#ifdef __wanna_smearing__
  smeared_Wilson_loop_file.open(smeared_Wilson_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  smeared_Wilson_loop_file << setprecision(precision_digits);
  smeared_Wilson_loop_file << scientific;
#endif
#endif

#ifdef __wanna_multilevel__
  multilevel_correlators_file.open(multilevel_correlators_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  multilevel_correlators_file << setprecision(precision_digits);
  multilevel_correlators_file << scientific;
#endif

#ifdef __wanna_Jarzynski_SF__
  Jarzynski_SF_file.open(Jarzynski_SF_filename
#ifdef __wanna_append__
  , ios_base::app
#endif
  );
  Jarzynski_SF_file << setprecision(precision_digits);
  Jarzynski_SF_file << scientific;
#endif

  bool locked_mode=
#ifdef __wanna_Jarzynski_SF__
    true
#else
    false
#endif
  ;

  total_tempus=thermalization_time;
#ifdef __wanna_zero_mom_loops_and_correlators__
  dc temp_loop[how_many_irreps];
#endif
  for (measurement=0;measurement<(unsigned long int) how_many_measurements;measurement++) { 
    for (tempus=0;tempus<(unsigned long int) updates_between_measurements;tempus++) { 
      
#ifdef __wanna_multilevel__
      if ((measurement+1)%how_many_slab_updates!=0) { 
	locked_mode=true;
      } 
      else {
	locked_mode=false;
      }
#endif

      update(locked_mode);
      total_tempus++;
      if (total_tempus%reunitarization_period==0) {
        reunitarize();
      }
      output_stream << plaquette()
#ifndef __wanna_Jarzynski_SF__
                    << "   " << time_plaq() << " " << space_plaq()
#endif
                    << endl;
    }
    
#ifdef __wanna_multilevel__
    slab_average(measurement, multilevel_correlators_file);
#endif

#ifdef __wanna_Jarzynski_SF__
    measure_Jarzynski_SF(measurement, Jarzynski_SF_file);
#endif
    
#ifdef __wanna_Wilson_loops__
    measure_Wilson_loops(measurement, Wilson_loop_file,0);
#ifdef __wanna_smearing__
    smear();
    measure_Wilson_loops(measurement, smeared_Wilson_loop_file,1);
#endif
#endif

#ifdef __wanna_zero_mom_loops_and_correlators__
    for (irrep=0; irrep<how_many_irreps; irrep++) 
    for (x=0; x<nx; x++) {
      zero_mom_loops_now[irrep][x]=dc(0.,0.);
    }
    
    average_Polyakov_loop=measure_Polyakov_loops(measurement, zero_mom_loops_now, zero_mom_corr_now
#ifdef __wanna_nonzero_mom_correlators__
    , nonzero_mom_corr_now
#endif
    );
    
    average_Polyakov_file << measurement+measurements_already_done << "   ";
    for (irrep=0; irrep<how_many_irreps; irrep++) {
      temp_loop[irrep]=dc(0.,0.);
      for (x=0; x<nx; x++) {
        temp_loop[irrep]+=zero_mom_loops_now[irrep][x];
      }
      temp_loop[irrep]/=nx;
      average_Polyakov_file << "   " << real(temp_loop[irrep]) << " " << imag(temp_loop[irrep]);
    } 
    average_Polyakov_file << endl;
        
    temp+=average_Polyakov_loop;
    
    for (steps_along_x=0;steps_along_x<nx;steps_along_x++) {
      zero_mom_Polyakov_correlators_file << measurement+measurements_already_done << " " << steps_along_x << " " << real(zero_mom_corr_now[0][steps_along_x])/*-norm(temp)*/ << " " << imag(zero_mom_corr_now[0][steps_along_x])/*-norm(temp)*/ << endl; 
    }
    
#ifdef __wanna_nonzero_mom_correlators__
    nonzero_mom_Polyakov_correlators_file << measurement+measurements_already_done << "   ";
    for (dist=0;dist<distances;dist++) {
      nonzero_mom_Polyakov_correlators_file << real(nonzero_mom_corr_now[dist]) << " " << imag(nonzero_mom_corr_now[dist]);
      if (dist<distances-1) {
	nonzero_mom_Polyakov_correlators_file << "  ";
      }
      else {
	nonzero_mom_Polyakov_correlators_file << endl;
      }
    }
#endif
    
#endif

    measurements_done_until_now=measurement+measurements_already_done+1;
    if (time_to_stop()) {
      terminate_run();
    }

  }
  
#ifdef __wanna_zero_mom_loops_and_correlators__
  temp/=how_many_measurements;
#endif

#ifdef __wanna_Wilson_loops__
  Wilson_loop_file.close();
#ifdef __wanna_smearing__
  smeared_Wilson_loop_file.close();
#endif
#endif

#ifdef __wanna_multilevel__
  multilevel_correlators_file.close();
#endif

#ifdef __wanna_zero_mom_loops_and_correlators__  
  for (irrep=0;irrep<how_many_irreps;irrep++) {
    delete [] zero_mom_loops_now[irrep];
    delete [] zero_mom_corr_now[irrep];
  }
  delete [] zero_mom_loops_now;
  delete [] zero_mom_corr_now;
  average_Polyakov_file.close();
  zero_mom_Polyakov_correlators_file.close();
#ifdef __wanna_nonzero_mom_correlators__
  delete [] nonzero_mom_corr_now;
  nonzero_mom_Polyakov_correlators_file.close();
#endif
#endif

#ifdef __wanna_Jarzynski_SF__
  delete [] Jarzynski_SF_filename;
#endif
}
