void measure_Jarzynski_SF(int configuration, std::ofstream& output_file) {

  double previous_action, current_action, action_difference;
// double exponential_average;
  int i, J_therm;
  int step_index;
  const bool locked_mode=true;
  
  action_difference=0.;

  for (i=0; i<ufielddimension; i++) {
    stored_ufield[i]=ufield[i];
  }

  for (step_index=0; step_index<nsteps; step_index++) {
    previous_action=compute_action_from_first_boundary();
    previous_action+=compute_action_from_second_boundary();
    update_first_boundary(step_index);
    update_second_boundary(step_index);
    current_action=compute_action_from_first_boundary();
    current_action+=compute_action_from_second_boundary();
    action_difference+=current_action-previous_action;
    if (step_index<nsteps_minus_one) {
      for (J_therm=0;J_therm<Jarzynski_thermalization;J_therm++) {
        update(locked_mode);
      }
      reunitarize();
    }
  }

//   exponential_average=exp(-action_difference);
//   output_file << configuration+measurements_already_done << " " << exponential_average << endl;

  output_file << configuration+measurements_already_done << " " << action_difference << " " << exp(-action_difference) 
  << endl;

  for (i=0; i<ufielddimension; i++) {
    ufield[i]=stored_ufield[i];
//  Note that this also restores the second boundary to its initial value
  }

}
