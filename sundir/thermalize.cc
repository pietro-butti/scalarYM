void thermalize() {

  int tempus;
  int total_tempus=0;
  bool locked_mode=
#ifdef __wanna_Jarzynski_SF__
    true
#else
    false
#endif
  ;

  output_stream << "# average"
#ifndef __wanna_Jarzynski_SF__
                << "        timelike       spacelike"
#endif
                << endl;
  for (tempus=0;tempus<thermalization_time;tempus++) { 
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
// cout << "# Test: " << plaquette() << " "; gauge_transform(); cout << plaquette() << " "; check_unitarity(); cout << plaquette() << endl; 
    if (time_to_stop()) {
      terminate_run();
    }
  }

}
