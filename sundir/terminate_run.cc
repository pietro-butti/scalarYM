void terminate_run() {

#ifdef __wanna_save_last_conf__
  save_last_conf();
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

#ifdef __wanna_Jarzynski_SF__
  Jarzynski_SF_file.close();
#endif


#ifdef __wanna_zero_mom_loops_and_correlators__
  if (average_Polyakov_file.is_open()) {
    average_Polyakov_file.close();
  }
  if (zero_mom_Polyakov_correlators_file.is_open()) {
    zero_mom_Polyakov_correlators_file.close();
  }
#ifdef __wanna_nonzero_mom_correlators__
  if (nonzero_mom_Polyakov_correlators_file.is_open()) {
    nonzero_mom_Polyakov_correlators_file.close();
  }
#endif
#endif

  
  time(&rawtime);
  ptm=gmtime(&rawtime);
  timeinfo=localtime(&rawtime);
//   printf("# Helsinki time:  %2d:%02d\n", (ptm->tm_hour+HEL_time)%24, ptm->tm_min);
  output_stream << "# Local time and date at run end: " << asctime(timeinfo);
//   cout << "# Current local time and date at run end: " << asctime(timeinfo) << endl;
    
  measfile.open(measfilename,fstream::out);
  measfile << " " << measurements_done_until_now << endl;
  measfile.close();
  
#if output_stream_index != 0
  plaquettefile.close();
#endif
  
}
