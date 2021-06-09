void update(bool locked_mode) {

#ifdef __wanna_trace_deformation__
  for (unsigned short int metro_counter=0; metro_counter<how_many_metro; metro_counter++) {  
    metropolis(locked_mode);
  }  
#else
  for (unsigned short int hb_counter=0;hb_counter<how_many_hb;hb_counter++) {  
    heat_bath(locked_mode);
// gauge_transform();
    for (unsigned short int or_counter=0;or_counter<how_many_or;or_counter++) { 
      overrelaxation(locked_mode);
// gauge_transform();
    }
  }
#endif


}
