# sundir
This program simulates SU(N) Yang-Mills (2<=N<=9) in dimension d=3 or d=4.

In this folder sundir (./sundir) you will find:
    - parameters.h: 'dim' and 'col'
                    output_stream_index ---> output_stream_index=0 directs plaquette output to screen; otherwise to plaquette file







-----------------------------------------------------------------------
the main function for SCALAR usage works as follows:
- reading of argv (main argument): argv[1] = input_file_string  (basement name of the run)
                                   argv[2] = thread_number      (basement code of the run)

- reading of input file and filling of global input variables

- producing of interface mask (output direction is controlled through output_stream_index)

- call simulation(mode)
- call deallocate_arrays()


simulation.cc =========================================================
field variable are stored in an array of dc (complex double).
        - ufield is the pointer to the array of dc (> dc *ufield;)
        - the memory is allocated in the main ( 
                                                >     ufielddimension = nlinks*Ncolsquare
                                                >     ufield = new dc[ufielddimension];
                                              )



initialize.cc ---------------------------------------------------------
// The following convention for the site labelling
// is adopted:
//
//   site =          z +
//                y*nz +
//             x*ny*nz +
//          t*nx*ny*nz
//  
// Note that, in particular, on a 10^4 lattice, 
// the site of index txyz has coordinates (t,x,y,z).

link variables are initialized:
    0 ---> hot start (random dc numbers (1-2r) and then projected to SU(N))
    1 ---> cold start (identity matrix)



thermalize.cc ---------------------------------------------------------
Cycling over thermalization_time
  for (tempus=0; tempus<thermalization_time; tempus++) 
  { 
    update(locked_mode);                                           ==> update.cc (how_many_hb heat-bath updates, how_many_or overrelaxation update PER HB UPDATE ) E.G. hb=2, or=3 => (hb or or or) (hb or or or)        
    total_tempus++;                                                                                              
    if (total_tempus%reunitarization_period==0) {reunitarize();}   ==> project over (N)                                   
    output_stream << plaquette()                                   ==> stream plaquette                         
                                                                                                                                
    if (time_to_stop()) {                                               
      terminate_run();                                               
    }                                                                                              
  }








# sundir
