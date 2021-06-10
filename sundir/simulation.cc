void simulation() {

  initialize();

  dc *proposal = new dc[ufielddimension];
  for(int ii=0; ii<ufielddimension; ii++) proposal[ii] = ufield[ii];



  Metropolis_for_one_link(2,1,.05, proposal);





	




}




/*
void simulation(int mode) {

  initialize();
  if (mode==0) {
    thermalize();
  }
  else {
    read_last_conf();
  }

 
  measurements();
 
  terminate_run();

  // QUI BISOGNEREBBE AGGIUNGERE LA PARTE PER SALVARE LE CONFIGURAZIONI

}
*/
