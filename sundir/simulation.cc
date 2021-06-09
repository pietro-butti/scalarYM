void simulation(int mode) {

  initialize();
  if (mode==0) {
    thermalize();
  }
  else {
    read_last_conf();
  }
	
  update(0);
  /*
  measurements();
  */
  terminate_run();

  // QUI BISOGNEREBBE AGGIUNGERE LA PARTE PER SALVARE LE CONFIGURAZIONI

}
