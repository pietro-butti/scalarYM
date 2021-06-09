void save_last_conf() {

  char outputfilename[maximum_filename_length];
  ofstream outputfile;

  sprintf(outputfilename,"datadir/saved_conf_%s.dat", rundetails);
  outputfile.open(outputfilename, ofstream::binary);
  if (outputfile.is_open()) {
    outputfile.write((char*)ufield, sizeof(dc)*ufielddimension);
  }
  else { 
    cout << "# Unable to open file" << endl;
    exit(0);
  }
  outputfile.close();
  
}
