void read_last_conf() {

  char inputfilename[maximum_filename_length];

  sprintf(inputfilename,"datadir/saved_conf_%s.dat", rundetails);

  ifstream inputfile;
  inputfile.open(inputfilename, ifstream::binary);
  if (inputfile.is_open()) {
    inputfile.read((char*)ufield, sizeof(dc)*ufielddimension);
    inputfile.close();
  }
  else { 
    cout << "# Unable to open file" << endl;
    exit(0);
  }

}
