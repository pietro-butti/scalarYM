double analyze_loop_susceptibility(char *input_filename, int how_many_header_lines, int which_bin, int irrep_label, int row_where_measurement_is, int rows_in_each_measurement) {

  int row_counter=0;
  int irrep_counter, misurate;
  int meas_counter;
  int dummy_int=0;
  char comment_line[maximum_filename_length];
  double dummy_re, dummy_im, dummy=0.;
  double dati[how_many_measurements];
  double squared_dati[how_many_measurements];
  double value=0.;
  double squared_value=value=0.;

  if (how_many_measurements%quantibins) {
    cout << "Error: the number of bins should be a divisor of the total number of measurements" << endl;
    exit(0);
  }

  ifstream input_file(input_filename);

  while (row_counter<how_many_header_lines) {
    input_file.getline(comment_line,maximum_filename_length);
    row_counter++;
  }

  misurate=0;
  for (meas_counter=0;meas_counter<how_many_measurements;meas_counter++) {
    if ((which_bin!=quantibins) && (meas_counter/binsize==which_bin)) {
      for (row_counter=0;row_counter<rows_in_each_measurement;row_counter++) {
        input_file.getline(comment_line,maximum_filename_length);
      }
    }
    else {
      for (row_counter=0;row_counter<row_where_measurement_is-1;row_counter++) {
        input_file.getline(comment_line,maximum_filename_length);
      }
      input_file >> dummy_int;
      for (irrep_counter=0; irrep_counter<irrep_label+1; irrep_counter++) {
        input_file >> dummy_re;
        input_file >> dummy_im;
      }
      dummy=sqrt(dummy_re*dummy_re + dummy_im*dummy_im);
      
// cout << dummy << endl;
      if (which_bin==quantibins) {
        dati[misurate]=dummy;
        squared_dati[misurate]=dummy*dummy;
      }
      value+=dummy;
      squared_value+=dummy*dummy;
      misurate++;
      for (irrep_counter=irrep_label+1; irrep_counter<how_many_irreps; irrep_counter++) {
        input_file >> dummy_re;
        input_file >> dummy_im;
      }
      for (row_counter=row_where_measurement_is-1;row_counter<rows_in_each_measurement;row_counter++) {
        input_file.getline(comment_line,maximum_filename_length);
      }
    }
  }

  value/=misurate;
  squared_value/=misurate;
  value=squared_value-value*value;

  input_file.close();

  return value;

}
