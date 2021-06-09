double analyze_trace_modulus(char *input_filename, int how_many_header_lines, int which_bin, int irrep_label, int row_where_measurement_is, int rows_in_each_measurement, double *autocorrelation_time, double histogram_minvalue, double histogram_maxvalue, int number_of_histograms, double *histograms, bool *out_of_range_flag) {

  
  
  int row_counter=0;
  int irrep_counter, misurate, m;
  int meas_counter, j;
  int dummy_int=0;
  char comment_line[maximum_filename_length];
  double dummy_re, dummy_im, dummy=0.;
  double dati[how_many_measurements];
  double act, autoc0, q, value=0.;

  double step=(histogram_maxvalue-histogram_minvalue)/((double) number_of_histograms);
  const double one_meas_histogram_increase=1./(step*how_many_measurements);
  int histogram_label;
    
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
        histogram_label=floor((dummy-histogram_minvalue)/step);
        if ( (histogram_label<0) || (histogram_label>=number_of_histograms) ) {
	  *out_of_range_flag=true;
	}
	else {
	  histograms[histogram_label]+=one_meas_histogram_increase;
	}
      }
      value+=dummy;
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

  input_file.close();

  if (which_bin==quantibins) {

    for (meas_counter=0;meas_counter<how_many_measurements;meas_counter++)
    dati[meas_counter]-=value;
    act=0.5;
    autoc0=0.0;
    for (meas_counter=0;meas_counter<how_many_measurements;meas_counter++) {
      autoc0+=dati[meas_counter]*dati[meas_counter];
    }
    autoc0/=how_many_measurements;
    m=1;

    for (meas_counter=1;meas_counter<how_many_measurements;meas_counter++) {
      q=0.0;
      for (j=0;j<how_many_measurements-meas_counter;j++) {
        q+=dati[j]*dati[j+meas_counter];
      }
      q/=(how_many_measurements-meas_counter);
      if (meas_counter<4*act+1) { 
        act+=q/autoc0;
        m++;
      }
    }

    *autocorrelation_time=act;
//     cout << "Autocorrelation time: " << act << endl;

  }

  return value;

}
