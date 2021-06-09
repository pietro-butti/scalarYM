double analyze_one_observable(char *input_filename, int how_many_header_lines, int which_bin, int how_many_columns, int column_where_measurement_is, int row_where_measurement_is, int rows_in_each_measurement, double *autocorrelation_time) {

  int row_counter=0;
  int column_counter, misurate, m;
  int meas_counter, j;
  char comment_line[maximum_filename_length];
  double dummy;
  double dati[how_many_measurements];
  double act, autoc0, q, value=0.;


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
      for (column_counter=0; column_counter<column_where_measurement_is; column_counter++) {
        input_file >> dummy;
      }
// cout << dummy << endl;
      if (which_bin==quantibins) {
        dati[misurate]=dummy;
      }
      value+=dummy;
      misurate++;
      for (column_counter=column_where_measurement_is; column_counter<how_many_columns; column_counter++) {
        input_file >> dummy;
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
