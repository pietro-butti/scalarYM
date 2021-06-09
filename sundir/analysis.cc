#include "headers.h"
#include "parameters.h"

using namespace std;

typedef complex<double> dc;

double analyze_one_observable(char *input_filename, int how_many_header_lines, int which_bin, int how_many_columns, int column_where_measurement_is, int row_where_measurement_is, int rows_in_each_measurement, double *autocorrelation_time);

char *input_filename=(char *) malloc(maximum_filename_length*sizeof(char));
char *output_filename=(char *) malloc(maximum_filename_length*sizeof(char));


int main() {

  int contabins=how_many_bins;
  
// cout << "contabins=" << contabins << endl; exit(0);
  
  double media;
  double valore;
  double errore=0.;
  double autocorrelation_time;
  double other_time;
  
  int how_many_header_lines=0;
  
  int how_many_columns=4;
  int column_where_measurement_is=3;
  int row_where_measurement_is;
  int rows_in_each_measurement=distances;
  
#if dim==3
    sprintf(input_filename,"datadir/multilevel_correlators_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf.dat", Ncol, nt, nx, ny, beta);
#endif

#if dim==4
    sprintf(input_filename,"datadir/multilevel_correlators_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf.dat", Ncol, nt, nx, ny, nz, beta);
#endif

/*     sprintf(input_filename,"datafile.dat"); */
  
  cout << "# Analysis of the file " << input_filename << endl << "#" << endl;
  
  cout << "# r/a  average_value errorbar  integrated_autocorrelation_time #" << endl;

  for (int dist=0; dist<distances; dist++) {

    row_where_measurement_is=dist+1;
    
    media=0.;
    errore=0.;
     
    media = analyze_one_observable(input_filename, 
      how_many_header_lines,
      contabins,
      how_many_columns,
      column_where_measurement_is,
      row_where_measurement_is,
      rows_in_each_measurement, 
      &autocorrelation_time) ;

    for (contabins=0;contabins<how_many_bins;contabins++) {

      valore=analyze_one_observable(input_filename, 
        how_many_header_lines,
        contabins,
        how_many_columns,
        column_where_measurement_is,
        row_where_measurement_is,
        rows_in_each_measurement, 
        &other_time) ;
            
      valore-=media;
      valore*=valore;
      errore+=valore;
    }
    errore*=(how_many_bins-1.)/how_many_bins;
    errore=sqrt(errore);
    printf("%2d    %12.10lf %12.10lf   %12.10lf", dist, media, errore, autocorrelation_time);
    printf("\n");
  }
}

///////////////////////////////////////////////////
///////////////////////////////////////////////////
///////////////////////////////////////////////////

double analyze_one_observable(char *input_filename, int how_many_header_lines, int which_bin, int how_many_columns, int column_where_measurement_is, int row_where_measurement_is, int rows_in_each_measurement, double *autocorrelation_time) {

/*

This function can be used to calculate the average value from a set of measurements of an observable in a data file, or to extract the binned averages to evaluate the jackknife error. 

In input: 

-  input_filename: the name of the input file, whose structure is:

     comment lines
     measurements (blocks of lines)

   each measurement can consists of several lines, and each line can contain several data

- how_many_header_lines: number of header lines (to be skipped)

- which_bin: the number of the bin that we are considering; the program calculates the average over all measurements, except those in the which_bin'th bin; which_bin ranges from 0 to how_many_bins

- how_many_columns: the total number of columns in the row where the measurement of interest is. The other rows may have a different number of columns

- column_where_measurement_is: the column of the data file where the data of interest are written; it is understood that column_where_measurement_is=2 means that the data is in the second (not in the third) column of the file.

- row_where_measurement_is: since each measurement is a block which may include several lines, row_where_measurement_is indicates in which row of each block the datum of interest is.

- rows_in_each_measurement: number of rows in each measurement.

- *autocorrelation_time: pointer to the double variable in which the integrated autocorrelation time should be stored. This is done only when the function is called in averaging (not in binning) mode, i.e. only when which_bin=how_many_bins.


An example of use is:



  int contabins=how_many_bins;
  double media;
  double valore;
  double errore=0.;
  double autocorrelation_time

  media = analyze_one_observable("datafile.dat", 1, contabins, 8, 5, 1, 2, &autocorrelation_time) ;
  cout << endl << "   Media: "<< media << endl;

  cout << endl <<endl <<endl <<endl;
  for (contabins=0;contabins<how_many_bins;contabins++) {
    valore=analyze_one_observable("datafile.dat", 1, contabins, 8, 5, 1, 2, &autocorrelation_time);
    cout << endl << "   media di questo bin: "<< valore << endl;
    valore-=media;
    valore*=valore;
    errore+=valore;
  }
  errore*=(how_many_bins-1.)/how_many_bins;
  errore=sqrt(errore);
  cout << endl << endl << "Eccolo: " << errore << endl;



where the datafile.dat file could be of the form:

# commento
1 2 0. 0 1. 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 2 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 3 0 4 5
0 0 0 0 0 0 0 0
1 2 0. 0 4 0 4 5
0 0 0. 0 0 0 0 0
1 2 0 0. 5 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 6 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 7 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 8 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 9 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 10 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 11 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 12 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 13 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 14 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 15 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 16 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 17 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 18 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 19 0 4 5
0 0 0 0 0 0 0 0
1 2 0 0 20 0 4 5
0 0 0 0 0 0 0 0

*/


  int row_counter=0;
  int column_counter, misurate, m;
  int meas_counter, j;
  char comment_line[maximum_filename_length];
  double dummy;
  double dati[how_many_measurements];
  double act, autoc0, q, value=0.;


  if (how_many_measurements%how_many_bins) {
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
// cout << misurate << endl;
    if ((which_bin!=how_many_bins) && (meas_counter/binsize==which_bin)) {
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
      if (which_bin==how_many_bins) {
        dati[misurate]=dummy;
// if ((which_bin==2) && (misurate==2)) {printf("%12.10e\n",dummy); exit(0);}
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

  if (which_bin==how_many_bins) {

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

  }

  return value;

}


