#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdio>
#include <ctime>
#include <complex>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#define pi 3.1415926535898
#define twopi 6.28318530717958648
#define maximum_filename_length 2000
#define precision_digits 12
#define tol 1.5e-31

double beta;
int Ncol;
int dim;
int nt;
int nspacelike;
int nx;
int ny;
int nz;
int how_many_measurements;
const int how_many_bins=10;
const int quantibins=how_many_bins;
int binsize;



#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_math.h>
gsl_multifit_fdfsolver *s;

int min_L;

#define cut_large_r 100//7
#define MAX_ITER 500
#define FIT_TOL 1.e-4

using namespace std;

typedef complex<double> dc;

#define __smeared__
char *input_filename=(char *) malloc(maximum_filename_length*sizeof(char));
char *output_filename_1=(char *) malloc(maximum_filename_length*sizeof(char));
char *output_filename_2=(char *) malloc(maximum_filename_length*sizeof(char));
char *output_filename_3=(char *) malloc(maximum_filename_length*sizeof(char));

int maxloopsize;

double value_here, previous_value, error_here, previous_error;

size_t ndata_for_constant_fit;
double *r1;

size_t ndata_for_cornell_fit;
const size_t number_of_cornell_fit_parameters = 3;
double *r2;

double analyze_one_observable(char *input_filename, int how_many_header_lines, int which_bin, int how_many_columns, int column_where_measurement_is, int row_where_measurement_is, int rows_in_each_measurement, double *autocorrelation_time);



struct data {
  size_t num;
  double * y;
  double * errorbar;
};

double *asymptotic_Vr;
double *error_asymptotic_Vr;
  
int cornell_f(const gsl_vector *x, void *data, gsl_vector *f);
int cornell_df(const gsl_vector *x, void *data, gsl_matrix *J);
int cornell_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J); 

int constant_f(const gsl_vector *x, void *data, gsl_vector *f);
int constant_df(const gsl_vector *x, void *data, gsl_matrix *J);
int constant_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J);


int main(int argc, char **argv) {

  if (argc!=8) {
    cout << "Usage: ./potential_analysis.x <beta> <Ncol> <dim> <nt> <nspacelike> <min_L> <how_many_measurements>" << endl; 
    exit(0);
  };

  beta=(double) atof(argv[1]);
  Ncol=atoi(argv[2]);
  dim=atoi(argv[3]);
  nt=atoi(argv[4]);
  nspacelike=atoi(argv[5]);
  min_L=atoi(argv[6]);
  how_many_measurements=atoi(argv[7]);
  
  if ((dim!=3)&&(dim!=4)) {
    cout << "The number of spacetime dimensions must be either 3 or 4" << endl;
  }
  nx=nspacelike;
  ny=nspacelike;
  nz=nspacelike;
  maxloopsize=((nspacelike+1)/2);
  ndata_for_constant_fit = maxloopsize-min_L+1;
  r1 = new double[ndata_for_constant_fit];
  ndata_for_cornell_fit = maxloopsize;
  r2 = new double[ndata_for_cornell_fit];
  asymptotic_Vr = new double[maxloopsize+1];
  error_asymptotic_Vr = new double[maxloopsize+1];
  if ((how_many_measurements%how_many_bins)!=0) {
    cout << "The number of measurements must be a multiple of the number of bins" << endl;
  }
  binsize=how_many_measurements/how_many_bins;

  char prefix[200];
#ifdef __smeared__
  sprintf(prefix,"smeared_");
#else
  sprintf(prefix,"");
#endif

  int contabins=how_many_bins;
  
// cout << "contabins=" << contabins << endl; exit(0);
  
  double media;
  double valore;
  double errore=0.;
  double autocorrelation_time, other_time;
  
  int how_many_header_lines=0;
  
  int how_many_columns=5;
  int column_where_measurement_is=4;
  int row_where_measurement_is=1;
  int rows_in_each_measurement=(nspacelike/2)*(nt/2);
  
  double average_Vr[maxloopsize+1][maxloopsize+1];
  double error_Vr[maxloopsize+1][maxloopsize+1];
  double jack_Vr[how_many_bins][maxloopsize+1][maxloopsize+1];
  
  for (int L1=0; L1<=maxloopsize; L1++)
  for (int L2=0; L2<=maxloopsize; L2++) {
    average_Vr[L1][L2]=0.;
    error_Vr[L1][L2]=0.;
    for (contabins=0;contabins<how_many_bins;contabins++) {
      jack_Vr[contabins][L1][L2]=0.;
    }
  }
  
//   double histogram_maxvalue=1.;
//   double histogram_minvalue=0.;
//   int number_of_histograms=100; 
//   double *histograms=(double *) malloc(number_of_histograms*sizeof(double));
//   bool out_of_range_flag=false;

  if (dim==3) {
    sprintf(input_filename,"datadir/%sWilson_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, beta);
    sprintf(output_filename_1,"%sWilson_mean_values_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, beta);
    sprintf(output_filename_2,"%spotential_at_finite_L_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, beta);
    sprintf(output_filename_3,"%spotential_Ncol_%d_nt_%d_nx_%d_ny_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, beta);
  }
  else {
    sprintf(input_filename,"datadir/%sWilson_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, nz, beta);
    sprintf(output_filename_1,"%sWilson_mean_values_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, nz, beta);
    sprintf(output_filename_2,"%spotential_at_finite_L_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, nz, beta);
    sprintf(output_filename_3,"%spotential_Ncol_%d_nt_%d_nx_%d_ny_%d_nz_%d_beta_%12.10lf.dat", prefix, Ncol, nt, nx, ny, nz, beta);
  }
  
/*     sprintf(input_filename,"datafile.dat"); */
  

  ofstream output_1, output_2, output_3;
  output_1.open(output_filename_1);
  output_2.open(output_filename_2);
  output_3.open(output_filename_3);
  output_1 << setprecision (precision_digits);
  output_2 << setprecision (precision_digits);
  output_3 << setprecision (precision_digits);
  cout << "# Analysis of the file " << input_filename << endl << "#" << endl;
//   cout << "# r  L  average_value error    V(r) at finite L  error   V(r) error #" << endl;


// Note: In the raw data file the first column is L (the unsmeared time direction)
// and the second is r (smeared directions), but, for convenience (because we 
// have to calculate differences at fixed r) here the convention is reversed:
// see the definition of row_where_measurement_is below, which assumes that,
// in the input file L1 varies fastest (=r).
  for (int L1=1; L1<=maxloopsize; L1++) {
    previous_value=0.;
    previous_error=0.;
    for (int L2=1; L2<=maxloopsize; L2++) {

      row_where_measurement_is=1+(L1-1)+maxloopsize*(L2-1);
    
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

        valore=analyze_one_observable(input_filename, how_many_header_lines, contabins, how_many_columns, column_where_measurement_is, row_where_measurement_is, rows_in_each_measurement, &other_time) ;

        valore-=media;
        valore*=valore;
        errore+=valore;
      }
      errore*=(how_many_bins-1.)/how_many_bins;
      errore=sqrt(errore);
      printf("%2d %2d    ", L1, L2);
//       printf("%12.10lf %12.10lf", media, errore);
      output_1 << L1 << " "<< L2 << " "<< media << " " << errore << endl;
      value_here=log(media);
      error_here=errore/media;    
      if (L2>1) {
	average_Vr[L1][L2]=previous_value-value_here;
//         printf("    %12.10lf %12.10lf", previous_value-value_here, sqrt(previous_error*previous_error+error_here*error_here));
        printf(" %12.10lf", average_Vr[L1][L2]); 
        output_2 << L1 << " "<< L2 << " "<< average_Vr[L1][L2] << " " << sqrt(previous_error*previous_error+error_here*error_here) << endl; 
        if (L2==maxloopsize) {
//           printf("    %12.10lf %12.10lf", previous_value-value_here, sqrt(previous_error*previous_error+error_here*error_here));
	  if (L1>=cut_large_r) {
	    output_3 << "# ";
	  }
          output_3 << L1 << " "<< average_Vr[L1][L2] << " " << sqrt(previous_error*previous_error+error_here*error_here) << endl;  
        }   
      }
      previous_value=value_here;
      previous_error=error_here;
      printf("\n");
 
    }
  }
  output_1.close();
  output_2.close();
  output_3.close();
  
    
  for (contabins=0;contabins<how_many_bins;contabins++) 
  for (int L1=1; L1<=maxloopsize; L1++) {
    previous_value=0.;
    previous_error=0.;
    for (int L2=1; L2<=maxloopsize; L2++) {

      row_where_measurement_is=1+(L1-1)+maxloopsize*(L2-1);
    
      media=0.;
      errore=0.;
     
      media = analyze_one_observable(input_filename, how_many_header_lines, contabins, how_many_columns, column_where_measurement_is, row_where_measurement_is, rows_in_each_measurement, &other_time) ;
      
//       printf("%d  %2d %2d    ", contabins, L1, L2);
      value_here=log(media);
      if (L2>1) {
	jack_Vr[contabins][L1][L2]=previous_value-value_here;
//         printf(" %12.10lf", previous_value-value_here); 
      }
      previous_value=value_here;
      previous_error=error_here;
//       printf("\n");
 
    }
  }
  
  
  for (int L1=0; L1<=maxloopsize; L1++)
  for (int L2=0; L2<=maxloopsize; L2++) {
    for (contabins=0;contabins<how_many_bins;contabins++) {
      jack_Vr[contabins][L1][L2]-=average_Vr[L1][L2];
      jack_Vr[contabins][L1][L2]*=jack_Vr[contabins][L1][L2];
      error_Vr[L1][L2]+=jack_Vr[contabins][L1][L2];
    }
    error_Vr[L1][L2]*=(how_many_bins-1.)/how_many_bins;
    error_Vr[L1][L2]=sqrt(error_Vr[L1][L2]);
    if ( (L1>0) && (L2>1) ) {
      printf("%d %d  %12.10lf %12.10lf\n",L1,L2,average_Vr[L1][L2],error_Vr[L1][L2]);
    }
  }
  
//////////////////////////////////////////////////  
  
  
  const gsl_multifit_fdfsolver_type *T;
//   gsl_multifit_fdfsolver *s;
  int status;
  unsigned int iter = 0;
  const size_t number_of_constant_fit_parameters = 1;

  gsl_matrix *covar = gsl_matrix_alloc (number_of_constant_fit_parameters, number_of_constant_fit_parameters);
  double y[ndata_for_constant_fit], errorbar[ndata_for_constant_fit];
  struct data d = { ndata_for_constant_fit, y, errorbar};
  gsl_multifit_function_fdf f;
  double x_init[number_of_constant_fit_parameters] = { 0.1 };
//   double x_init[number_of_parameters] = { 0.1, 0.0, -0.25 };
  gsl_vector_view x = gsl_vector_view_array (x_init, number_of_constant_fit_parameters);
  const gsl_rng_type * type;
  gsl_rng * r;

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);
  
  f.f = &constant_f;
  f.df = &constant_df;
  f.fdf = &constant_fdf;
  f.n = ndata_for_constant_fit;
  f.p = number_of_constant_fit_parameters;
  f.params = &d;

  for (int L1=0; L1<=maxloopsize; L1++) {
    asymptotic_Vr[L1]=0.;
    error_asymptotic_Vr[L1]=0.;
  }

  for (int L1=1; L1<=maxloopsize; L1++) {

    printf("%12.10lf    ", (double) L1);
  /* This is the data to be fitted */
    for (int i=0; i< (int) ndata_for_constant_fit; i++) {
      r1[i]=min_L+i;
      y[i]=average_Vr[L1][min_L+i];
      errorbar[i]=error_Vr[L1][min_L+i];
//       printf ("data: %g  %g %g\n", r1[i], y[i], errorbar[i]);
    }

    T = gsl_multifit_fdfsolver_lmsder;
    s = gsl_multifit_fdfsolver_alloc (T, ndata_for_constant_fit, number_of_constant_fit_parameters);
    gsl_multifit_fdfsolver_set (s, &f, &x.vector);

//     print_state (iter, s, number_of_constant_fit_parameters);
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate (s);

//       printf ("status = %s\n", gsl_strerror (status));
//       print_state (iter, s, number_of_constant_fit_parameters);

      if (status)
        break;

      status = gsl_multifit_test_delta (s->dx, s->x, FIT_TOL, FIT_TOL);
    } while (status == GSL_CONTINUE && iter < MAX_ITER);

    
    gsl_multifit_covar (s->J, 0.0, covar);

    { 
      double chi = gsl_blas_dnrm2(s->f);
      double dof = ndata_for_constant_fit - number_of_constant_fit_parameters;
      double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

      for (int i=0; i< (int) number_of_constant_fit_parameters; i++) {
        printf ("%12.10lf %12.10lf   ", gsl_vector_get(s->x, i), c*sqrt(gsl_matrix_get(covar,i,i)));
      }
      printf("  chisq/dof = %g\n",  pow(chi, 2.0) / dof);
      
      asymptotic_Vr[L1]=gsl_vector_get(s->x,0);
      error_asymptotic_Vr[L1]=c*sqrt(gsl_matrix_get(covar,0,0));
      
//       printf ("sigma  = %.5f +/- %.5f\n", FIT(0), c*ERR(0));
//       printf ("V0     = %.5f +/- %.5f\n", FIT(1), c*ERR(1));
//       printf ("gamma  = %.5f +/- %.5f\n", FIT(2), c*ERR(2));

    }

//     printf ("status = %s\n", gsl_strerror (status));

  }
  gsl_multifit_fdfsolver_free (s);
  gsl_matrix_free (covar);
  gsl_rng_free (r);

  
  
////////////////////////////////////////////////////////////////

    
  const gsl_multifit_fdfsolver_type *T2;
  gsl_multifit_fdfsolver *s2;

  gsl_matrix *covar2 = gsl_matrix_alloc (number_of_cornell_fit_parameters, number_of_cornell_fit_parameters);
  double y2[ndata_for_cornell_fit], errorbar2[ndata_for_cornell_fit];
  struct data d2 = { ndata_for_cornell_fit, y2, errorbar2};
  gsl_multifit_function_fdf f2;
  double x_init2[number_of_cornell_fit_parameters] = { 0.1, 0.0, 0.0 };
  gsl_vector_view x2 = gsl_vector_view_array (x_init2, number_of_cornell_fit_parameters);

  gsl_rng_env_setup();

  type = gsl_rng_default;
  r = gsl_rng_alloc (type);

  f2.f = &cornell_f;
  f2.df = &cornell_df;
  f2.fdf = &cornell_fdf;
  f2.n = ndata_for_cornell_fit;
  f2.p = number_of_cornell_fit_parameters;
  f2.params = &d2;

  
//   printf("%12.10lf    ", beta);
  /* This is the data to be fitted */
  for (int i=0; i < (int) ndata_for_cornell_fit; i++) {
    r2[i]=i+1.;
    y2[i]=asymptotic_Vr[i+1];
    errorbar2[i]=error_asymptotic_Vr[i+1];
//     printf ("data: %g  %g %g\n", r2[i], y2[i], errorbar2[i]);
  }

  T2 = gsl_multifit_fdfsolver_lmsder;
  s2 = gsl_multifit_fdfsolver_alloc (T2, ndata_for_cornell_fit, number_of_cornell_fit_parameters);
  gsl_multifit_fdfsolver_set (s2, &f2, &x2.vector);

  iter=0;
//   print_state (iter, s2, number_of_cornell_fit_parameters);
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate (s2);
//     printf ("status = %s\n", gsl_strerror (status));
//     print_state (iter, s2, number_of_cornell_fit_parameters);

    if (status)
      break;

    status = gsl_multifit_test_delta (s2->dx, s2->x, FIT_TOL, FIT_TOL);
  } while (status == GSL_CONTINUE && iter < MAX_ITER);

    
  gsl_multifit_covar (s2->J, 0.0, covar2);

  { 
    double chi = gsl_blas_dnrm2(s2->f);
    double dof = ndata_for_cornell_fit - number_of_cornell_fit_parameters;
    double c = GSL_MAX_DBL(1, chi / sqrt(dof)); 

    for (int i=0; i < (int) number_of_cornell_fit_parameters; i++) {
      printf ("%s\t = %12.10lf %12.10lf   \n",(i==0)?"sigma":(i==1)?"V0":"gamma", gsl_vector_get(s2->x, i), c*sqrt(gsl_matrix_get(covar2,i,i)));
    }
    printf("  chisq/dof = %g\n",  pow(chi, 2.0) / dof);
    
  }

  gsl_multifit_fdfsolver_free (s2);
  gsl_matrix_free (covar2);
  gsl_rng_free (r);

  delete r1;
  delete r2;
  delete asymptotic_Vr;
  delete error_asymptotic_Vr;
  return 0;
  
}

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



/* cornellfit.cc -- model functions for Cornell potential */
    
int cornell_f (const gsl_vector * x, void *data, gsl_vector * f) {

  size_t num = ((struct data *)data)->num;
  double *y = ((struct data *)data)->y;
  double *errorbar = ((struct data *) data)->errorbar;
     
  double stringtension = gsl_vector_get (x, 0);
  double V0 = gsl_vector_get (x, 1);
  double gamma = gsl_vector_get (x, 2);
     
  size_t i;
  
  for (i = 0; i < num; i++) {
/* Model Yi = stringtension * r2[i] + V0 + gamma/r2[i] */
    double t1 = r2[i];
    double Yi = stringtension * t1 + V0 + gamma/t1;
    gsl_vector_set (f, i, (Yi - y[i])/errorbar[i]);
  }
     
  return GSL_SUCCESS;

}
     
int cornell_df (const gsl_vector * x, void *data, gsl_matrix * J) {

  size_t num = ((struct data *)data)->num;
  double *errorbar = ((struct data *) data)->errorbar;
     
//        double stringtension = gsl_vector_get (x, 0);
//        double V0 = gsl_vector_get (x, 1);
//        double gamma = gsl_vector_get (x, 2);
     
  size_t i;
     
  for (i = 0; i < num; i++) {
/* Jacobian matrix J(i,j) = dfi / dxj, */
/* where fi = (Yi - yi)/errorbar[i],      */
/*       Yi = stringtension * r2[i] + V0 + gamma/r2[i]  */
/* and the xj are the parameters (stringtension,V0,gamma) */
    double t1 = r2[i];
    double s = errorbar[i];
    gsl_matrix_set (J, i, 0, t1/s ); 
    gsl_matrix_set (J, i, 1, 1./s );
    gsl_matrix_set (J, i, 2, 1./(t1*s) );
  }
  return GSL_SUCCESS;
}
     
int cornell_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
  cornell_f (x, data, f);
  cornell_df (x, data, J);
     
  return GSL_SUCCESS;
}



/* constantfit.cc -- model functions for fit to a constant */
     
int constant_f (const gsl_vector * x, void *data, gsl_vector * f) {

  size_t num = ((struct data *)data)->num;
  double *y = ((struct data *)data)->y;
  double *errorbar = ((struct data *) data)->errorbar;
     
  double constantvalue = gsl_vector_get (x, 0);
     
  size_t i;
     
  for (i = 0; i < num; i++) {
/* Model Yi = constantvalue */
    double Yi = constantvalue;
    gsl_vector_set (f, i, (Yi - y[i])/errorbar[i]);
  }
  return GSL_SUCCESS;
}
     
int constant_df (const gsl_vector * x, void *data, gsl_matrix * J) {

  size_t num = ((struct data *)data)->num;
  double *errorbar = ((struct data *) data)->errorbar;

  size_t i;
     
  for (i = 0; i < num; i++) {
/* Jacobian matrix J(i,j) = dfi / dxj, */
/* where fi = (Yi - yi)/errorbar[i],      */
/*       Yi = constantvalue  */
/* and the xj are the parameters (constantvalue) */
    double s = errorbar[i];
    gsl_matrix_set (J, i, 0, 1./s ); 
  }
  return GSL_SUCCESS;
}
    
int constant_fdf (const gsl_vector * x, void *data, gsl_vector * f, gsl_matrix * J) {
  constant_f (x, data, f);
  constant_df (x, data, J);
  return GSL_SUCCESS;
}

