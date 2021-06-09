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
