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
