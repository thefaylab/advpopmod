// A random effects version of the original SIMPLE.TPL example

DATA_SECTION
  init_int nstream
  init_int nobs
  init_matrix Y(1,nobs,1,nstream)
  !! cout << Y << endl;

PARAMETER_SECTION
  init_number beta
  init_bounded_number sigma_Y(0.000001,50)  // standard deviation of Y
  init_bounded_number sigma_b(0.000001,50)
  random_effects_vector b(1,nstream)
  vector pred_Y(1,nobs)

  objective_function_value f

PRELIMINARY_CALCS_SECTION

PROCEDURE_SECTION
  int i,j;
  
  f = 0;
  
  for (i=1;i<=nstream;i++)
   {
    f -= -log(sigma_b) -0.5*square(b(i)/sigma_b);
    for (j=1;j<=nobs;j++)
     f -= -log(sigma_Y) -0.5*square((Y(j,i)-(beta+b(i)))/sigma_Y);
   }  

