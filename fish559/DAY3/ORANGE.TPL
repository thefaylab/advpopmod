// Orange data used in Ch. 8.2, Pinheiro & Bates (2000)

DATA_SECTION

  init_int n			// Number of data points
  init_vector y(1,n)		// Response vector
  init_vector t(1,n)		// Primary covariate
  init_int M			// Number of groups		
  init_vector ngroup(1,M)	// Group indicator
  init_int m			// Number of parameters in nonlinear regression model		

PARAMETER_SECTION

  init_bounded_vector beta(1,3,-50,2000,1)       // Fixed effects parameters
  init_bounded_number log_sigma(-6.0,5.0,1)	// log(residual variance)
  init_bounded_number log_sigma_u(-10,5,3)		// 0.5*log(variance component)
  random_effects_vector u(1,M,2)			// Unscaled random effects
  objective_function_value g
  number sigma_u;
  number sigma;


PROCEDURE_SECTION
  int i,j,ii;
  dvariable Pred;

  g = 0.0;
  ii = 0;
  sigma = mfexp(log_sigma);
  sigma_u = mfexp(log_sigma_u);
  for(i=1;i<=M;i++)
  {
   g -= -log(sigma_u) -0.5*square(u(i)/sigma_u);
   for (j=1;j<=n/M;j++)
    {
     ii++;
     Pred = (beta(1)+u(i))/(1+exp(-1*(t(ii)-beta(2))/beta(3)));
     g -= -log(sigma) -0.5*square((y(ii)-Pred)/sigma);
    }
  }



