//#############################
//## yelloweye maturity at age
//## MAR 580 Fall 2015
//## Lab 10
//## Gavin Fay
//## Last modified 11/17/15
//#############################

DATA_SECTION
  init_int Nobs
  init_matrix TheData(1,Nobs,1,3)
  int iobs

  !!CLASS ofstream MCMCreport("MCMCreport.out",ios::trunc);

PARAMETER_SECTION
  init_number a
  init_bounded_number delta(1,10)

  sdreport_number age50
  sdreport_number b

  vector p(1,Nobs)

  objective_function_value objfun

PROCEDURE_SECTION
  // transform the parameters
  b = (log(0.95)-log(0.05))/delta;
  age50 = -1.*a/b;

  // calculate the prior for age50
  // (we don't need one for delta because it is uniform 
  // and it's range constraint is placed above)
  objfun += log(0.3) + 0.5*square(age50-15)/square(0.3);
  
  //calculate the probabilities for the binomial
  p = a + b*column(TheData,1);
  p = elem_div(mfexp(p),(1.+mfexp(p)));
  
  // calculate the likelihood
  for (iobs=1;iobs<=Nobs;iobs++)
   objfun -= TheData(iobs,3)*log(p(iobs)) + (TheData(iobs,2)-TheData(iobs,3))*log(1.e-07+1.-p(iobs));

  //cout << objfun << endl;
  if (mceval_phase())
   {
    MCMCreport << objfun << " " << a << " " << b << " " << age50 << " " << delta << endl;
   }

REPORT_SECTION
  report << age50 << " " << delta << endl;
  report << p << endl;
