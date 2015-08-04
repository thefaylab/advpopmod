DATA_SECTION
  init_int m;
  init_vector T(1,m);
  init_int T_max;
  init_int L;

  init_matrix B(1,m,1,T_max+1-L);
  init_matrix R(1,m,1,T_max+1-L);

  init_vector phi0(1,m);

PARAMETER_SECTION
  init_bounded_number mu(-10,10,2);
  init_bounded_number log_tau(-7,4,2);
  init_bounded_vector B0(1,m,.1,1000);
  init_bounded_vector log_sigR(1,m,-7,4);
  random_effects_vector eps(1,m);

  objective_function_value obj_fun;
  sdreport_vector h(1,m);
  sdreport_vector R0(1,m);
  sdreport_number tau; 

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(11000);
  arrmblsize = 50000000; // in bytes
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(23888889);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(3000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(1000);
  
PRELIMINARY_CALCS_SECTION
  cout << setprecision(2);
 
GLOBALS_SECTION
  #include <df1b2fun.h>
//  #include <fvar.hpp>
//  #include <dflocmin.cpp>

PROCEDURE_SECTION
  int k;  
  dvariable beta;

  obj_fun = 0.0;
  R0 = elem_div(B0, phi0);
  tau = mfexp(log_tau);

  for (k=1; k<=m; k++) {
    stock_sep(k, eps(k), B0(k), log_sigR(k), mu, log_tau);
    beta = mu + tau * eps(k);
    h(k) = (mfexp(beta) + 0.2) / (1.0 + mfexp(beta));
  }


SEPARABLE_FUNCTION void stock_sep(int k_sep, const dvariable& eps_sep, const dvariable& B0_sep, const dvariable& log_sigR_sep, const dvariable& mu_sep, const dvariable& log_tau_sep)
  int i;
  dvariable tau = mfexp(log_tau_sep);
  dvariable R0 = B0_sep / phi0(k_sep);
  dvariable sigR = mfexp(log_sigR_sep);
  dvariable BH;

  obj_fun += .5*square(eps_sep);

  dvariable beta = mu_sep + tau * eps_sep;
  dvariable h = (mfexp(beta) + 0.2) / (1.0 + mfexp(beta));

  for (i=1; i<=T(k_sep)+1-L; i++) {
    BH = 4*R0*h*B(k_sep,i) / (B0_sep*(1-h) + (5*h-1)*B(k_sep,i));
    obj_fun += square(log(R(k_sep,i)) - log(BH) + square(sigR)/2) / (2*square(sigR)) + log(sigR);
  }

REPORT_SECTION
  report << "# mu\n" << mu << endl;
  report << "# log_tau\n" << log_tau << endl;
  report << "# tau\n" << tau << endl;
  report << "# h\n" << h << endl;
  report << "# log_sigR\n" << log_sigR << endl;
  report << "# m\n" << m << endl;
  report << "# T\n" << T << endl;
  report << "# B\n" << B << endl;
  report << "# R\n" << R << endl;
  report << "# phi0\n" << phi0 << endl;
  report << "# B0\n" << B0 << endl;
  report << "# R0\n" << R0 << endl;
  report << "# h\n" << h << endl;
  report << "# eps\n" << eps << endl;
