//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Gavin Fay
//
//  Lecture 3: Streams
//////////////////////////////////////////////

DATA_SECTION
  init_int ndata;
  init_int nstreams;
  init_matrix thedata(1,ndata,1,2);
  ivector streams(1,ndata) 
  vector y(1,ndata);
  !! streams = (ivector)column(thedata,1);
  !! y = column(thedata,2);

PARAMETER_SECTION
  init_number ln_beta;
  init_number ln_sigma;
  init_number ln_sigma_b(2);

  vector ypred(1,ndata);
  sdreport_number sigma;
  sdreport_number sigma_b;
  sdreport_number beta;

  random_effects_vector b(1,nstreams,2);

  objective_function_value obj_fun;

PROCEDURE_SECTION
  //distribution for the random effects
  // Normal (0,1) random variables
  obj_fun = log(1.) + 0.5*norm2(b)/square(1.);
  // we can also write this as obj_fun = 0.5*norm2(b);

  //transformations
  sigma = mfexp(ln_sigma);
  sigma_b = mfexp(ln_sigma_b);
  beta = mfexp(ln_beta);

  //predictions
  ypred = mfexp(ln_beta) + sigma_b*b(streams);

  //normal likelihood for the observations
  obj_fun += ndata*ln_sigma + 0.5*norm2(y-ypred)/square(sigma);

REPORT_SECTION
  report << "random effects" << endl;
  report << sigma_b*b << endl;
