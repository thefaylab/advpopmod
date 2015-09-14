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
  vector y;
  !! streams = (ivector)column(thedata,1);
  !! y = column(thedata,2);

PARAMETER_SECTION
  init_number ln_beta;
  init_number ln_sigma;

  vector ypred(1,ndata);
  sdreport_number sigma;
  sdreport_number beta;

  objective_function_value obj_fun;

PROCEDURE_SECTION
  //transformations
  sigma = mfexp(ln_sigma);
  beta = mfexp(ln_beta);

  //predictions
  ypred = mfexp(ln_beta);

  //normal likelihood for the observations
  obj_fun += ndata*ln_sigma + norm2(y-ypred)/square(sigma);

