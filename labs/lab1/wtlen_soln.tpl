//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Gavin Fay
//
//  Lab 1: Weight length relationships
//////////////////////////////////////////////
DATA_SECTION
  init_int num_apar;   //number of a parameters
  init_int ndata;      //number of data points   
  init_ivector subject(1,ndata);
  init_vector len(1,ndata);
  init_vector weight(1,ndata);
  //declare variables that will be functions of the data
  vector loglen(1,ndata);
  vector logwt(1,ndata);
  //declare counter variables
  int i;
  int ilen;
  //change subject to 1 if only 1 a parameter
  !!if (num_apar==1) for(int i=1;i<=ndata;i++) subject(i) = 1;
   
PARAMETER_SECTION
  //init_number dummy;             //dummy parameter to use during debugging
  init_vector ln_a(1,num_apar);    //natural logs of the a parameter
  init_number b;                   //exponent of the w-l relationship
  
  vector log_wtpred(1,ndata);      //predicted values for the data

  matrix newpred(1,20,1,num_apar); //matrix to fill in predictions for lengths 1-20

  objective_function_value obj_fun;   //objective function

PRELIMINARY_CALCS_SECTION
  //perform some calculations on the data prior to estimation
  //cout << subject << endl;
  //cout << len << endl;
  //cout << weight << endl;
  loglen = log(len);
  logwt = log(weight);
  //cout << loglen << endl;
  //cout << logwt << endl;
 

PROCEDURE_SECTION
  //predict the ln_weights
  //for (i=1;i<=ndata;i++)
  // log_wtpred(i) = ln_a(subject(i)) + b*loglen(i);
  
  log_wtpred = ln_a(subject) + b*loglen;

  //calculate the Residual sums of squares
  //obj_fun = square(dummy);   // evaluate the objective function when using the dummy
  obj_fun = norm2(logwt-log_wtpred);   //residual sums of squares of the logged weights
  
  
REPORT_SECTION
  // predict the weights for lengths 1 through 20
  for (ilen=1;ilen<=20;ilen++)
   newpred(ilen) = mfexp(ln_a+b*log(ilen));

  //write to report file
  report << "weights for lengths 1 through 20" << endl;
  for (ilen=1;ilen<=20;ilen++)
   report << ilen << " " << newpred(ilen) << endl;
