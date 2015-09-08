//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Gavin Fay
//
//  Lab 1: Linear regression example
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
  //declare counter variable
  int i;
  //change subject to 1 if only 1 a parameter
  !!if (num_apar==1) for(int i=1;i<=ndata;i++) subject(i) = 1;
   
PARAMETER_SECTION
  //init_number dummy;
  init_vector ln_a(1,num_apar);
  init_number b;
  
  vector log_wtpred(1,ndata);
  
  sdreport_matrix newpred(1,20,1,num_apar);

  objective_function_value obj_fun;

PRELIMINARY_CALCS_SECTION
  //perform some calculations on the data prior to estimation
  cout << subject << endl;
  //cout << len << endl;
  //cout << weight << endl;
  loglen = log(len);
  logwt = log(weight);
  //cout << loglen << endl;
  //cout << logwt << endl;
 

PROCEDURE_SECTION
  for (i=1;i<=ndata;i++)
   log_wtpred(i) = ln_a(subject(i)) + b*loglen;

  obj_fun = norm2(logwt-log_wtpred);
  //obj_fun = square(dummy);

