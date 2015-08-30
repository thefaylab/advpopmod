//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Gavin Fay
//
//  Lab 1: Linear regression example
//////////////////////////////////////////////
DATA_SECTION
  init_int ndata;
  init_vector x(1,ndata);
  init_vector y(1,ndata);

PARAMETER_SECTION
  init_number a;
  init_number b;
  
  vector ypred(1,ndata);
  objective_function_value obj_fun;

PROCEDURE_SECTION
  ypred = a+b*x;
  obj_fun = norm2(y-ypred);

