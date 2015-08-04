DATA_SECTION
  init_int NData;
  init_vector x(1,NData);
  init_vector y(1,NData);
  vector yobs(1,NData);
  
 !!  int II;
 !!  for (II=1;II<=NData;II++)
 !!   yobs(II) = log(y(II)/(1-y(II)));

PARAMETER_SECTION  
  init_number logA50;
  init_number logA95;
  init_number logSigma;

  number A50;
  number A95;
  number Sigma;

  vector ypred(1,NData);
  sdreport_number A502
  objective_function_value obj_fun;
  
INITIALIZATION_SECTION
  logA50 3
  logA95 4
  logSigma -1

PROCEDURE_SECTION
  int II;
  
  A50 = mfexp(logA50);
  A502 = A50;
  A95 = mfexp(logA95);
  Sigma = mfexp(logSigma);

  for (II = 1;II<=NData;II++)
   {
    ypred(II) = 1.0/(1+mfexp(-1*log(19)*(x(II)-A50)/(A95-A50)));
    ypred(II) = log(ypred(II)/(1-ypred(II)));
   }
  
  obj_fun = NData*log(Sigma) + norm2(yobs-ypred)/(2*square(Sigma));
//  cout << obj_fun << endl;

  if (mceval_phase())
   cout << A50 << " " << A95 << " " << Sigma << " " << obj_fun << endl;
    
