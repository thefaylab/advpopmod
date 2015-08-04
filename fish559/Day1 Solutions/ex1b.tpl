DATA_SECTION
  init_int NData;
  init_matrix TheData(1,NData,1,2);
  vector Age(1,NData);
  vector Length(1,NData);
  !! Age = column(TheData,1);
  !! Length = column(TheData,2);
  
PARAMETER_SECTION
  init_number logLinf;
  init_number logKappa
  init_number t0
  init_number logSigma;
  number Linf;
  number Sigma;
  number Kappa;
  
  vector Pred(1,NData);
  objective_function_value obj_fun

INITIALIZATION_SECTION
 logLinf 4.78
 logKappa -3
 logSigma 2


PROCEDURE_SECTION
 
 // Extract the parameters
 Linf = exp(logLinf);
 Kappa = exp(logKappa);
 Sigma = exp(logSigma);
 
 // Predictions
 Pred = Linf*(1-exp(-Kappa*(Age-t0)));
 
 // Objective function
 obj_fun = NData*log(Sigma)+ 0.5*norm2(Length-Pred)/square(Sigma);
 cout << Linf << " " << Kappa << " " << t0 << " " << Sigma << " " << obj_fun << endl;
 
REPORT_SECTION
 int II;
 dvariable Iage,Pred;
 report << Linf << " " << Kappa << " " << t0 << endl;
 for (II=1;II<=20;II++)
  {
   Iage = II*1.0;
	   Pred = Linf*(1-exp(-Kappa*(Iage-t0)));
   report << Iage << " " << Pred << endl;
   
  }
