DATA_SECTION
  init_int NData;
  init_matrix TheData(1,NData,1,2);
  vector Age(1,NData);
  vector Length(1,NData);
  !! Age = column(TheData,1);
  !! Length = column(TheData,2);
  
PARAMETER_SECTION
  init_number logLinf;
  init_number loga50;
  init_number logDelta;
  init_number logSigma;
//  init_number logKappa
//  init_number t0
  number Linf;
  number a50;
  number Delta;
  number Sigma;
  number Kappa;
  
  vector Pred(1,NData);
  objective_function_value obj_fun

INITIALIZATION_SECTION
 logLinf 4.78
 loga50 2.30
 logDelta 2.1
 

PROCEDURE_SECTION
 
 // Extract the parameters
 Linf = exp(logLinf);
 a50 = exp(loga50);
 Delta = exp(logDelta);
 Sigma = exp(logSigma);
 
 // Predictions
 Pred = Linf/(1+exp(-log(19)*(Age-a50)/Delta));
 
 // Objective function
 obj_fun = NData*log(Sigma)+ 0.5*norm2(Length-Pred)/square(Sigma);
 cout << Linf << " " << a50 << " " << Delta << " " << Sigma << " " << obj_fun << endl;
 
REPORT_SECTION
 int II;
 dvariable Iage,Pred;
 report << Linf << " " << a50 << " " << Delta << endl;
 for (II=1;II<=20;II++)
  {
   Iage = II*1.0;
   Pred = Linf/(1+exp(-log(19)*(Iage-a50)/Delta));
   report << Iage << " " << Pred << endl;
   
  }
