DATA_SECTION
  init_int Nline
  init_int DatCol;
  init_matrix TheData(1,Nline,1,3)
  !! cout << TheData << endl;
  vector Years(1,Nline);
  vector Cites(1,Nline);
  !! Years = column(TheData,1);
  !! Cites = column(TheData,DatCol);
  
PARAMETER_SECTION
  init_number intercept;
  init_number slope;
  init_bounded_number scale(-500,500,2);
  init_number offset(3);
  init_bounded_number period(0,200,2);

  vector Pred(1,Nline) 

  objective_function_value obj_fun

PROCEDURE_SECTION
  int Idata;
  dvariable Term2;
  
  obj_fun = 0;
  for (Idata=1;Idata<=Nline;Idata++)
   {
    Term2 = (Years(Idata)-Years(1)-offset*period)/period;
    Pred(Idata) = exp(intercept + slope*(Years(Idata)-Years(1))+scale*sin(Term2));
//    obj_fun = obj_fun + Pred(Idata) - Cites(Idata)*log(Pred(Idata));
    obj_fun = obj_fun - log_density_poisson(Cites(Idata),Pred(Idata));
   } 
  cout << obj_fun << endl; 

//===================================================================================================

RUNTIME_SECTION
 maximum_function_evaluations 1000,10000,100000

REPORT_SECTION
 int Idata; 
 for (Idata=1;Idata<=Nline;Idata++)
  report << Years(Idata) << " " << Cites(Idata) << " " << Pred(Idata) << endl;

 