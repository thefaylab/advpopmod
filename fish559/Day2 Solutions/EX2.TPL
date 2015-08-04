DATA_SECTION
  init_int Nprey
  init_int NData
  init_matrix N(1,NData,1,7)
  
  
  !! ad_comm::change_datafile_name("EX2.CTL");
  init_int Model_Num
  
  int PhaseBeta;
  !! PhaseBeta = 1;
  !! if (Model_Num == 1) PhaseBeta = -1;
  int PhaseGamma;
  !! PhaseGamma = 1;
  !! if (Model_Num == 1 | Model_Num == 2) PhaseGamma = -1;
  
PARAMETER_SECTION
  init_bounded_vector alpha(1,Nprey,0,1000)
  init_bounded_vector beta(1,Nprey,0,1000,PhaseBeta)
  init_bounded_vector gamma(1,Nprey,0,5,PhaseGamma)

  matrix Pred(1,NData,1,Nprey) 

  objective_function_value obj_fun

PROCEDURE_SECTION
  int Idata,Iprey;
  
  obj_fun = 0;
  for (Idata=1;Idata<=NData;Idata++)
   for (Iprey=1;Iprey<=Nprey;Iprey++)
    {
     // Model 1
     if (Model_Num == 1) Pred(Idata,Iprey) = alpha(Iprey)*N(Idata,1);
     
     // Model 2
     if (Model_Num == 2) Pred(Idata,Iprey) = alpha(Iprey)*N(Idata,1)/(1.0+beta(Iprey)*N(Idata,1+Iprey));
     
     // Model 3
     if (Model_Num == 3) Pred(Idata,Iprey) = alpha(Iprey)*N(Idata,1)*pow(N(Idata,1+Iprey),gamma(Iprey)-1.0)/(1.0+beta(Iprey)*pow(N(Idata,1+Iprey),gamma(Iprey)));
     
     // Model 4
     if (Model_Num == 4) Pred(Idata,Iprey) = alpha(Iprey)*N(Idata,1)/(1+beta(Iprey)*N(Idata,1+Iprey)+gamma(Iprey)*N(Idata,1));
     
     // Residuals
     obj_fun = obj_fun + square(log(Pred(Idata,Iprey))-log(N(Idata,4+Iprey)));
    }

//===================================================================================================

REPORT_SECTION
 int Ipred, IpreyV,Iprey,Nprey;
 dvariable PredOut;
 dvar_vector PredDens(1,3);
 dvar_vector PreyDens(0,100);
 PredDens(1) = 0.5;
 PredDens(2) = 1.0;
 PredDens(3) = 2.0;

 cout << "IN REPORT" << endl; 
 Nprey = 100;
 for (Iprey=0;Iprey<=Nprey;Iprey++) PreyDens(Iprey) = float(Iprey)*0.02;
 
 for (Ipred=1;Ipred<=3;Ipred++)
  for (IpreyV=1;IpreyV<=3;IpreyV++)
   for (Iprey=0;Iprey<=Nprey;Iprey++)
    {
          
     // Model 1
     if (Model_Num == 1) PredOut = alpha(IpreyV)*PredDens(Ipred);
     
     // Model 2
     if (Model_Num == 2) PredOut = alpha(IpreyV)*PredDens(Ipred)/(1.0+beta(IpreyV)*PreyDens(Iprey));
     
     // Model 3
     if (Model_Num == 3) PredOut = alpha(IpreyV)*PredDens(Ipred)*pow(PreyDens(Iprey),gamma(IpreyV)-1.0)/(1.0+beta(IpreyV)*pow(PreyDens(Iprey),gamma(IpreyV)));
     
     // Model 4
     if (Model_Num == 4) PredOut = alpha(IpreyV)*PredDens(Ipred)/(1+beta(IpreyV)*PreyDens(Iprey)+gamma(IpreyV)*PredDens(Ipred));
     report << Ipred << " " << IpreyV << " " << Iprey << " " << PredOut << " " << PredDens(Ipred) << " " << PreyDens(Iprey) << endl;
     
    }

 cout << "OUT REPORT" << endl; 
