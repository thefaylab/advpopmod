DATA_SECTION
 init_int Nyear;                                              // Number of years
 init_int Nclass;                                             // Number of size-classes
 init_vector Length(1,Nclass)                                 // lengths
 init_vector WL(1,Nclass)                                     // weights
 init_matrix X(1,Nclass,1,Nclass)                             // Size-transition matrix
 init_number S50;                                             // Length-at-50% selectivity (fishery)  
 init_number S95;                                             // Length-at-95% selectivity (fishery)  
 init_number SS50;                                            // Length-at-50% selectivity (survey)  
 init_number SS95;                                            // Length-at-95% selectivity (survey)  
 init_number M;                                               // Natural mortality
 init_vector CWObs(1,Nyear);                                  // Catches
 init_matrix CALObs(1,Nyear,1,Nclass);                        // Observed CAL
 init_number Neff;                                            // Effectuve sample size
 init_vector BioIndex(1,Nyear);                               // Biomass index  
 init_number BioSig;                                          // Sigma for biomass
 
 init_int Check;
 !!if (Check != 123456) { cout << "Input file in error" << endl; exit(1); }

PARAMETER_SECTION
 init_number LogRbar(1);                                      // Mean recruitment  
 init_vector LogNinit(1,Nclass,1)                             // Initial Ns
 init_bounded_vector logFullF(1,Nyear,-20,5)                  // Fully-selected F
 init_vector Eps(1,Nyear)                                     // Epsilon

 matrix N(1,Nyear+1,1,Nclass);                                // Numbers-at-length
 matrix F(1,Nyear,1,Nclass);                                  // Fishing mortality-at-length
 matrix Z(1,Nyear,1,Nclass);                                  // Total mortality-at-length
 matrix CAL(1,Nyear,1,Nclass);                                // catch-at-length
 vector CW(1,Nyear);                                          // catch-in-weight
 vector BioPred(1,Nyear);                                     // Predicted biomass
 vector S(1,Nclass);
 vector SurveyS(1,Nclass);

 number Penal;                                                // Penalty on rec-devs
 number LikeCatch;                                            // catches
 number LikeBio;                                              // Relative index of abundance
 number LikeCAL;                                              // Catch-at-length
 sdreport_number Rbar;
 
 objective_function_value obj_fun;

PRELIMINARY_CALCS_SECTION
 int Iyear,Ilen;
 float Total;
 
 // Selectivity

 // Catch-at-length
 for (Iyear=1;Iyear<=Nyear;Iyear++)
  {
   Total = 0;
   for (Ilen=1;Ilen<=Nclass;Ilen++) Total += CALObs(Iyear,Ilen);
   for (Ilen=1;Ilen<=Nclass;Ilen++) CALObs(Iyear,Ilen) = CALObs(Iyear,Ilen)/ Total;
  }
 
PROCEDURE_SECTION
 Rbar = mfexp(LogRbar);

 CAL.initialize();
 CW.initialize();

 // Creat F and Z matrices
 Get_FZ_matrix();
 
 // Create N matrix
 Get_N_matrix();

 // Get the loglikelihood
 LogLikelihood();

 obj_fun = LikeCatch + LikeCAL + LikeBio + Penal;
 //cout << obj_fun << " " << LikeCatch << " " << LikeCAL << " " << LikeBio << " " << Penal << endl;

 if (mceval_phase())
  cout << BioPred << endl;


// =========================================================================================================

FUNCTION Get_FZ_matrix                                    
 int Ilen,Iyear;
 


// =========================================================================================================

FUNCTION Get_N_matrix
 int Ilen,Jlen,Iyear;
 dvariable CALtot;

 for (Ilen=1;Ilen<=Nclass;Ilen++) N(1,Ilen) = exp(LogNinit(Ilen));
 
 for (Iyear=1;Iyear<=Nyear;Iyear++)
  {

    
  
  }
 
// =========================================================================================================

FUNCTION LogLikelihood
 int Iyear,Ilen; 
 dvariable SS,Top,Bot,q;
  

REPORT_SECTION
 int Iyear;
 
 for (Iyear =1;Iyear<=Nyear;Iyear++)
  report << Iyear << " " << BioIndex(Iyear) << " " << BioPred(Iyear)  << endl;
  report << N << endl;
 report << CALObs << endl;
 report << CAL << endl;
 