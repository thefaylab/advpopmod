DATA_SECTION
 init_int Nyear;                                 // Number of years
 init_int Nage;                                  // Last age-class
 init_number Mval;                               // Natural mortality
 init_vector Wght(0,Nage);                       // Weight-at-age
 init_number SigCatch;                           // CV for catch data
 init_number SigCPUE;                            // Sigma for log(CPUE)
 init_number Omega;                              // Weight for proportion data
 init_matrix CatchCpue(1,Nyear,0,2)              // Catch and cpue data
 init_matrix Propn(1,Nyear,-1,Nage)              // Proportion data
 vector Catch(1,Nyear);                          // Catch data
 vector CPUE(1,Nyear);                           // CPUE data
 
 !! Catch = column(CatchCpue,1);                 // Extract the catch data
 !! CPUE = column(CatchCpue,2);                  // Extract the CPUE data
 
// !! cout << "Done Data" << endl;

INITIALIZATION_SECTION
 logN 3
 Sel50 3
 Sel95 8
 logFish -1
 logq -1

PARAMETER_SECTION
 init_vector logN(1,Nyear+Nage,1)                // Log of recruits
 init_bounded_number Sel50(0,Nage,1)             // Age-at-50%-selectivity
 init_bounded_number Sel95(0,Nage,1)             // Age-at-95%-selectivity
 init_vector logFish(1,Nyear,1)                  // Log of fishing mortality
 init_number logq(1)                             // Log of catchability
 init_number dummy(10);

 matrix N(1,Nyear+1,0,Nage)                      // Numbers-at-age
 vector S(0,Nage)                                // Selectivity-at-age
 matrix F(1,Nyear,0,Nage)                        // F-matrix
 matrix Z(1,Nyear,0,Nage)                        // Z-matrix
 
 vector CPred(1,Nyear);                          // Predicted catch (weight);
 vector CPUEPred(1,Nyear);                       // Predicted CPUE (weight);
 matrix PropnPred(1,Nyear,0,Nage);               // Predicted catch-at-age
 vector Bio(1,Nyear);
 
 number Like1;                                   // Catch data
 number Like2;                                   // CPUE data
 number Like3;                                   // Proportion data
 objective_function_value objn;

 sdreport_number Ratio1;                         // Output statistic
 likeprof_number Ratio3;                         // Output statistic
 likeprof_number Ratio2;                         // Output statistic
 likeprof_number Sel1;

 number weight;
 !! weight = 0.0001;
 
 
// !! cout << "Done Parameters" << endl;
 
// ==========================================================================

PROCEDURE_SECTION
 int Year;
 dvariable SS;

 Bratio = Bio(20)/Bio(1);

 // Set up the selectivity pattern
 Select();

 // Project the model forward and compute various outputs
 Numbers();
 
 // Compute the likelihood
 Likelihood();

 objn = Like1+Like2+Like3;
 objn += weight*square(Bio(Nyear)/Bio(1)-0.1);
 cout << weight << endl;
  cout << objn << endl;
 
 // Statistic to summarize
 Ratio1 = Bio(Nyear)/Bio(1);
 Ratio2 = Bio(Nyear)/Bio(1);
 Ratio3 = mfexp(logFish(Nyear))/mfexp(logFish(1));
 Sel1 = Sel50;

 if (mceval_phase())
  {
   cout << objn << " " << Ratio1 << " " << Sel50 << " " << Sel95 << " " << logq << " " << logN << " "  << logFish << endl;
  }


// ---------------------------------------------------------------------------- 

FUNCTION Select
 int Age;
 
 for (Age=0;Age<=Nage;Age++)
  {
   S(Age) = 1.0/(1+exp(-log(19)*(Age-Sel50)/(Sel95-Sel50)));
  }

// ---------------------------------------------------------------------------- 

FUNCTION Numbers
 int Age,Year;
 
 // Clear the N matrix
 N.initialize();

 // Compute the F matrix
 for (Year=1;Year<=Nyear;Year++)
  for (Age=0;Age<=Nage;Age++)
   F(Year,Age) = mfexp(logFish(Year))*S(Age);
 Z = F + Mval;  

 // Insert the recruitments
 for (Age=0;Age<=Nage;Age++)
  N(1,Age) = mfexp(logN(Nage-Age+1));
 for (Year=2;Year<=Nyear;Year++)
  N(Year,0) = mfexp(logN(Nage+Year)); 

 // Project the whole N matrix
 for (Year=1;Year<=Nyear;Year++)
  for (Age=0;Age<Nage;Age++)
   N(Year+1,Age+1) = N(Year,Age)*mfexp(-Z(Year,Age));
   
 // Compute the predicted exploitable biomass, catch-at-age and catch
 PropnPred.initialize();
 for (Year=1;Year<=Nyear;Year++)
  {
   Bio(Year) = 0; CPred(Year) = 0;
   for (Age=0;Age<=Nage;Age++)
    {
     PropnPred(Year,Age) = F(Year,Age)/Z(Year,Age)*N(Year,Age)*(1.0-mfexp(-Z(Year,Age)));
     CPred(Year) += Wght(Age)*PropnPred(Year,Age);
     Bio(Year) += Wght(Age)*S(Age)*N(Year,Age)*mfexp(-Z(Year,Age)/2.0);
    } 
   CPUEPred(Year) = mfexp(logq)*Bio(Year); 
   PropnPred(Year) /= sum(PropnPred(Year));
  }

// ---------------------------------------------------------------------------- 

FUNCTION Likelihood
 int Year,Age;
 
 // Catch data
 Like1 = 0;
 for (Year=1;Year<=Nyear;Year++)
  Like1 += square( (Catch(Year)-CPred(Year))/CPred(Year));
 Like1 = Like1 / (2.0*square(SigCatch));
 
 // CPUE data
 Like2 = 0;
 for (Year=1;Year<=Nyear;Year++)
  Like2 += square( log(CPUE(Year)) - log(CPUEPred(Year)));
 Like2 = Like2 / (2.0*square(SigCPUE)); 
 
 // Catch-at-age data
 Like3 = 0;
 for (Year=1;Year<=Nyear;Year++)
  for (Age=0;Age<=Nage;Age++)
   if (Propn(Year,Age) >0)
    Like3 += Propn(Year,Age)*log(PropnPred(Year,Age)/Propn(Year,Age));
 Like3 = -1*Omega*Like3;  
 
// ==========================================================================

REPORT_SECTION
 report << objn << " " << Like1 << " " << Like2 << " " << Like3 << endl;
 report << S << endl;
 report << mfexp(logFish) << endl;
 report << N << endl;
 report << PropnPred << endl;
 report << CPred << endl;
 report << CPUEPred << endl;
 
 weight = weight * 10;


