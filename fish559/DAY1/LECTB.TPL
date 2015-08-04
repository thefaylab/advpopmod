DATA_SECTION
 init_int Fyear;
 init_int Lyear;
 init_vector Catch(Fyear,Lyear);
 init_vector CPUE(Fyear,Lyear);
 !! cout << Catch << endl;

PARAMETER_SECTION
 init_number logr;
 init_number logq;
 init_number logK;
 number r;
 number q;
 number K;
 vector Bio(Fyear,Lyear+1);
 
 objective_function_value objn;

INITIALIZATION_SECTION
 logr -0.6
 logq -9
 logK 8.5 

PROCEDURE_SECTION
 int Year;
 dvariable SS;
 
 r = mfexp(logr);
 q = mfexp(logq);
 K = mfexp(logK);

 Bio(Fyear) = K;
 for (Year=Fyear;Year<=Lyear;Year++)
  {
   Bio(Year+1) = Bio(Year) + r*Bio(Year)*(1.0-Bio(Year)/K) - Catch(Year);
   if (Bio(Year+1) < 1) Bio(Year+1) = 1;
  } 
 SS = 0;
 for (Year=Fyear;Year<=Lyear;Year++)
  if (CPUE(Year) > 0)
   SS += square(log(q*Bio(Year)/CPUE(Year)));
 objn = SS;  
 cout << SS << endl;
 

