//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Steve Cadrin
//
//  Lab 4: Biomass Dynamics
//  time series observation error model assuming logistic growth
//////////////////////////////////////////////
DATA_SECTION
  init_int nyear;
  init_vector year(1,nyear);  // calendar year
  init_vector yield(1,nyear); // yield (catch biomass)
  init_vector cpue(1,nyear);

PARAMETER_SECTION
  init_number log_r;     // intrinsic growth rate
  init_number log_K;  // carrying capacity (same units as catch)
  init_number log_B1;                 // biomass in year 1 (same units as catch)

  sdreport_number r;     // intrinsic growth rate
  sdreport_number K;  // carrying capacity (same units as catch)
  sdreport_number B1;                 // biomass in year 1 (same units as catch)
  sdreport_number q;
  sdreport_number sigma;

  vector biomass(1,nyear+1);       // time series of biomass estimates (same units as catch)
  vector logcpue(1,nyear);         // Polachek et al. 1993 assume multiplicative observation error in cpue 
  vector logcpue_pred(1,nyear);  // model predictions of log CPUE series
  vector resid(1,nyear)          // log residual
  vector effort(1,nyear)         // effort derived from catch and CPUE
  vector F(1,nyear)              // fishing mortality
  number MSY                     // maximum sustainable yield
  number Bmsy                    // biomass associated with MSY
  number Fmsy                    // F associated with MSY
  vector rel_biomass(1,nyear+1)             // B/Bmsy
  vector rel_F(1,nyear)                   // F/Fmsy
  vector vuln_bio(1,nyear)

  number survival
  number eps

  objective_function_value obj_fun;

PRELIMINARY_CALCS_SECTION
// output data series to screen
  cout << year << endl;
  cout << yield << endl;
  cout << cpue << endl;
// log transform cpue
  logcpue = log(cpue);
// output transformed data
  cout << logcpue << endl;
// derive effort from catch and CPUE
  effort = elem_div(yield, cpue);    //this is a value assigned to a dvar_vector, should occur in the procedure section
// output effort series
  cout << effort << endl;

PROCEDURE_SECTION

  // constant added to cpue predictions
  eps = 1.e-07;

  r = mfexp(log_r);
  K = mfexp(log_K);
  B1 = mfexp(log_B1);  

// begin time series of biomass estimates with estimated parameter B1
// if time series includes the beginning of the fishery, B1~K may be a reasonable assumption
  biomass(1) = B1;

// for loop to derive time series of biomass estimates as a function of biomass in the previous year, catch in the previous year, r, and K
  //cout << r << " " << " " << K << " " << biomass(1) << endl;
  for (int i=1;i<=nyear;i++)
   {
    biomass(i+1) = biomass(i)+r*biomass(i)*(1.-biomass(i)/K);  //biomass prediction (without catch)
    vuln_bio(i) = biomass(i+1); //temporary storage for exploitable biomass
    survival = 1.-yield(i)/biomass(i+1);    //calculate survival from catch,  needs to be > 0, we will use posfun() to achieve this
    dvariable fpen = 0.;  // set tempvariable fpen to 0
    biomass(i+1) *= posfun(survival,0.01,fpen);   // keeps survival above 0.01, if original val is below then large value added to fpen
    obj_fun += 100*fpen;                          // add fpen to objective_function
    F(i) = (vuln_bio(i)-biomass(i+1))/vuln_bio(i);  // calculate F
    //cout << i << " " << biomass(i+1) << " " << obj_fun << endl;
   }


  // observation model

  // MLE for catchability (Polacheck et al 1993 equation 11)
  q = mfexp(sum(log(elem_div(cpue,eps+biomass(1,nyear))))/nyear);

  // predictions
  logcpue_pred = log(eps+q*biomass(1,nyear));
  // calculate model residuals
  resid = logcpue-logcpue_pred;
  
  // MLE for observation error standard deviation (Polacheck et al 1993 equation 10c)
  sigma = sqrt(norm2(resid)/nyear);

  // normal likelihood function (polacheck et al 1993 eqn 10)
  //obj_fun = nyear*logsigma+0.5*norm2(resid)/square(mfexp(logsigma));
  obj_fun = nyear*log(sigma) + 0.5*nyear;
  cout << obj_fun << endl;


  //biological reference points
  // derive MSY reference points from logistic growth parameters
  MSY =  r*K/4;
  Fmsy = r/2;
  Bmsy = K/2;
  rel_biomass = biomass/Bmsy;
  rel_F = F/Fmsy;


REPORT_SECTION
//  report model estimates
  report << "CPUE residuals" << endl;
  report << resid << endl;
  report << "Biomass estimates" << endl;
  report << biomass << endl;
  report << "Fishing mortality estimates" << endl;
  report << F << endl;
  report << "MSY" << endl;
  report << MSY << endl;
  report << "Bmsy" << endl;
  report << Bmsy << endl;
  report << "Fmsy" << endl;
  report << Fmsy << endl;
  report << "Relative biomass estimates (B/Bmsy)" << endl;
  report << rel_biomass << endl;
  report << "Relative F estimates (F/Bmsy)" << endl;
  report << rel_F << endl;



