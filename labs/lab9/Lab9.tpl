//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Steve Cadrin
//  November 5, 2015
//
//  Lab 8: Multi-State Mark Recapture 
//////////////////////////////////////////////
DATA_SECTION
  init_int A00;     // frequency of encounter history
  init_int A0A;     // frequency of encounter history
  init_int A0B;     // frequency of encounter history
  init_int A0C;     // frequency of encounter history
  init_int AA0;     // frequency of encounter history

PARAMETER_SECTION
  init_number logit_phi;   //estimating survival in logit space
  init_number logit_p;     //estimating capture probability in logit space
  init_number logit_psi;   //estimating movement from A probability in logit space
  sdreport_number phi;
  sdreport_number p;
  sdreport_number psi;
  number p_A00;
  number p_A0A;
  number p_A0B;
  number p_A0C;
  number p_AA0;
  number exp_frq_A00;
  number exp_frq_A0A;
  number exp_frq_A0B;
  number exp_frq_A0C;
  number exp_frq_AA0;

  objective_function_value obj_fun;

PROCEDURE_SECTION
  //transform the estimated params to the 'real' parameters
  //cout << logit_phi << endl;
  //cout << logit_p << endl;
  //cout << logit_psi << endl;
  phi = mfexp(logit_phi)/(1.+mfexp(logit_phi));
  p = mfexp(logit_p)/(1.+mfexp(logit_p)); 
  psi = mfexp(logit_psi)/(1.+mfexp(logit_psi)); 
  
  // probabilities for each recapture history
  p_A00 = 1-phi*p-phi*(1-p)*phi*p;
  // probability of no recaptures from simple mark-recapture with one state (see lab 2)
  p_A0A = phi*(1-p)*phi*p*(psi*psi+(1-psi)*(1-psi));
  // survived_t1 and not_captured_t2 and survived_t2 and captured_t3 and ((moved_t1 and moved_t2)or (stayed_t1 and stayed_t2)     
  p_A0B = phi*(1-p)*phi*p*(psi*(1-psi)+(1-psi)*(psi));
  // survived_t1 and not_captured_t2 and survived_t2 and captured_t3 and ((moved_t1 and stayed_t2)or (stayed_t1 and moved_t2)     
  p_A0C = phi*(1-p)*phi*p*(psi*(1-psi)+(1-psi)*(psi));
  // survived_t1 and not_captured_t2 and survived_t2 and captured_t3 and ((moved_t1 and stayed_t2)or (stayed_t1 and moved_t2)     
  p_AA0 = phi*psi*p*((1-phi)+(1-p));
  // survived_t1 and stayed_t1 and captured_t2 and (died_t2 or not_captured_t3)
  
  //calculate the objective function
  obj_fun = -A00*log(p_A00)-A0A*log(p_A0A)-A0B*log(p_A0B)-A0C*log(p_A0C)-AA0*log(p_AA0);

  exp_frq_A00=p_A00*(A00+A0A+A0B+A0C+AA0);
  exp_frq_A0A=p_A0A*(A00+A0A+A0B+A0C+AA0);
  exp_frq_A0B=p_A0B*(A00+A0A+A0B+A0C+AA0);
  exp_frq_A0C=p_A0C*(A00+A0A+A0B+A0C+AA0);
  exp_frq_AA0=p_AA0*(A00+A0A+A0B+A0C+AA0);
  
REPORT_SECTION
  report << "phi" << endl;
  report << phi << endl;
  report << "p" << endl;
  report << p << endl;
  report << "psi" << endl;
  report << psi << endl;
  report << "probability of A00" << endl;
  report << p_A00 << endl;
  report << "probability of A0A" << endl;
  report << p_A0A << endl;
  report << "probability of A0B" << endl;
  report << p_A0B << endl;
  report << "probability of A0C" << endl;
  report << p_A0C << endl;
  report << "probability of AA0" << endl;
  report << p_AA0 << endl;
  report << "expected frequency A00" << endl;
  report << exp_frq_A00 << endl;
  report << "expected frequency A0A" << endl;
  report << exp_frq_A0A << endl;
  report << "expected frequency A0B" << endl;
  report << exp_frq_A0B << endl;
  report << "expected frequency A0C" << endl;
  report << exp_frq_A0C << endl;
  report << "expected frequency AA0" << endl;
  report << exp_frq_AA0 << endl;

