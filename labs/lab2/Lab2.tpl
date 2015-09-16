//////////////////////////////////////////////
//  MAR 580: Advanced Population Modeling
//  Fall 2015
//  Gavin Fay
//  September 8, 2015
//
//  Lab 2: Mark Recapture
//////////////////////////////////////////////
DATA_SECTION
  init_int nphi;   //number of phi params
  init_int np;     // number of p params
  init_int ndata;   //number of data points, not really used but model is generalized
  init_vector recaps_111(1,ndata);
  init_vector recaps_110(1,ndata);
  init_vector recaps_101(1,ndata);
  init_vector recaps_100(1,ndata);

  int i;

PARAMETER_SECTION
  init_vector logit_phi(1,nphi);   //estimating parameters in logit space
  init_vector logit_p(1,np);    //estimating parameters in logit space
  
  sdreport_vector phi(1,2);
  sdreport_vector p(1,2);

  vector p_111(1,ndata);
  vector p_110(1,ndata);
  vector p_101(1,ndata);
  vector p_100(1,ndata);

  objective_function_value obj_fun;

PROCEDURE_SECTION

  //transform the estimated params to the 'real' parameters
  //cout << logit_phi << endl;
  //cout << logit_p << endl;
  for (i=1;i<=2;i++)
   {
    phi(i) = mfexp(logit_phi(1))/(1.+mfexp(logit_phi(1)));
    p(i) = mfexp(logit_p(1))/(1.+mfexp(logit_p(1))); 
   }
  if (nphi==2) phi(2) = mfexp(logit_phi(2))/(1.+mfexp(logit_phi(2)));
  if (np==2) p(2) = mfexp(logit_p(2))/(1.+mfexp(logit_p(2)));
  

  // get the probabilities for each recapture history
  p_111 = phi(1)*p(1)*phi(2)*p(2);
  p_110 = phi(1)*p(1)*(1.-phi(2)*p(2));
  p_101 = phi(1)*(1.-p(1))*phi(2)*p(2);
  p_100 = 1.-p_111-p_110-p_101;

  //obj_fun = -(a*log(phi*p*phi*p)+b*log(phi*p*(1-phi*p))+c*log(phi*(1-p)*phi*p)+d*log(1-phi*p-phi*(1-p)*phi*p));

  //calculate the objective function
  obj_fun -= sum(elem_prod(recaps_111,log(p_111)));
  obj_fun -= sum(elem_prod(recaps_110,log(p_110)));
  obj_fun -= sum(elem_prod(recaps_101,log(p_101)));
  obj_fun -= sum(elem_prod(recaps_100,log(p_100)));

REPORT_SECTION
  report << "phi" << endl;
  report << phi << endl;
  report << "p" << endl;
  report << p << endl;

