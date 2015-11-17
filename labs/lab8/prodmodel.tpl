// SCHAEFER MODEL FOR MS_PROD OUTPUT
// SINGLE SPECIES
// GAVIN FAY
// 9 Feb 2012

DATA_SECTION

 init_int Nsp
 init_int rFase
 init_vector r_init(1,Nsp)
 init_int KFase
 init_vector K_init(1,Nsp)
 init_int zFase
 init_vector z_init(1,Nsp)
 init_int thetaFase
 init_vector theta_init(1,Nsp)
 init_int Fyear
 init_int Lyear
 init_matrix ObsCat(Fyear,Lyear,1,Nsp)
 init_int NBio
 !!int ncol=Nsp+1;
 init_matrix ObsBio(1,NBio,1,ncol)
 //init_matrix ObsCV(1,NBio,1,ncol)

 matrix logcpue(1,NBio,1,Nsp)

 int i
 int iyr
 int iobs


PARAMETER_SECTION

 init_number dummy(-1);

 //init_bounded_vector log_r(1,Nsp,-10.1,0.4,rFase)
 init_vector log_r(1,Nsp,rFase) 
 init_bounded_vector log_K(1,Nsp,-10,20.1,KFase)
 //init_bounded_vector log_z(1,Nsp,0.00001,10,zFase)
 init_vector log_z(1,Nsp,zFase)
// init_bounded_vector logit_theta(1,Nsp,-10,0.,thetaFase)
 init_vector logit_theta(1,Nsp,thetaFase)


 number eps

 vector r(1,Nsp)
 vector K(1,Nsp)
 vector z(1,Nsp)
 vector theta(1,Nsp)
 vector q(1,Nsp)
 vector sigma(1,Nsp) 
 number gamma
 number temp

 vector logcpue_pred(Fyear,Lyear)
 vector resid(Fyear,Lyear)
 matrix bio_pred(Fyear,Lyear,1,Nsp)

 sdreport_matrix Bio(Fyear,Lyear+1,1,Nsp)

 objective_function_value objfun

PRELIMINARY_CALCS_SECTION

 log_r = log(r_init);
 log_K = log(K_init);
 log_z = log(z_init-1.);
 //logit_theta = log(elem_div(theta_init,(1.-theta_init)));
 logit_theta = log(theta_init);

PROCEDURE_SECTION

  // constant added to cpue predictions
  eps = 1.e-07;

 objfun += square(dummy);

 //parameter values
 r = mfexp(log_r);
 K = mfexp(log_K);
 z = 1.+eps+mfexp(log_z);
 theta = mfexp(logit_theta); //elem_div(mfexp(logit_theta),(1.+mfexp(logit_theta)));

 //cout << r << endl;
 //cout << K << endl;
 //cout << theta << endl;

 //get_biomass
 Bio(Fyear) = elem_prod(theta,K);
 for (i=1;i<=Nsp;i++)
  {
   gamma = pow(z(i),(z(i)/(z(i)-1)))/(z(i)-1);
 for (iyr=Fyear+1;iyr<=Lyear+1;iyr++)
  {
   dvariable fpen1=0.;
   //dvariable temp = gamma*y(i)*Bio(iyr-1,i)-gamma*y(i)*K(i)*pow(Bio(iyr-1,i)/K(i),z(i));
   //Bio(iyr,i) = posfun(Bio(iyr-1,i)+temp,1.,fpen1);
   Bio(iyr,i) = posfun(Bio(iyr-1,i)*(1.+r(i)*(1.-pow((Bio(iyr-1,i)/K(i)),z(i)-1))/(z(i)-1)),1.,fpen1);
   dvariable sr = 1.-ObsCat(iyr-1,i)/Bio(iyr,i);
   dvariable kcat=ObsCat(iyr-1,i);
   objfun+=1000*fpen1;
   if(sr< 0.001)
    {
     dvariable fpen=0.;
     kcat=Bio(iyr,i)*(1.-posfun(sr,0.001,fpen));
     objfun+=10000*fpen;
     // cout << " kludge "<<iy <<" "<<kcat<<" "<<cat(iy)<<" "<<fpen<<endl;
    }
    Bio(iyr,i)-=kcat;
    bio_pred(iyr-1,i) = 0.5*(Bio(iyr-1,i)+Bio(iyr,i));
 //   cout << i << " " << iyr << " " << Bio(iyr,i) << endl;
  }
 }
 

 //objective function value
 for (i=1;i<=Nsp;i++)
  { 
  // MLE for catchability (Polacheck et al 1993 equation 11)
  q(i) = mfexp(sum(log(elem_div(column(ObsBio,2),eps+column(bio_pred,i)))/NBio));

  // predictions
  logcpue_pred = log(eps+q(i)*column(bio_pred,i));
  // calculate model residuals
  resid = log(column(ObsBio,2))-logcpue_pred;
  
  // MLE for observation error standard deviation (Polacheck et al 1993 equation 10c)
  sigma(i) = sqrt(norm2(resid)/NBio);
  cout << i << " " << q(i) << " " << sigma(i) << endl;
  objfun += NBio*log(sigma(i)) + NBio/2;

  }
  cout << objfun << endl;


REPORT_SECTION

 report << objfun << endl; 
 report << r << endl;
 report << K << endl;
 report << z << endl;
 report << theta << endl;
 report << q << endl;
 report << sigma << endl;
 for (iyr=Fyear;iyr<=Lyear+1;iyr++)
  report << iyr << " " << Bio(iyr) << endl;

 