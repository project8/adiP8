#include "TMath.h"

Double_t djdx(Int_t n, Double_t x)
{
  Double_t dj = TMath::BesselI(n, x) - TMath::BesselI(n, x - 1e-6);
  return dj / (1e-6);
}


Double_t s(Double_t w, Int_t n, Double_t beta0, Double_t beta1, Double_t theta)
{
  Double_t qe = 1.602e-19;
  Double_t w00 = 100e-9;        //100 rad-GHz
  Double_t w0 = w00 / sqrt(1 - beta0 * beta0);
  Double_t pi = TMath::Pi();
  Double_t c = 3e8;             //meters/sec  

//   return qe*qe*beta0*beta0*w0*w0*n*n/(2*TMath:Power(pi,3),c)*
//     TMath::Power((sin(theta) - beta1)/(beta0*cos(theta)),2)*
//     TMath::Power(TMath::BesselI(n,n*beta0*cos(theta)),2) + 
//     TMath::Power(djdx(n,n*beta0*cos(theta)),2) / 
//     (TMath::Power(w*(1-beta1*sin(theta)) - n*w00*sqrt(1-beta0*beta0-beta1*beta1),2)); 

  return 1.0 / (TMath::Power(w * (1 - beta1 * sin(theta)) - n * w00 * sqrt(1 - beta0 * beta0 - beta1 * beta1), 2));


}



Double_t ptot(Double_t b, Double_t beta0, Double_t theta_pitch)
{
  Double_t betat = beta0 * sin(theta_pitch);  //parallel v
  Double_t e = 1.602e-19;
  Double_t pi = TMath::Pi();
  Double_t c = 3e8;             //meters/sec
  Double_t m = 9.109e-31;       //kg
  Double_t w = b * e / m;       // rad-Hz

  return 1.0 / (4 * pi * 8.854e-12) * 2 * e * e * w * w / (3 * c) * betat * betat / (1 - beta0 * beta0);
}

Double_t w(Double_t b, Double_t beta0)
{
  Double_t e = 1.602e-19;
  Double_t pi = TMath::Pi();
  Double_t c = 3e8;             //meters/sec
  Double_t m = 9.109e-31;       //kg
  Double_t w = b * e / (m / sqrt(1 - beta0 * beta0));  // rad-Hz

  return w;
}

Double_t dpdd(Double_t w, Double_t betat, Double_t theta_pitch, Double_t theta)
{
  Double_t e = 1.602e-19;
  Double_t pi = TMath::Pi();
  Double_t c = 3e8;             //meters/sec
  Double_t beta0 = betat * cos(theta_pitch);  //parallel v
  Double_t betap = betat * sin(theta_pitch);  //perpendicular velocity
  Double_t gp = 1 - betap * cos(theta);
  return ptot(w, beta0, betat) / (4 * pi) * (3.0 / 4.0) * TMath::Power((1 - betat * betat), 2) * (4 * gp * gp * ((1 + betap * betap) * (1 + TMath::Power(cos(theta), 2)) - 4 * betap * cos(theta))
                                                                                                  - (1 - betap * betap + 3 * beta0 * beta0) * beta0 * beta0 * TMath::Power(sin(theta), 4)) /
      4 * TMath::Power(gp * gp - TMath::Power(beta0 * sin(theta), 2), 7.0 / 2.0);
}

Double_t sq(Double_t x)         // squared
{
  return TMath::Power(x, 2);
}

Double_t dpdd2(Double_t b, Double_t beta, Double_t theta_pitch, Double_t th)
{

  const Double_t e = 1.602e-19; //coulomb
  const Double_t m = 9.109e-31; //kg
  const Double_t w = b * e / m; // rad-Hz
  const Double_t pi = TMath::Pi();
  const Double_t oo4pe = 1 / (4 * pi * 8.85e-12);  // 1/4 pi epsilon0 = m/F
  const Double_t c = 3e8;       //meters/sec
  const Double_t bpar = beta * cos(theta_pitch);  //parallel v
  const Double_t bperp = beta * sin(theta_pitch);  //perpendicular velocity
  const Double_t gpar = 1 - bpar * cos(th);

  return oo4pe * e * e * w * w / (8 * pi * c) * (1 - sq(beta)) * (sq(bperp)) * (4 * sq(gpar) * ((1 + sq(bpar)) * (1 + sq(cos(th)))
                                                                                                - 4 * bpar * cos(th)) - (1 - sq(bpar) + 3 * sq(bperp)) * sq(bperp) * TMath::Power(sin(th), 4)) /
      (4 * TMath::Power(sq(gpar) - sq(bperp * sin(th)), 7.0 / 2.0));
}
