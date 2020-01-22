#include "rootheader.h"

#ifndef CONSTANTS_H_
#define CONSTANTS_H_
const double PI=TMath::Pi();
const double roo2=TMath::Sqrt(2);
const double pimass=0.13957;
#endif

// References for some fitting functions
// PRL Vol 89 Number 27 December 2002. C. Adler et al. 
// "Coherent rho0 Production in Ultraperipheral Heavy-Ion Collisions."

double fitGausPol1(double *x,double *par);
double fitGausPol2(double *x,double *par);
double fitGausExp(double *x,double *par);

double BreitWpol1(double *x,double *par);
double BreitWpol2(double *x,double *par);
double BreitWexp(double *x,double *par);

double BreitSodPol1(double *x,double *par);
double BreitSodPol2(double *x,double *par);
double BreitSodExp(double *x,double *par);

double modSodPol1(double *x,double *par);
double modSodPol2(double *x,double *par);
double modSodExp(double *x,double *par);

double RossStodPol1(double *x,double *par);
double RossStodPol2(double *x,double *par);
double RossStodExp(double *x,double *par);

// backgound functions
double fPol1(double *x,double *par);
double fPol2(double *x,double *par);
double fExp(double *x, double *par);

// these are the curves with offset for plotting purposes
double fGauss(double *x,double *par);
double fBreitW(double *x,double *par);
double fBreitSod(double *x,double *par);
double fmodSod(double *x,double *par);
double fRossStod(double *x,double *par);

// these are just the curves without offset for calculation purposes
double curveBreitW(double *x,double *par);
double curveBreitSod(double *x,double *par);
double curveModSod(double *x,double *par);
double curveRossStod(double *x,double *par);

double v2cosfit(double *x,double *par);
double expDecay(double *x,double *par);

TF1 *bkgFit(vector<double> &vals,double (*fitType) (double *,double *),int someNum);
TF1 *funcFit(vector<double> &vals,double (*fitType) (double *,double *),int someNum);



