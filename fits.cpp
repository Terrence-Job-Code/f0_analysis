#include "rootheader.h"
#include "fitvalues.h"
#include "fits.h"

//#ifndef CONSTANTS_H_
//#define CONSTANTS_H_
//const double PI=TMath::Pi();
//const double roo2=TMath::Sqrt(2);
//const double pimass=0.13957;
//#endif

#ifndef FITFUNCTIONS_H_
#define FITFUNCTIONS_H_

// need to make it so par[0] gives the amplitude on all of them

// gaussian fit is incorrect at the moment, par[2] needs to be out of the sqrt
double fitGausPol1(double *x,double *par){
    // par[0] is the area
    // par[1] is the mass
    // par[2] is the width
    double arg=0;
    if(par[2]!=0) arg = (x[0] - par[1])/par[2];
    double gauss = (par[0]/(TMath::Sqrt(2*PI)*par[2]))*TMath::Exp(-0.5*arg*arg); 
    double fitval = gauss + par[3] + par[4]*x[0];
    return fitval; 
}
double fitGausPol2(double *x,double *par){
    double arg=0;
    if(par[2]!=0) arg = (x[0] - par[1])/par[2];
    double gauss = (par[0]/(TMath::Sqrt(2*PI)*par[2]))*TMath::Exp(-0.5*arg*arg); 
    double fitval = gauss + par[3] + par[4]*x[0] + par[5]*x[0]*x[0];
    return fitval; 
}
double fitGausExp(double *x,double *par){
    double arg=0;
    if(par[2]!=0) arg = (x[0] - par[1])/par[2];
    double gauss = (par[0]/(TMath::Sqrt(2*PI)*par[2]))*TMath::Exp(-0.5*arg*arg); 
    double fitval = gauss + par[3] + par[4]*TMath::Exp(-par[5]*x[0]);
    return fitval; 
}
double BreitWpol1(double *x,double *par){
    // par[0] = Amplitude
    // par[1] = mass
    // par[2] = width

    double gam = TMath::Sqrt(par[1]*par[1]*(par[1]*par[1]+par[2]*par[2]));
    double k=(2*roo2*par[1]*par[2]*gam)/(PI*TMath::Sqrt(par[1]*par[1]+gam));
    double sim1=x[0]*x[0] - par[1]*par[1];
    double sim2=par[1]*par[2];

    double f = (par[0]*k)/(sim1*sim1 + sim2*sim2) + par[3] + par[4]*x[0];

    return f;
}
double BreitWpol2(double *x,double *par){

    double gam = TMath::Sqrt(par[1]*par[1]*(par[1]*par[1]+par[2]*par[2]));
    double k=(2*roo2*par[1]*par[2]*gam)/(PI*TMath::Sqrt(par[1]*par[1]+gam));
    double sim1=x[0]*x[0] - par[1]*par[1];
    double sim2=par[1]*par[2];

    double f = (par[0]*k)/(sim1*sim1 + sim2*sim2) + par[3] + par[4]*x[0] + par[5]*x[0]*x[0];

    return f;
   
}
double BreitWexp(double *x,double *par){
    
    double gam = TMath::Sqrt(par[1]*par[1]*(par[1]*par[1]+par[2]*par[2]));
    double k=(2*roo2*par[1]*par[2]*gam)/(PI*TMath::Sqrt(par[1]*par[1]+gam));
    double sim1=x[0]*x[0] - par[1]*par[1];
    double sim2=par[1]*par[2];

    double f = (par[0]*k)/(sim1*sim1 + sim2*sim2) + par[3] + par[4]*TMath::Exp(-par[5]*x[0]);

    return f;

}
double BreitSodPol1(double *x,double *par){
    // par[0] is amplitude of breit term
    // par[1] is f0 mass
    // par[2] is f0 width
    // par[3] is amplitude of sod term
    double term1 = x[0]*par[1]*par[2];
    double term2 = par[1]*par[1]-x[0]*x[0];
    double term3 = par[1]*par[2];
    // BW + Sod + pol
    double f = (par[0]*term1)/(term2*term2 + term3*term3) + (par[3]*term2)/(term2*term2 + term3*term3) + par[4] + par[5]*x[0];
    return f;
}
double BreitSodPol2(double *x,double *par){
    // par[0] is amplitude of breit term
    // par[1] is f0 mass
    // par[2] is f0 width
    // par[3] is amplitude of sod term
    double term1 = x[0]*par[1]*par[2];
    double term2 = par[1]*par[1]-x[0]*x[0];
    double term3 = par[1]*par[2];
    // BW + Sod + pol
    double f = (par[0]*term1)/(term2*term2 + term3*term3) + (par[3]*term2)/(term2*term2 + term3*term3) + par[4] + par[5]*x[0] + par[6]*x[0]*x[0];
    return f;
}
double BreitSodExp(double *x,double *par){
    // par[0] is amplitude of breit term
    // par[1] is f0 mass
    // par[2] is f0 width
    // par[3] is amplitude of sod term
    double term1 = x[0]*par[1]*par[2];
    double term2 = par[1]*par[1]-x[0]*x[0];
    double term3 = par[1]*par[2];
    // BW + Sod + pol
    double f = (par[0]*term1)/(term2*term2 + term3*term3) + (par[3]*term2)/(term2*term2 + term3*term3) + par[4] + par[5]*x[0];
    return f;
}
double modSodPol1(double *x,double *par){
    // par[0] is first term amplitude (amp of "A")
    // par[1] is mass
    // par[2] is width
    // par[3] is interference term (amp of "B") : check reference
    double a = par[0]*TMath::Sqrt(x[0]*par[1]*par[2]);
    double c = x[0]*x[0] - par[1]*par[1];
    double d = par[1]*par[2];

    double f = (a*a + 2*a*c*par[3])/(c*c + d*d) + par[3]*par[3] + par[4] + par[5]*x[0];

    return f;
}
double modSodPol2(double *x,double *par){
    // par[0] is first term amplitude (amp of "A")
    // par[1] is mass
    // par[2] is width
    // par[3] is interference term (amp of "B") : check reference
    double a = par[0]*TMath::Sqrt(x[0]*par[1]*par[2]);
    double c = x[0]*x[0] - par[1]*par[1];
    double d = par[1]*par[2];

    double f = (a*a + 2*a*c*par[3])/(c*c + d*d) + par[3]*par[3] + par[4] + par[5]*x[0] + par[6]*x[0]*x[0];

    return f;
}
double modSodExp(double *x,double *par){
    // par[0] is first term amplitude (amp of "A")
    // par[1] is mass
    // par[2] is width
    // par[3] is interference term (amp of "B") : check reference
    double a = par[0]*TMath::Sqrt(x[0]*par[1]*par[2]);
    double c = x[0]*x[0] - par[1]*par[1];
    double d = par[1]*par[2];

    double f = (a*a + 2*a*c*par[3])/(c*c + d*d) + par[3]*par[3] + par[4] + par[5]*TMath::Exp(-par[6]*x[0]);

    return f;
}
double RossStodPol1(double *x,double *par){
    // par[0] is amplitude
    // par[1] is mass
    // par[2] is width
    // par[3] is n
    double term1 = x[0]*x[0] - 4*pimass*pimass;
    double term2 = par[1]*par[1]-4*pimass*pimass;
    if(term2==0) term2=-9999; // want a different error checker
    double term3 = term1/term2;
    double gamma = par[2]*(par[1]/x[0])*TMath::Power(term3,1.5);

    double term4 = x[0]*par[1]*gamma;
    double term5 = par[1]*par[1] - x[0]*x[0];
    double term6 = par[1]*gamma;

    double f = par[0]*(term4/(term5*term5 + term6*term6))*TMath::Power((par[1]/x[0]),par[3]) + par[4] + par[5]*x[0];

    return f;
}
double RossStodPol2(double *x,double *par){
    // par[0] is amplitude
    // par[1] is mass
    // par[2] is width
    // par[3] is n
    double term1 = x[0]*x[0] - 4*pimass*pimass;
    double term2 = par[1]*par[1]-4*pimass*pimass;
    if(term2==0) term2=-9999; // want a different error checker
    double term3 = term1/term2;
    double gamma = par[2]*(par[1]/x[0])*TMath::Power(term3,1.5);

    double term4 = x[0]*par[1]*gamma;
    double term5 = par[1]*par[1] - x[0]*x[0];
    double term6 = par[1]*gamma;

    double f = par[0]*(term4/(term5*term5 + term6*term6))*TMath::Power((par[1]/x[0]),par[3]) + par[4] + par[5]*x[0] +par[6]*x[0]*x[0];

    return f;
}
double RossStodExp(double *x,double *par){
    // par[0] is amplitude
    // par[1] is mass
    // par[2] is width
    // par[3] is n
    double term1 = x[0]*x[0] - 4*pimass*pimass;
    double term2 = par[1]*par[1]-4*pimass*pimass;
    if(term2==0) term2=-9999; // want a different error checker
    double term3 = term1/term2;
    double gamma = par[2]*(par[1]/x[0])*TMath::Power(term3,1.5);

    double term4 = x[0]*par[1]*gamma;
    double term5 = par[1]*par[1] - x[0]*x[0];
    double term6 = par[1]*gamma;

    double f = par[0]*(term4/(term5*term5 + term6*term6))*TMath::Power((par[1]/x[0]),par[3]) + par[4] + par[5]*TMath::Exp(-par[6]*x[0]);

    return f;
}
double fPol1(double *x,double *par){
    double f = par[0] + par[1]*x[0];
    return f;
}
double fPol2(double *x,double *par){
    double f = par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
    return f;
}
double fExp(double *x,double *par){
    double f = par[0] + par[1]*TMath::Exp(-par[2]*x[0]);
    return f;
}
// all fits have an offset parameter added for plotting purposes
double fGauss(double *x,double *par){
    double arg = 0;
    if(par[2]!=0) arg = (x[0] - par[1])/par[2];
    double fitval = (par[0]/(TMath::Sqrt(2*TMath::Pi())*par[2]))*TMath::Exp(-0.5*arg*arg) + par[3];
    return fitval;
}
double fBreitW(double *x,double *par){
    // par[0] = Amplitude
    // par[1] = mass
    // par[2] = width

    double gam = TMath::Sqrt(par[1]*par[1]*(par[1]*par[1]+par[2]*par[2]));
    double k=(2*roo2*par[1]*par[2]*gam)/(PI*TMath::Sqrt(par[1]*par[1]+gam));
    double sim1=x[0]*x[0] - par[1]*par[1];
    double sim2=par[1]*par[2];


    double f = (par[0]*k)/(sim1*sim1 + sim2*sim2) + par[3];

    return f;
}
double fBreitSod(double *x,double *par){
    // par[0] is amplitude of breit term
    // par[1] is f0 mass
    // par[2] is f0 width
    // par[3] is amplitude of sod term
    double term1 = x[0]*par[1]*par[2];
    double term2 = par[1]*par[1]-x[0]*x[0];
    double term3 = par[1]*par[2];
    // BW + Sod + pol
    double f = (par[0]*term1)/(term2*term2 + term3*term3) + (par[3]*term2)/(term2*term2 + term3*term3) + par[4];
    return f;
}
double fmodSod(double *x,double *par){
    // par[0] is first term amplitude (amp of "A")
    // par[1] is mass
    // par[2] is width
    // par[3] is interference term (amp of "B") : check reference
    // I'm assuming that A and B are real, but there is no reason they need to be
    // so I should look into these being complex
    double a = par[0]*TMath::Sqrt(x[0]*par[1]*par[2]);
    double c = x[0]*x[0] - par[1]*par[1];
    double d = par[1]*par[2];

    double f = (a*a + 2*a*c*par[3])/(c*c + d*d) + par[3]*par[3] + par[4];

    return f;
}
double fRossStod(double *x,double *par){
    // par[0] is amplitude
    // par[1] is mass
    // par[2] is width
    // par[3] is n
    double term1 = x[0]*x[0] - 4*pimass*pimass;
    double term2 = par[1]*par[1]-4*pimass*pimass;
    if(term2==0) term2=-9999; // want a different error checker
    double term3 = term1/term2;
    double gamma = par[2]*(par[1]/x[0])*TMath::Power(term3,1.5);

    double term4 = x[0]*par[1]*gamma;
    double term5 = par[1]*par[1] - x[0]*x[0];
    double term6 = par[1]*gamma;

    double f = par[0]*(term4/(term5*term5 + term6*term6))*TMath::Power((par[1]/x[0]),par[3]) +par[4];

    return f;

}

double curveBreitW(double *x,double *par)
{
    // par[0] = Amplitude
    // par[1] = mass
    // par[2] = width

    double gam = TMath::Sqrt(par[1]*par[1]*(par[1]*par[1]+par[2]*par[2]));
    double k=(2*roo2*par[1]*par[2]*gam)/(PI*TMath::Sqrt(par[1]*par[1]+gam));
    double sim1=x[0]*x[0] - par[1]*par[1];
    double sim2=par[1]*par[2];


    double f = (par[0]*k)/(sim1*sim1 + sim2*sim2);

    return f;

}

double curveBreitSod(double *x,double *par)
{
    // par[0] is amplitude of breit term
    // par[1] is f0 mass
    // par[2] is f0 width
    // par[3] is amplitude of sod term
    double term1 = x[0]*par[1]*par[2];
    double term2 = par[1]*par[1]-x[0]*x[0];
    double term3 = par[1]*par[2];
    // BW + Sod + pol
    double f = (par[0]*term1)/(term2*term2 + term3*term3) + (par[3]*term2)/(term2*term2 + term3*term3);
    return f;

}

double curveModSod(double *x,double *par)
{
    // par[0] is first term amplitude (amp of "A")
    // par[1] is mass
    // par[2] is width
    // par[3] is interference term (amp of "B") : check reference
    // I'm assuming that A and B are real, but there is no reason they need to be
    // so I should look into these being complex
    double a = par[0]*TMath::Sqrt(x[0]*par[1]*par[2]);
    double c = x[0]*x[0] - par[1]*par[1];
    double d = par[1]*par[2];

    double f = (a*a + 2*a*c*par[3])/(c*c + d*d) + par[3]*par[3];

    return f;

}

double curveRossStod(double *x,double *par)
{
    // par[0] is amplitude
    // par[1] is mass
    // par[2] is width
    // par[3] is n
    double term1 = x[0]*x[0] - 4*pimass*pimass;
    double term2 = par[1]*par[1]-4*pimass*pimass;
    if(term2==0) term2=-9999; // want a different error checker
    double term3 = term1/term2;
    double gamma = par[2]*(par[1]/x[0])*TMath::Power(term3,1.5);

    double term4 = x[0]*par[1]*gamma;
    double term5 = par[1]*par[1] - x[0]*x[0];
    double term6 = par[1]*gamma;

    double f = par[0]*(term4/(term5*term5 + term6*term6))*TMath::Power((par[1]/x[0]),par[3]);

    return f;

}

double v2cosfit(double *x,double *par)
{
    double arg = par[0]*(1+2*par[1]*TMath::Cos(2*x[0]));
    return arg;
}
double expDecay(double *x,double *par)
{
    double arg = par[0]*TMath::Exp(-par[1]*x[0]);
    return arg;
}

// the fit pointers only differ by the name
TF1 *bkgFit(vector<double> &vals,double (*fitType) (double *,double *),int someNum)
{
    char name[50];
    sprintf(name,"bkg%i",someNum);

    TF1 *bkg = new TF1(name,fitType,0,3,vals.size());

    for(int i=0;i<vals.size();i++)
    {
        bkg->SetParameter(i,vals[i]); 
    }

    bkg->SetRange(0.9,1.1);

    return bkg;
}

TF1 *funcFit(vector<double> &vals,double (*fitType) (double *, double *),int someNum)
{
    char name[50];
    sprintf(name,"func%i",someNum);

    TF1 *func = new TF1(name,fitType,0,3,vals.size());

    for(int i=0;i<vals.size();i++)
    {
        func->SetParameter(i,vals[i]);
    }

    func->SetRange(0.9,1.1);

    return func;
}


#endif
