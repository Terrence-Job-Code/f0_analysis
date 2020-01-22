#include "rootheader.h"
#include "fitvalues.h"
#include "auxFun.cpp"

#ifndef PRIMARYFIT_H_
#define PRIMARYFIT_H_
class primaryFit
{
    private:
        // consider putting m_ in front of all member data, and maybe something for all functions

        //const int numCents=9;
        TFile *file;
        int cent;
        TH3D *m_ULhists[numCents];
        TH3D *m_LShists[numCents];
     
        // for reference
        // centrality number and percentage
        // 9: 0-5% these are the most central collisions
        // 8: 5-10%
        // 7: 10-20%
        // 6: 20-30%
        // 5: 30-40%
        // 4: 40-50%
        // 3: 50-60%
        // 2: 60-70%
        // 1: 70-80%

        vector < vector<double> > m_v2;
        vector < vector<double> > m_v2Err;
        // will only have two elements, it will be the same for all of them
        vector<double> m_weightedv2pT;
        vector<double> m_weightedv2pTerr;

        // Could include struct as a data type, would reduce inputs a lot 
        fitContainer m_infoBox;
        

    public:
        const double f0mass=0.9894; // from pdg table
        const double f0width_KK=15.3e-3; // from KK decay, I think: check pdg table
        //const double pimass=0.13957;
        //const double roo2=TMath::Sqrt(2); 
        //const double PI=TMath::Pi();
        
        primaryFit(TString fileName,int cent);
        ~primaryFit();

        // fills the class data member fit container
        void fillConstr(fitContainer intermediateBox); // consider a better name 

        TH1D *funFit(int *bins, int count1, int count2, int count3, int countPass);
        void basicCurve(int runNum); 
        void legend(); // creates the legend  

        void Analysis(struct fitContainer infoBox);
        void multHistPlot(vector < vector<TH1D*> > &hists,vector < vector<TF1*> > &TheFit,vector < vector<TF1*> > &TheBck,TString type,fitContainer infoBox);

        void YieldPlot(vector < vector<double> > &yield,vector < vector<double> > &yieldErr,vector < vector<double> > &xAxis, 
                vector < vector<double> > &xAxisErr,fitContainer infoBox);
        
        void weightv2PT(); 
        
        void AreaCalculator( double (*fitType) (double *,double *), vector < vector<double> > &storArea, 
            vector < vector<double> > &storAreaErr, vector < vector < vector<double> > > &functionParams);

        double NumericalIntegCalc(TF1 *func,double lowerLim,double upperLim);
        double YieldCalc(TH1D *histo,FitFuncPtr bkgFunc,vector<double> &bkgvals);
        double YieldErrorCalc(TH1D *histo,FitFuncPtr bkgFunc,vector<double> &bkgvals);

        // these will likely be changed to return histograms for the sake of easy comparison
        void f0splits();
        void singleParticlepT();

        // hard coded functions mostly unrelated to the rest of the class
        void compareRuns_01();
        void compareRuns_02();    
        void SinglePartDistrs();
        void f0pTtests();
        void correctionFit();
};
#endif
