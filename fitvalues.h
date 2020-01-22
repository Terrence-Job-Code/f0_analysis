#include "rootheader.h"
#include "fits.h"

using namespace std;

#ifndef FITVALUES_H_
#define FITVALUES_H_

TString dir="mainPlots/";

const int numCents=9;
const int numPtDisBins=6;
const int numPhiDisBins=5;
const int numBinBounds=4;
const int numStringsFinFun=7;

const int numPol1Par=2;
const int numPol2Par=3;
const int numExpPar=3;
const int numGauPar=4;
const int numBreitWPar=4;
const int numBreitSodPar=5;
const int numModSodPar=5;
const int numRossStodPar=5;

typedef double (*FitFuncPtr) (double *, double *);

//////////  ONLY GAUSSIAN HAS AREA PARAMETER, REST ARE AMPLITUDES, HAVE TO USE A FUNCTION TO CALCULATE AREA FOR YIELDS

// pt distribution over all phi
int ptDisBin[numPtDisBins][numBinBounds]={
    {1,2,1,20},
    {3,4,1,20},
    {5,6,1,20},
    {7,8,1,20},
    {9,12,1,20},
    {13,20,1,20}
};

const int numPhiDistrs=2; // low and hight pt range

// phi distributions 0<pt<2 GeV , 2<pt<5 GeV
int phiDisBin[numPhiDistrs][numPhiDisBins][numBinBounds]={
    {
        {1,8,1,4},
        {1,8,5,8},
        {1,8,9,12},
        {1,8,13,16},
        {1,8,17,20}
    },
    {
        {9,20,1,4},
        {9,20,5,8},
        {9,20,9,12},
        {9,20,13,16},
        {9,20,17,20}
    }
};

///// BEGINNING OF RUN11 FITS ///////

// there is no reason for the bounds to be specific to the gaussian fit
// will leave these alone for the moment, might not need them
vector < vector<double> > YaxisGauPtPol1={
    {5e4,5e6},
    {5e4,5e6},
    {5e2,5e6},
    {-5e3,3e6},
    {5e4,5e5},
    {-1.5e3,7e3}
};
vector < vector<double> > YaxisGauPtPol2={
    {5e4,5e6},
    {5e4,5e6},
    {5e2,5e6},
    {-5e3,3e6},
    {5e4,5e5},
    {-1.5e3,7e3}
};
vector < vector<double> > YaxisGauPtExp={
    {5e4,5e6},
    {5e4,5e6},
    {5e2,5e6},
    {-5e3,3e6},
    {5e4,5e5},
    {-1.5e3,7e3}
};
vector < vector<double> > YaxisGauPhiPol1lpt={
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6}
};
vector < vector<double> > YaxisGauPhiPol2lpt={
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6}
};
vector < vector<double> > YaxisGauPhiExplpt={
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6},
    {5e4,3.5e6}
};
vector < vector<double> > YaxisGauPhiPol1hpt={
    {5e3,2e5},
    {2e4,1.6e5},
    {2e4,1.2e5},
    {1e4,1e5},
    {1e4,1e5}
};
vector < vector<double> > YaxisGauPhiPol2hpt={
    {5e3,2e5},
    {2e4,1.6e5},
    {2e4,1.2e5},
    {1e4,1e5},
    {1e4,1e5}
};
vector < vector<double> > YaxisGauPhiExphpt={
    {5e3,2e5},
    {2e4,1.6e5},
    {2e4,1.2e5},
    {1e4,1e5},
    {1e4,1e5}
};

// Gaussian fit information
const int numGauPol1Par=5;
const int numGauPol2Par=6;
const int numGauExpPar=6;

// par[0] is area
// par[1] is mass
// par[2] is width
// for exponential
// par[3] is offset
// par[4] is amplitude
// par[5] is decay

vector < vector<FitFuncPtr> > GaussFuncs={
    {fitGausPol1,fPol1,fGauss},
    {fitGausPol2,fPol2,fGauss},
    {fitGausExp,fExp,fGauss}
};

vector < vector<int> > GaussNumParams={
    {numGauPol1Par,numPol1Par,numGauPar},
    {numGauPol2Par,numPol2Par,numGauPar},
    {numGauExpPar,numExpPar,numGauPar}
};

vector < vector<double> > run11parBoundsGauMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.97,1.0,5e-3,20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {3e3,4e3,0.98,0.99,5e-3, 20e-3,-1e8,7e7,-1e8,-7.1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsGauMassPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.97,1.00,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsGauMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsGauMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsGauMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsGauMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsGauMassPhiPol1hpt={
    {460,1e5,0.95,1.05,5e-3,20e-3,-1e6,1e6,-1e6,1e6},
    {0,1e5,0.95,1.05,5e-3,50e-3,-1e6,1e6,-1e6,1e6},
    {0,1e5,0.95,1.05,5e-3,50e-3,-1e6,1e6,-1e6,1e6},
    {0,1e5,0.95,1.05,5e-3,20e-3,-1e6,1e6,-1e6,1e6},
    {4,1e5,0.95,1.05,5e-3,50e-3,-1e6,1e6,-1e6,1e6}
};
vector < vector<double> > run11parBoundsGauMassPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsGauMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11sGauMassPtPol1={
    {"aGauPt1Pol1","bGauPt1Pol1","cGauPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt2Pol1","bGauPt2Pol1","cGauPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt3Pol1","bGauPt3Pol1","cGauPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt4Pol1","bGauPt4Pol1","cGauPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt5Pol1","bGauPt5Pol1","cGauPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt6Pol1","bGauPt6Pol1","cGauPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run11sGauMassPtPol2={
    {"aGauPt1Pol2","bGauPt1Pol2","cGauPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt2Pol2","bGauPt2Pol2","cGauPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt3Pol2","bGauPt3Pol2","cGauPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt4Pol2","bGauPt4Pol2","cGauPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt5Pol2","bGauPt5Pol2","cGauPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt6Pol2","bGauPt6Pol2","cGauPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11sGauMassPtExp={
    {"aGauPt1exp","bGauPt1exp","cGauPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt2exp","bGauPt2exp","cGauPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt3exp","bGauPt3exp","cGauPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt4exp","bGauPt4exp","cGauPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt5exp","bGauPt5exp","cGauPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt6exp","bGauPt6exp","cGauPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11sGauMassPhiPol1lpt={
    {"aGauPhi1Pol1lpt","bGauPhi1Pol1lpt","cGauPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1lpt","bGauPhi2Pol1lpt","cGauPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1lpt","bGauPhi3Pol1lpt","cGauPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1lpt","bGauPhi4Pol1lpt","cGauPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1lpt","bGauPhi5Pol1lpt","cGauPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run11sGauMassPhiPol2lpt={
    {"aGauPhi1Pol2lpt","bGauPhi1Pol2lpt","cGauPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2lpt","bGauPhi2Pol2lpt","cGauPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2lpt","bGauPhi3Pol2lpt","cGauPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2lpt","bGauPhi4Pol2lpt","cGauPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2lpt","bGauPhi5Pol2lpt","cGauPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sGauMassPhiExplpt={
    {"aGauPhi1Explpt","bGauPhi1Explpt","cGauPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Explpt","bGauPhi2Explpt","cGauPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Explpt","bGauPhi3Explpt","cGauPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Explpt","bGauPhi4Explpt","cGauPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Explpt","bGauPhi5Explpt","cGauPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11sGauMassPhiPol1hpt={
    {"aGauPhi1Pol1hpt","bGauPhi1Pol1hpt","cGauPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1hpt","bGauPhi2Pol1hpt","cGauPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1hpt","bGauPhi3Pol1hpt","cGauPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1hpt","bGauPhi4Pol1hpt","cGauPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1hpt","bGauPhi5Pol1hpt","cGauPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run11sGauMassPhiPol2hpt={
    {"aGauPhi1Pol2hpt","bGauPhi1Pol2hpt","cGauPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2hpt","bGauPhi2Pol2hpt","cGauPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2hpt","bGauPhi3Pol2hpt","cGauPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2hpt","bGauPhi4Pol2hpt","cGauPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2hpt","bGauPhi5Pol2hpt","cGauPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sGauMassPhiExphpt={
    {"aGauPhi1Exphpt","bGauPhi1Exphpt","cGauPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Exphpt","bGauPhi2Exphpt","cGauPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Exphpt","bGauPhi3Exphpt","cGauPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Exphpt","bGauPhi4Exphpt","cGauPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Exphpt","bGauPhi5Exphpt","cGauPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};

//vector < vector<TString> > run11v2GaussNames={
//    {"dN/dPhi vs Phi (lin) Bkg 0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
//    {"dN/dPhi vs Phi (quad) Bkg 0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
//    {"dN/dPhi vs Phi (exp) Bkg 0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
//    {"dN/dPhi vs Phi (lin) Bkg 2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
//    {"dN/dPhi vs Phi (quad) Bkg 2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
//    {"dN/dPhi vs Phi (exp) Bkg 2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"}
//};

// Breit Wigner fit information
const int numBreitWpol1Par=5;
const int numBreitWpol2Par=6;
const int numBreitWexpPar=6;

vector < vector<FitFuncPtr> > BreitWigFuncs={
    {BreitWpol1,fPol1,fBreitW},
    {BreitWpol2,fPol2,fBreitW},
    {BreitWexp,fExp,fBreitW}
};

vector < vector<int> > BreitWigNumParams={
    {numBreitWpol1Par,numPol1Par,numBreitWPar},
    {numBreitWpol2Par,numPol2Par,numBreitWPar},
    {numBreitWexpPar,numExpPar,numBreitWPar}
};

vector < vector<double> > run11parBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e5,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2.5e3,1e4,0.97,1.00,15e-3,40e-3,1e6,1e8,-1e8,-5e7,1e6,5e7},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,60e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11sBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run11sBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11sBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11sBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run11sBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11sBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run11sBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};

// Breit Soding fit information
const int numBreitSodPol1Par=6;
const int numBreitSodPol2Par=7;
const int numBreitSodExpPar=7;

vector < vector<FitFuncPtr> > BreitSodFuncs={
    {BreitSodPol1,fPol1,fBreitSod},
    {BreitSodPol2,fPol2,fBreitSod},
    {BreitSodExp,fExp,fBreitSod}
};

vector < vector<int> > BreitSodNumParams={
    {numBreitSodPol1Par,numPol1Par,numBreitSodPar},
    {numBreitSodPol2Par,numPol2Par,numBreitSodPar},
    {numBreitSodExpPar,numExpPar,numBreitSodPar}
};

vector < vector<double> > run11parBoundsBreitSodMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,-8e6},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,-4e6},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitSodMassPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e5,0.95,1.05,5e-3,100e-3,-1e3,1e3,1e5,1e8,-1e8,-1e7,5e5,5e7},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e7,-5e4}
};
vector < vector<double> > run11parBoundsBreitSodMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,1e6,1e8,1,3},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,1e6,1e8,1,3},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,1e6,1e8,1,3},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,1e6,1e8,1,3},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,1e6,1e8,1,3},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,0,1e8,1,3}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiExplpt={
    {-1e3,1e4,0.95,1.05,5e-3,100e-3,-1e5,1e5,-1e5,1e8,-1e8,1e8,1,10},
    {-1e3,2e4,0.95,1.05,5e-3,100e-3,-1e5,1e5,-1e5,1e8,-1e8,1e8,1,10},
    {-1e3,2e4,0.95,1.05,5e-3,100e-3,-1e5,1e5,-1e5,1e8,-1e8,1e8,1,10},
    {-1e3,2e4,0.95,1.05,5e-3,100e-3,-1e5,1e5,-1e5,1e8,-1e8,1e8,1,10},
    {-1e3,2e4,0.95,1.05,5e-3,100e-3,-1e5,1e5,-1e5,1e8,-1e8,1e8,1,10}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e2,1e4,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e2,1e4,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsBreitSodMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e6,1e8,-1e8,1e8,1,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,-1e6,1e8,-1e8,1e8,1,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,-1e6,1e8,-1e8,1e8,1,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,-1e6,1e8,-1e8,1e8,1,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e6,1e8,-1e8,1e8,1,10}
};

vector < vector<TString> > run11sBreitSodMassPtPol1={
    {"aBreitSodPt1Pol1","bBreitSodPt1Pol1","cBreitSodPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPt2Pol1","bBreitSodPt2Pol1","cBreitSodPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPt3Pol1","bBreitSodPt3Pol1","cBreitSodPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPt4Pol1","bBreitSodPt4Pol1","cBreitSodPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPt5Pol1","bBreitSodPt5Pol1","cBreitSodPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPt6Pol1","bBreitSodPt6Pol1","cBreitSodPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"}
};

vector < vector<TString> > run11sBreitSodMassPtPol2={
    {"aBreitSodPt1Pol2","bBreitSodPt1Pol2","cBreitSodPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"},
    {"aBreitSodPt2Pol2","bBreitSodPt2Pol2","cBreitSodPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"},
    {"aBreitSodPt3Pol2","bBreitSodPt3Pol2","cBreitSodPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"},
    {"aBreitSodPt4Pol2","bBreitSodPt4Pol2","cBreitSodPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"},
    {"aBreitSodPt5Pol2","bBreitSodPt5Pol2","cBreitSodPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"},
    {"aBreitSodPt6Pol2","bBreitSodPt6Pol2","cBreitSodPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11sBreitSodMassPtExp={
    {"aBreitSodPt1exp","bBreitSodPt1exp","cBreitSodPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPt2exp","bBreitSodPt2exp","cBreitSodPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPt3exp","bBreitSodPt3exp","cBreitSodPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPt4exp","bBreitSodPt4exp","cBreitSodPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPt5exp","bBreitSodPt5exp","cBreitSodPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPt6exp","bBreitSodPt6exp","cBreitSodPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11sBreitSodMassPhiPol1lpt={
    {"aBreitSodPhi1Pol1lpt","bBreitSodPhi1Pol1lpt","cBreitSodPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi2Pol1lpt","bBreitSodPhi2Pol1lpt","cBreitSodPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi3Pol1lpt","bBreitSodPhi3Pol1lpt","cBreitSodPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi4Pol1lpt","bBreitSodPhi4Pol1lpt","cBreitSodPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi5Pol1lpt","bBreitSodPhi5Pol1lpt","cBreitSodPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"}
};
vector < vector<TString> > run11sBreitSodMassPhiPol2lpt={
    {"aBreitSodPhi1Pol2lpt","bBreitSodPhi1Pol2lpt","cBreitSodPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi2Pol2lpt","bBreitSodPhi2Pol2lpt","cBreitSodPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi3Pol2lpt","bBreitSodPhi3Pol2lpt","cBreitSodPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi4Pol2lpt","bBreitSodPhi4Pol2lpt","cBreitSodPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi5Pol2lpt","bBreitSodPhi5Pol2lpt","cBreitSodPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sBreitSodMassPhiExplpt={
    {"aBreitSodPhi1Explpt","bBreitSodPhi1Explpt","cBreitSodPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi2Explpt","bBreitSodPhi2Explpt","cBreitSodPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi3Explpt","bBreitSodPhi3Explpt","cBreitSodPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi4Explpt","bBreitSodPhi4Explpt","cBreitSodPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi5Explpt","bBreitSodPhi5Explpt","cBreitSodPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11sBreitSodMassPhiPol1hpt={
    {"aBreitSodPhi1Pol1hpt","bBreitSodPhi1Pol1hpt","cBreitSodPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi2Pol1hpt","bBreitSodPhi2Pol1hpt","cBreitSodPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi3Pol1hpt","bBreitSodPhi3Pol1hpt","cBreitSodPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi4Pol1hpt","bBreitSodPhi4Pol1hpt","cBreitSodPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"},
    {"aBreitSodPhi5Pol1hpt","bBreitSodPhi5Pol1hpt","cBreitSodPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1"}

};
vector < vector<TString> > run11sBreitSodMassPhiPol2hpt={
    {"aBreitSodPhi1Pol2hpt","bBreitSodPhi1Pol2hpt","cBreitSodPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi2Pol2hpt","bBreitSodPhi2Pol2hpt","cBreitSodPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi3Pol2hpt","bBreitSodPhi3Pol2hpt","cBreitSodPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi4Pol2hpt","bBreitSodPhi4Pol2hpt","cBreitSodPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"},
    {"aBreitSodPhi5Pol2hpt","bBreitSodPhi5Pol2hpt","cBreitSodPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sBreitSodMassPhiExphpt={
    {"aBreitSodPhi1Exphpt","bBreitSodPhi1Exphpt","cBreitSodPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi2Exphpt","bBreitSodPhi2Exphpt","cBreitSodPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi3Exphpt","bBreitSodPhi3Exphpt","cBreitSodPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi4Exphpt","bBreitSodPhi4Exphpt","cBreitSodPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"},
    {"aBreitSodPhi5Exphpt","bBreitSodPhi5Exphpt","cBreitSodPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Sodinger","mass GeV","counts",dir+"saves.pdf","Breit Amp","mass","width","Sod Amp","offset","exp amp","decay"}
};


// modified Soding fit information
const int numModSodPol1Par=6;
const int numModSodPol2Par=7;
const int numModSodExpPar=7;

vector < vector<FitFuncPtr> > ModSodFuncs={
    {modSodPol1,fPol1,fmodSod},
    {modSodPol2,fPol2,fmodSod},
    {modSodExp,fExp,fmodSod}
};

vector < vector<int> > ModSodNumParams={
    {numModSodPol1Par,numPol1Par,numModSodPar},
    {numModSodPol2Par,numPol2Par,numModSodPar},
    {numModSodExpPar,numExpPar,numModSodPar}
};

vector < vector<double> > run11parBoundsModSodMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e5,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsModSodMassPtPol2={
    {0,1e3,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {10,1e3,0.95,1.05,5e-3,100e-3,-1e3,1e3,1e7,5e7,-1e8,-1e7,1e7,5e7},
    {0,1e3,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {10,1e3,0.97,1.0,5e-3,100e-3,-1e3,1e3,4e7,6e7,-1.1e8,-8e7,1e7,5e7},
    {0,1e3,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {10,1e3,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsModSodMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsModSodMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,5e6,7e6,-6.5e6,-1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,5e6,7e6,-6.5e6,-1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,5e6,7e6,-6.5e6,-1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,5e6,7e6,-6.5e6,-1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,5e6,7e6,-6.5e6,-1e5}
};
vector < vector<double> > run11parBoundsModSodMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,1e6,1e8,-1e8,-1e7,1e7,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsModSodMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsModSodMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,1e4,1e6,-1.4e5,-1e4},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsModSodMassPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,10,-4.45e5,-4e5,1e6,2e6,-6.7e5,-5e5},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e3,1e3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsModSodMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11sModSodMassPtPol1={
    {"aModSodPt1Pol1","bModSodPt1Pol1","cModSodPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPt2Pol1","bModSodPt2Pol1","cModSodPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPt3Pol1","bModSodPt3Pol1","cModSodPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPt4Pol1","bModSodPt4Pol1","cModSodPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPt5Pol1","bModSodPt5Pol1","cModSodPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPt6Pol1","bModSodPt6Pol1","cModSodPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"}
};

vector < vector<TString> > run11sModSodMassPtPol2={
    {"aModSodPt1Pol2","bModSodPt1Pol2","cModSodPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"},
    {"aModSodPt2Pol2","bModSodPt2Pol2","cModSodPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"},
    {"aModSodPt3Pol2","bModSodPt3Pol2","cModSodPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"},
    {"aModSodPt4Pol2","bModSodPt4Pol2","cModSodPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"},
    {"aModSodPt5Pol2","bModSodPt5Pol2","cModSodPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"},
    {"aModSodPt6Pol2","bModSodPt6Pol2","cModSodPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11sModSodMassPtExp={
    {"aModSodPt1exp","bModSodPt1exp","cModSodPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPt2exp","bModSodPt2exp","cModSodPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPt3exp","bModSodPt3exp","cModSodPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPt4exp","bModSodPt4exp","cModSodPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPt5exp","bModSodPt5exp","cModSodPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPt6exp","bModSodPt6exp","cModSodPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11sModSodMassPhiPol1lpt={
    {"aModSodPhi1Pol1lpt","bModSodPhi1Pol1lpt","cModSodPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi2Pol1lpt","bModSodPhi2Pol1lpt","cModSodPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi3Pol1lpt","bModSodPhi3Pol1lpt","cModSodPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi4Pol1lpt","bModSodPhi4Pol1lpt","cModSodPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi5Pol1lpt","bModSodPhi5Pol1lpt","cModSodPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"}
};
vector < vector<TString> > run11sModSodMassPhiPol2lpt={
    {"aModSodPhi1Pol2lpt","bModSodPhi1Pol2lpt","cModSodPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi2Pol2lpt","bModSodPhi2Pol2lpt","cModSodPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi3Pol2lpt","bModSodPhi3Pol2lpt","cModSodPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi4Pol2lpt","bModSodPhi4Pol2lpt","cModSodPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi5Pol2lpt","bModSodPhi5Pol2lpt","cModSodPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sModSodMassPhiExplpt={
    {"aModSodPhi1Explpt","bModSodPhi1Explpt","cModSodPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi2Explpt","bModSodPhi2Explpt","cModSodPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi3Explpt","bModSodPhi3Explpt","cModSodPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi4Explpt","bModSodPhi4Explpt","cModSodPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi5Explpt","bModSodPhi5Explpt","cModSodPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11sModSodMassPhiPol1hpt={
    {"aModSodPhi1Pol1hpt","bModSodPhi1Pol1hpt","cModSodPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi2Pol1hpt","bModSodPhi2Pol1hpt","cModSodPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi3Pol1hpt","bModSodPhi3Pol1hpt","cModSodPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi4Pol1hpt","bModSodPhi4Pol1hpt","cModSodPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"},
    {"aModSodPhi5Pol1hpt","bModSodPhi5Pol1hpt","cModSodPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1"}

};
vector < vector<TString> > run11sModSodMassPhiPol2hpt={
    {"aModSodPhi1Pol2hpt","bModSodPhi1Pol2hpt","cModSodPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi2Pol2hpt","bModSodPhi2Pol2hpt","cModSodPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi3Pol2hpt","bModSodPhi3Pol2hpt","cModSodPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi4Pol2hpt","bModSodPhi4Pol2hpt","cModSodPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"},
    {"aModSodPhi5Pol2hpt","bModSodPhi5Pol2hpt","cModSodPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sModSodMassPhiExphpt={
    {"aModSodPhi1Exphpt","bModSodPhi1Exphpt","cModSodPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi2Exphpt","bModSodPhi2Exphpt","cModSodPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi3Exphpt","bModSodPhi3Exphpt","cModSodPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi4Exphpt","bModSodPhi4Exphpt","cModSodPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"},
    {"aModSodPhi5Exphpt","bModSodPhi5Exphpt","cModSodPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Modified Sodinger","mass GeV","counts",dir+"saves.pdf","A amp","mass","width","B amp","offset","exp amp","decay"}
};



// Ross Stodolsky fit information
const int numRossStodPol1Par=6;
const int numRossStodPol2Par=7;
const int numRossStodExpPar=7;

vector < vector<FitFuncPtr> > RossStodFuncs={
    {modSodPol1,fPol1,fmodSod},
    {modSodPol2,fPol2,fmodSod},
    {modSodExp,fExp,fmodSod}
};

vector < vector<int> > RossStodNumParams={
    {numRossStodPol1Par,numPol1Par,numRossStodPar},
    {numRossStodPol2Par,numPol2Par,numRossStodPar},
    {numRossStodExpPar,numExpPar,numRossStodPar}
};

vector < vector<double> > run11parBoundsRossStodMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e5,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPtPol2={
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,4.5e7,5.2e7,-1.2e8,-7e7,3.5e7,1e8},
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,5e-3,100e-3,1,50,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsRossStodMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};
vector < vector<double> > run11parBoundsRossStodMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPhiPol2hpt={
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,1e2,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,1e2,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,1e2,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,30,-4e5,-3.5e5,1e6,5e6,-7e5,-1e5},
    {0,1e8,0.97,1.0,25e-3,100e-3,-1e2,1e2,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11parBoundsRossStodMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11sRossStodMassPtPol1={
    {"aRossStodPt1Pol1","bRossStodPt1Pol1","cRossStodPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPt2Pol1","bRossStodPt2Pol1","cRossStodPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPt3Pol1","bRossStodPt3Pol1","cRossStodPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPt4Pol1","bRossStodPt4Pol1","cRossStodPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPt5Pol1","bRossStodPt5Pol1","cRossStodPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPt6Pol1","bRossStodPt6Pol1","cRossStodPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"}
};

vector < vector<TString> > run11sRossStodMassPtPol2={
    {"aRossStodPt1Pol2","bRossStodPt1Pol2","cRossStodPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"},
    {"aRossStodPt2Pol2","bRossStodPt2Pol2","cRossStodPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"},
    {"aRossStodPt3Pol2","bRossStodPt3Pol2","cRossStodPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"},
    {"aRossStodPt4Pol2","bRossStodPt4Pol2","cRossStodPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"},
    {"aRossStodPt5Pol2","bRossStodPt5Pol2","cRossStodPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"},
    {"aRossStodPt6Pol2","bRossStodPt6Pol2","cRossStodPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11sRossStodMassPtExp={
    {"aRossStodPt1exp","bRossStodPt1exp","cRossStodPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPt2exp","bRossStodPt2exp","cRossStodPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPt3exp","bRossStodPt3exp","cRossStodPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPt4exp","bRossStodPt4exp","cRossStodPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPt5exp","bRossStodPt5exp","cRossStodPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPt6exp","bRossStodPt6exp","cRossStodPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11sRossStodMassPhiPol1lpt={
    {"aRossStodPhi1Pol1lpt","bRossStodPhi1Pol1lpt","cRossStodPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi2Pol1lpt","bRossStodPhi2Pol1lpt","cRossStodPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi3Pol1lpt","bRossStodPhi3Pol1lpt","cRossStodPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi4Pol1lpt","bRossStodPhi4Pol1lpt","cRossStodPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi5Pol1lpt","bRossStodPhi5Pol1lpt","cRossStodPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"}
};
vector < vector<TString> > run11sRossStodMassPhiPol2lpt={
    {"aRossStodPhi1Pol2lpt","bRossStodPhi1Pol2lpt","cRossStodPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi2Pol2lpt","bRossStodPhi2Pol2lpt","cRossStodPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi3Pol2lpt","bRossStodPhi3Pol2lpt","cRossStodPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi4Pol2lpt","bRossStodPhi4Pol2lpt","cRossStodPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi5Pol2lpt","bRossStodPhi5Pol2lpt","cRossStodPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sRossStodMassPhiExplpt={
    {"aRossStodPhi1Explpt","bRossStodPhi1Explpt","cRossStodPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi2Explpt","bRossStodPhi2Explpt","cRossStodPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi3Explpt","bRossStodPhi3Explpt","cRossStodPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi4Explpt","bRossStodPhi4Explpt","cRossStodPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi5Explpt","bRossStodPhi5Explpt","cRossStodPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11sRossStodMassPhiPol1hpt={
    {"aRossStodPhi1Pol1hpt","bRossStodPhi1Pol1hpt","cRossStodPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi2Pol1hpt","bRossStodPhi2Pol1hpt","cRossStodPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi3Pol1hpt","bRossStodPhi3Pol1hpt","cRossStodPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi4Pol1hpt","bRossStodPhi4Pol1hpt","cRossStodPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"},
    {"aRossStodPhi5Pol1hpt","bRossStodPhi5Pol1hpt","cRossStodPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1"}

};
vector < vector<TString> > run11sRossStodMassPhiPol2hpt={
    {"aRossStodPhi1Pol2hpt","bRossStodPhi1Pol2hpt","cRossStodPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi2Pol2hpt","bRossStodPhi2Pol2hpt","cRossStodPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi3Pol2hpt","bRossStodPhi3Pol2hpt","cRossStodPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi4Pol2hpt","bRossStodPhi4Pol2hpt","cRossStodPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"},
    {"aRossStodPhi5Pol2hpt","bRossStodPhi5Pol2hpt","cRossStodPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","poly 1","poly2"}
};
vector < vector<TString> > run11sRossStodMassPhiExphpt={
    {"aRossStodPhi1Exphpt","bRossStodPhi1Exphpt","cRossStodPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi2Exphpt","bRossStodPhi2Exphpt","cRossStodPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi3Exphpt","bRossStodPhi3Exphpt","cRossStodPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi4Exphpt","bRossStodPhi4Exphpt","cRossStodPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"},
    {"aRossStodPhi5Exphpt","bRossStodPhi5Exphpt","cRossStodPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Ross-Stodolsky","mass GeV","counts",dir+"saves.pdf","Ross Amp","mass","width","n","offset","exp amp","decay"}
};





////////////// BEGINNING OF RUN 16 FITS ///////////////////

/// RUN 16 GAUSSIAN FIT

vector < vector<double> > run16parBoundsGauMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16parBoundsGauMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16parBoundsGauMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsGauMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run16sGauMassPtPol1={
    {"run16aGauPt1Pol1","run16bGauPt1Pol1","run16cGauPt1Pol1","run16 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPt2Pol1","run16bGauPt2Pol1","run16cGauPt2Pol1","run16 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPt3Pol1","run16bGauPt3Pol1","run16cGauPt3Pol1","run16 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPt4Pol1","run16bGauPt4Pol1","run16cGauPt4Pol1","run16 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPt5Pol1","run16bGauPt5Pol1","run16cGauPt5Pol1","run16 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPt6Pol1","run16bGauPt6Pol1","run16cGauPt6Pol1","run16 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run16sGauMassPtPol2={
    {"run16aGauPt1Pol2","run16bGauPt1Pol2","run16cGauPt1Pol2","run16 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"run16aGauPt2Pol2","run16bGauPt2Pol2","run16cGauPt2Pol2","run16 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"run16aGauPt3Pol2","run16bGauPt3Pol2","run16cGauPt3Pol2","run16 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"run16aGauPt4Pol2","run16bGauPt4Pol2","run16cGauPt4Pol2","run16 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"run16aGauPt5Pol2","run16bGauPt5Pol2","run16cGauPt5Pol2","run16 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"run16aGauPt6Pol2","run16bGauPt6Pol2","run16cGauPt6Pol2","run16 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run16sGauMassPtExp={
    {"run16aGauPt1exp","run16bGauPt1exp","run16cGauPt1exp","run16 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPt2exp","run16bGauPt2exp","run16cGauPt2exp","run16 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPt3exp","run16bGauPt3exp","run16cGauPt3exp","run16 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPt4exp","run16bGauPt4exp","run16cGauPt4exp","run16 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPt5exp","run16bGauPt5exp","run16cGauPt5exp","run16 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPt6exp","run16bGauPt6exp","run16cGauPt6exp","run16 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run16sGauMassPhiPol1lpt={
    {"run16aGauPhi1Pol1lpt","run16bGauPhi1Pol1lpt","run16cGauPhi1Pol1lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi2Pol1lpt","run16bGauPhi2Pol1lpt","run16cGauPhi2Pol1lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi3Pol1lpt","run16bGauPhi3Pol1lpt","run16cGauPhi3Pol1lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi4Pol1lpt","run16bGauPhi4Pol1lpt","run16cGauPhi4Pol1lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi5Pol1lpt","run16bGauPhi5Pol1lpt","run16cGauPhi5Pol1lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run16sGauMassPhiPol2lpt={
    {"run16aGauPhi1Pol2lpt","run16bGauPhi1Pol2lpt","run16cGauPhi1Pol2lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi2Pol2lpt","run16bGauPhi2Pol2lpt","run16cGauPhi2Pol2lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi3Pol2lpt","run16bGauPhi3Pol2lpt","run16cGauPhi3Pol2lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi4Pol2lpt","run16bGauPhi4Pol2lpt","run16cGauPhi4Pol2lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi5Pol2lpt","run16bGauPhi5Pol2lpt","run16cGauPhi5Pol2lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16sGauMassPhiExplpt={
    {"run16aGauPhi1Explpt","run16bGauPhi1Explpt","run16cGauPhi1Explpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi2Explpt","run16bGauPhi2Explpt","run16cGauPhi2Explpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi3Explpt","run16bGauPhi3Explpt","run16cGauPhi3Explpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi4Explpt","run16bGauPhi4Explpt","run16cGauPhi4Explpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi5Explpt","run16bGauPhi5Explpt","run16cGauPhi5Explpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run16sGauMassPhiPol1hpt={
    {"run16aGauPhi1Pol1hpt","run16bGauPhi1Pol1hpt","run16cGauPhi1Pol1hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi2Pol1hpt","run16bGauPhi2Pol1hpt","run16cGauPhi2Pol1hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi3Pol1hpt","run16bGauPhi3Pol1hpt","run16cGauPhi3Pol1hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi4Pol1hpt","run16bGauPhi4Pol1hpt","run16cGauPhi4Pol1hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"run16aGauPhi5Pol1hpt","run16bGauPhi5Pol1hpt","run16cGauPhi5Pol1hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run16sGauMassPhiPol2hpt={
    {"run16aGauPhi1Pol2hpt","run16bGauPhi1Pol2hpt","run16cGauPhi1Pol2hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi2Pol2hpt","run16bGauPhi2Pol2hpt","run16cGauPhi2Pol2hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi3Pol2hpt","run16bGauPhi3Pol2hpt","run16cGauPhi3Pol2hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi4Pol2hpt","run16bGauPhi4Pol2hpt","run16cGauPhi4Pol2hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"run16aGauPhi5Pol2hpt","run16bGauPhi5Pol2hpt","run16cGauPhi5Pol2hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16sGauMassPhiExphpt={
    {"run16aGauPhi1Exphpt","run16bGauPhi1Exphpt","run16cGauPhi1Exphpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi2Exphpt","run16bGauPhi2Exphpt","run16cGauPhi2Exphpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi3Exphpt","run16bGauPhi3Exphpt","run16cGauPhi3Exphpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi4Exphpt","run16bGauPhi4Exphpt","run16cGauPhi4Exphpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"run16aGauPhi5Exphpt","run16bGauPhi5Exphpt","run16cGauPhi5Exphpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};

vector < vector<TString> > run16v2GaussNames={
    {"run 16 dN/dPhi vs Phi (gaussian lin)  0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
    {"run 16 dN/dPhi vs Phi (gaussian quad) 0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
    {"run 16 dN/dPhi vs Phi (gaussian exp)  0<pT<2GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
    {"run 16 dN/dPhi vs Phi (gaussian lin)  2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
    {"run 16 dN/dPhi vs Phi (gaussian quad) 2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"},
    {"run 16 dN/dPhi vs Phi (gaussian exp)  2<pT<5GeV","dN/dPhi","dphi",dir+"v2plotGauss"}
};

///// Run 16 Breit Wig

vector < vector<double> > run16parBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e5,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2.5e3,1e4,0.97,1.00,15e-3,40e-3,1e6,1e8,-1e8,-5e7,1e6,5e7},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiPol2lpt={
    {1e3,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16parBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run16sBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","run16 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","run16 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","run16 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","run16 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","run16 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","run16 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run16sBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","run16 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","run16 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","run16 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","run16 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","run16 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","run16 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run16sBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","run16 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","run16 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","run16 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","run16 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","run16 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","run16 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run16sBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run16sBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16sBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run16sBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run16sBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16sBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};

// Run14

vector < vector<double> > run14parBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1,1e5,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2,1e4,0.97,1.00,15e-3,40e-3,1e6,1e8,-1e8,-5e7,1e6,5e7},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiPol2lpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run14parBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,15e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run14sBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","run14 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","run14 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","run14 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","run14 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","run14 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","run14 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run14sBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","run14 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","run14 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","run14 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","run14 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","run14 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","run14 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run14sBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","run14 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","run14 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","run14 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","run14 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","run14 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","run14 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run14sBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","run14 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","run14 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","run14 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","run14 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","run14 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run14sBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","run14 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","run14 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","run14 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","run14 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","run14 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run14sBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","run14 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","run14 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","run14 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","run14 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","run14 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run14sBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","run14 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","run14 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","run14 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","run14 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","run14 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run14sBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","run14 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","run14 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","run14 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","run14 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","run14 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run14sBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","run14 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","run14 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","run14 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","run14 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","run14 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};





////////////// BEGINNING OF RUN11 + RUN16 FITS /////////////////

vector < vector<double> > run11full16partparBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e5,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2.5e3,1e4,0.97,1.00,15e-3,70e-3,1e6,1e8,-1e9,-5e7,1e6,5e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {100,1e4,0.98,1.0,25e-3,55e-3,-1e8,1e8,1e3,1e8,-5e6,-1e3},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e3,0.97,0.99,25e-3,40e-3,-1e7,1e7,1e5,1e7,-1e7,-1e5}
};
vector < vector<double> > run11full16partparBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11full16partsBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","runs 11 and 16 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","runs 11 and 16 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","runs 11 and 16 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","runs 11 and 16 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","runs 11 and 16 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","runs 11 and 16 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run11full16partsBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","runs 11 and 16 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","runs 11 and 16 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","runs 11 and 16 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","runs 11 and 16 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","runs 11 and 16 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","runs 11 and 16 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11full16partsBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","runs 11 and 16 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","runs 11 and 16 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","runs 11 and 16 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","runs 11 and 16 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","runs 11 and 16 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","runs 11 and 16 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11full16partsBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run11full16partsBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11full16partsBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11full16partsBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run11full16partsBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11full16partsBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","runs 11 and 16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","runs 11 and 16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","runs 11 and 16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","runs 11 and 16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","runs 11 and 16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};

// some test

    // run 11 test
vector < vector<double> > run11testparBoundsPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.97,1.0,5e-3,20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {3e3,4e3,0.98,0.99,5e-3, 20e-3,-1e8,7e7,-1e8,-7.1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.97,1.00,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11testparBoundsPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run11testparBoundsPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run11testparBoundsPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run11testStrPtPol1={
    {"aGauPt1Pol1","bGauPt1Pol1","cGauPt1Pol1","run11 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt2Pol1","bGauPt2Pol1","cGauPt2Pol1","run11 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt3Pol1","bGauPt3Pol1","cGauPt3Pol1","run11 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt4Pol1","bGauPt4Pol1","cGauPt4Pol1","run11 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt5Pol1","bGauPt5Pol1","cGauPt5Pol1","run11 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt6Pol1","bGauPt6Pol1","cGauPt6Pol1","run11 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run11testStrPtPol2={
    {"aGauPt1Pol2","bGauPt1Pol2","cGauPt1Pol2","run11 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt2Pol2","bGauPt2Pol2","cGauPt2Pol2","run11 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt3Pol2","bGauPt3Pol2","cGauPt3Pol2","run11 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt4Pol2","bGauPt4Pol2","cGauPt4Pol2","run11 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt5Pol2","bGauPt5Pol2","cGauPt5Pol2","run11 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt6Pol2","bGauPt6Pol2","cGauPt6Pol2","run11 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run11testStrPtExp={
    {"aGauPt1exp","bGauPt1exp","cGauPt1exp","run11 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt2exp","bGauPt2exp","cGauPt2exp","run11 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt3exp","bGauPt3exp","cGauPt3exp","run11 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt4exp","bGauPt4exp","cGauPt4exp","run11 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt5exp","bGauPt5exp","cGauPt5exp","run11 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt6exp","bGauPt6exp","cGauPt6exp","run11 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run11testStrPhiPol1lpt={
    {"aGauPhi1Pol1lpt","bGauPhi1Pol1lpt","cGauPhi1Pol1lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1lpt","bGauPhi2Pol1lpt","cGauPhi2Pol1lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1lpt","bGauPhi3Pol1lpt","cGauPhi3Pol1lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1lpt","bGauPhi4Pol1lpt","cGauPhi4Pol1lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1lpt","bGauPhi5Pol1lpt","cGauPhi5Pol1lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run11testStrPhiPol2lpt={
    {"aGauPhi1Pol2lpt","bGauPhi1Pol2lpt","cGauPhi1Pol2lpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2lpt","bGauPhi2Pol2lpt","cGauPhi2Pol2lpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2lpt","bGauPhi3Pol2lpt","cGauPhi3Pol2lpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2lpt","bGauPhi4Pol2lpt","cGauPhi4Pol2lpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2lpt","bGauPhi5Pol2lpt","cGauPhi5Pol2lpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11testStrPhiExplpt={
    {"aGauPhi1Explpt","bGauPhi1Explpt","cGauPhi1Explpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Explpt","bGauPhi2Explpt","cGauPhi2Explpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Explpt","bGauPhi3Explpt","cGauPhi3Explpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Explpt","bGauPhi4Explpt","cGauPhi4Explpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Explpt","bGauPhi5Explpt","cGauPhi5Explpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run11testStrPhiPol1hpt={
    {"aGauPhi1Pol1hpt","bGauPhi1Pol1hpt","cGauPhi1Pol1hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1hpt","bGauPhi2Pol1hpt","cGauPhi2Pol1hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1hpt","bGauPhi3Pol1hpt","cGauPhi3Pol1hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1hpt","bGauPhi4Pol1hpt","cGauPhi4Pol1hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1hpt","bGauPhi5Pol1hpt","cGauPhi5Pol1hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run11testStrPhiPol2hpt={
    {"aGauPhi1Pol2hpt","bGauPhi1Pol2hpt","cGauPhi1Pol2hpt","run11 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2hpt","bGauPhi2Pol2hpt","cGauPhi2Pol2hpt","run11 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2hpt","bGauPhi3Pol2hpt","cGauPhi3Pol2hpt","run11 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2hpt","bGauPhi4Pol2hpt","cGauPhi4Pol2hpt","run11 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2hpt","bGauPhi5Pol2hpt","cGauPhi5Pol2hpt","run11 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run11testStrPhiExphpt={
    {"aGauPhi1Exphpt","bGauPhi1Exphpt","cGauPhi1Exphpt","run11 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Exphpt","bGauPhi2Exphpt","cGauPhi2Exphpt","run11 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Exphpt","bGauPhi3Exphpt","cGauPhi3Exphpt","run11 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Exphpt","bGauPhi4Exphpt","cGauPhi4Exphpt","run11 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Exphpt","bGauPhi5Exphpt","cGauPhi5Exphpt","run11 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
    // run 16 test
vector < vector<double> > run16testparBoundsPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.97,1.0,5e-3,20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {3e3,4e3,0.98,0.99,5e-3, 20e-3,-1e8,7e7,-1e8,-7.1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.97,1.00,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16testparBoundsPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16testparBoundsPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16testparBoundsPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run16testStrPtPol1={
    {"aGauPt1Pol1","bGauPt1Pol1","cGauPt1Pol1","run16 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt2Pol1","bGauPt2Pol1","cGauPt2Pol1","run16 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt3Pol1","bGauPt3Pol1","cGauPt3Pol1","run16 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt4Pol1","bGauPt4Pol1","cGauPt4Pol1","run16 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt5Pol1","bGauPt5Pol1","cGauPt5Pol1","run16 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt6Pol1","bGauPt6Pol1","cGauPt6Pol1","run16 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run16testStrPtPol2={
    {"aGauPt1Pol2","bGauPt1Pol2","cGauPt1Pol2","run16 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt2Pol2","bGauPt2Pol2","cGauPt2Pol2","run16 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt3Pol2","bGauPt3Pol2","cGauPt3Pol2","run16 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt4Pol2","bGauPt4Pol2","cGauPt4Pol2","run16 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt5Pol2","bGauPt5Pol2","cGauPt5Pol2","run16 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt6Pol2","bGauPt6Pol2","cGauPt6Pol2","run16 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run16testStrPtExp={
    {"aGauPt1exp","bGauPt1exp","cGauPt1exp","run16 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt2exp","bGauPt2exp","cGauPt2exp","run16 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt3exp","bGauPt3exp","cGauPt3exp","run16 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt4exp","bGauPt4exp","cGauPt4exp","run16 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt5exp","bGauPt5exp","cGauPt5exp","run16 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt6exp","bGauPt6exp","cGauPt6exp","run16 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run16testStrPhiPol1lpt={
    {"aGauPhi1Pol1lpt","bGauPhi1Pol1lpt","cGauPhi1Pol1lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1lpt","bGauPhi2Pol1lpt","cGauPhi2Pol1lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1lpt","bGauPhi3Pol1lpt","cGauPhi3Pol1lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1lpt","bGauPhi4Pol1lpt","cGauPhi4Pol1lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1lpt","bGauPhi5Pol1lpt","cGauPhi5Pol1lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run16testStrPhiPol2lpt={
    {"aGauPhi1Pol2lpt","bGauPhi1Pol2lpt","cGauPhi1Pol2lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2lpt","bGauPhi2Pol2lpt","cGauPhi2Pol2lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2lpt","bGauPhi3Pol2lpt","cGauPhi3Pol2lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2lpt","bGauPhi4Pol2lpt","cGauPhi4Pol2lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2lpt","bGauPhi5Pol2lpt","cGauPhi5Pol2lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16testStrPhiExplpt={
    {"aGauPhi1Explpt","bGauPhi1Explpt","cGauPhi1Explpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Explpt","bGauPhi2Explpt","cGauPhi2Explpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Explpt","bGauPhi3Explpt","cGauPhi3Explpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Explpt","bGauPhi4Explpt","cGauPhi4Explpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Explpt","bGauPhi5Explpt","cGauPhi5Explpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run16testStrPhiPol1hpt={
    {"aGauPhi1Pol1hpt","bGauPhi1Pol1hpt","cGauPhi1Pol1hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1hpt","bGauPhi2Pol1hpt","cGauPhi2Pol1hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1hpt","bGauPhi3Pol1hpt","cGauPhi3Pol1hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1hpt","bGauPhi4Pol1hpt","cGauPhi4Pol1hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1hpt","bGauPhi5Pol1hpt","cGauPhi5Pol1hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run16testStrPhiPol2hpt={
    {"aGauPhi1Pol2hpt","bGauPhi1Pol2hpt","cGauPhi1Pol2hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2hpt","bGauPhi2Pol2hpt","cGauPhi2Pol2hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2hpt","bGauPhi3Pol2hpt","cGauPhi3Pol2hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2hpt","bGauPhi4Pol2hpt","cGauPhi4Pol2hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2hpt","bGauPhi5Pol2hpt","cGauPhi5Pol2hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16testStrPhiExphpt={
    {"aGauPhi1Exphpt","bGauPhi1Exphpt","cGauPhi1Exphpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Exphpt","bGauPhi2Exphpt","cGauPhi2Exphpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Exphpt","bGauPhi3Exphpt","cGauPhi3Exphpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Exphpt","bGauPhi4Exphpt","cGauPhi4Exphpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Exphpt","bGauPhi5Exphpt","cGauPhi5Exphpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};

// run 16 partial data momentum correction
vector < vector<double> > run16momCorrparBoundsPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.97,1.0,5e-3,20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {3e3,4e3,0.98,0.99,5e-3, 20e-3,-1e8,7e7,-1e8,-7.1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.97,1.00,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16momCorrparBoundsPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run16momCorrparBoundsPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run16momCorrparBoundsPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run16momCorrStrPtPol1={
    {"aGauPt1Pol1","bGauPt1Pol1","cGauPt1Pol1","run16 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt2Pol1","bGauPt2Pol1","cGauPt2Pol1","run16 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt3Pol1","bGauPt3Pol1","cGauPt3Pol1","run16 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt4Pol1","bGauPt4Pol1","cGauPt4Pol1","run16 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt5Pol1","bGauPt5Pol1","cGauPt5Pol1","run16 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt6Pol1","bGauPt6Pol1","cGauPt6Pol1","run16 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run16momCorrStrPtPol2={
    {"aGauPt1Pol2","bGauPt1Pol2","cGauPt1Pol2","run16 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt2Pol2","bGauPt2Pol2","cGauPt2Pol2","run16 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt3Pol2","bGauPt3Pol2","cGauPt3Pol2","run16 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt4Pol2","bGauPt4Pol2","cGauPt4Pol2","run16 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt5Pol2","bGauPt5Pol2","cGauPt5Pol2","run16 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt6Pol2","bGauPt6Pol2","cGauPt6Pol2","run16 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run16momCorrStrPtExp={
    {"aGauPt1exp","bGauPt1exp","cGauPt1exp","run16 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt2exp","bGauPt2exp","cGauPt2exp","run16 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt3exp","bGauPt3exp","cGauPt3exp","run16 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt4exp","bGauPt4exp","cGauPt4exp","run16 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt5exp","bGauPt5exp","cGauPt5exp","run16 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt6exp","bGauPt6exp","cGauPt6exp","run16 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run16momCorrStrPhiPol1lpt={
    {"aGauPhi1Pol1lpt","bGauPhi1Pol1lpt","cGauPhi1Pol1lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1lpt","bGauPhi2Pol1lpt","cGauPhi2Pol1lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1lpt","bGauPhi3Pol1lpt","cGauPhi3Pol1lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1lpt","bGauPhi4Pol1lpt","cGauPhi4Pol1lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1lpt","bGauPhi5Pol1lpt","cGauPhi5Pol1lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run16momCorrStrPhiPol2lpt={
    {"aGauPhi1Pol2lpt","bGauPhi1Pol2lpt","cGauPhi1Pol2lpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2lpt","bGauPhi2Pol2lpt","cGauPhi2Pol2lpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2lpt","bGauPhi3Pol2lpt","cGauPhi3Pol2lpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2lpt","bGauPhi4Pol2lpt","cGauPhi4Pol2lpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2lpt","bGauPhi5Pol2lpt","cGauPhi5Pol2lpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16momCorrStrPhiExplpt={
    {"aGauPhi1Explpt","bGauPhi1Explpt","cGauPhi1Explpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Explpt","bGauPhi2Explpt","cGauPhi2Explpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Explpt","bGauPhi3Explpt","cGauPhi3Explpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Explpt","bGauPhi4Explpt","cGauPhi4Explpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Explpt","bGauPhi5Explpt","cGauPhi5Explpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run16momCorrStrPhiPol1hpt={
    {"aGauPhi1Pol1hpt","bGauPhi1Pol1hpt","cGauPhi1Pol1hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1hpt","bGauPhi2Pol1hpt","cGauPhi2Pol1hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1hpt","bGauPhi3Pol1hpt","cGauPhi3Pol1hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1hpt","bGauPhi4Pol1hpt","cGauPhi4Pol1hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1hpt","bGauPhi5Pol1hpt","cGauPhi5Pol1hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run16momCorrStrPhiPol2hpt={
    {"aGauPhi1Pol2hpt","bGauPhi1Pol2hpt","cGauPhi1Pol2hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2hpt","bGauPhi2Pol2hpt","cGauPhi2Pol2hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2hpt","bGauPhi3Pol2hpt","cGauPhi3Pol2hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2hpt","bGauPhi4Pol2hpt","cGauPhi4Pol2hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2hpt","bGauPhi5Pol2hpt","cGauPhi5Pol2hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run16momCorrStrPhiExphpt={
    {"aGauPhi1Exphpt","bGauPhi1Exphpt","cGauPhi1Exphpt","run16 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Exphpt","bGauPhi2Exphpt","cGauPhi2Exphpt","run16 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Exphpt","bGauPhi3Exphpt","cGauPhi3Exphpt","run16 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Exphpt","bGauPhi4Exphpt","cGauPhi4Exphpt","run16 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Exphpt","bGauPhi5Exphpt","cGauPhi5Exphpt","run16 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};

vector < vector<double> > run10firstparBoundsPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e4,1e8,0.97,1.0,5e-3,20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {3e3,4e3,0.98,0.99,5e-3, 20e-3,-1e8,7e7,-1e8,-7.1e5},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPtPol2={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.97,1.00,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.98,0.99,5e-3,20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {1e4,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {50,1e8,0.98,0.99,5e-3, 20e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run10firstparBoundsPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run10firstparBoundsPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPhiPol2hpt={
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 25e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3, 20e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run10firstparBoundsPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 19e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3, 30e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

vector < vector<TString> > run10firstStrPtPol1={
    {"aGauPt1Pol1","bGauPt1Pol1","cGauPt1Pol1","run10 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt2Pol1","bGauPt2Pol1","cGauPt2Pol1","run10 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt3Pol1","bGauPt3Pol1","cGauPt3Pol1","run10 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt4Pol1","bGauPt4Pol1","cGauPt4Pol1","run10 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt5Pol1","bGauPt5Pol1","cGauPt5Pol1","run10 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPt6Pol1","bGauPt6Pol1","cGauPt6Pol1","run10 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};

vector < vector<TString> > run10firstStrPtPol2={
    {"aGauPt1Pol2","bGauPt1Pol2","cGauPt1Pol2","run10 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt2Pol2","bGauPt2Pol2","cGauPt2Pol2","run10 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt3Pol2","bGauPt3Pol2","cGauPt3Pol2","run10 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt4Pol2","bGauPt4Pol2","cGauPt4Pol2","run10 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt5Pol2","bGauPt5Pol2","cGauPt5Pol2","run10 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"},
    {"aGauPt6Pol2","bGauPt6Pol2","cGauPt6Pol2","run10 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run10firstStrPtExp={
    {"aGauPt1exp","bGauPt1exp","cGauPt1exp","run10 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt2exp","bGauPt2exp","cGauPt2exp","run10 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt3exp","bGauPt3exp","cGauPt3exp","run10 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt4exp","bGauPt4exp","cGauPt4exp","run10 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt5exp","bGauPt5exp","cGauPt5exp","run10 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPt6exp","bGauPt6exp","cGauPt6exp","run10 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run10firstStrPhiPol1lpt={
    {"aGauPhi1Pol1lpt","bGauPhi1Pol1lpt","cGauPhi1Pol1lpt","run10 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1lpt","bGauPhi2Pol1lpt","cGauPhi2Pol1lpt","run10 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1lpt","bGauPhi3Pol1lpt","cGauPhi3Pol1lpt","run10 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1lpt","bGauPhi4Pol1lpt","cGauPhi4Pol1lpt","run10 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1lpt","bGauPhi5Pol1lpt","cGauPhi5Pol1lpt","run10 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}
};
vector < vector<TString> > run10firstStrPhiPol2lpt={
    {"aGauPhi1Pol2lpt","bGauPhi1Pol2lpt","cGauPhi1Pol2lpt","run10 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2lpt","bGauPhi2Pol2lpt","cGauPhi2Pol2lpt","run10 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2lpt","bGauPhi3Pol2lpt","cGauPhi3Pol2lpt","run10 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2lpt","bGauPhi4Pol2lpt","cGauPhi4Pol2lpt","run10 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2lpt","bGauPhi5Pol2lpt","cGauPhi5Pol2lpt","run10 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run10firstStrPhiExplpt={
    {"aGauPhi1Explpt","bGauPhi1Explpt","cGauPhi1Explpt","run10 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Explpt","bGauPhi2Explpt","cGauPhi2Explpt","run10 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Explpt","bGauPhi3Explpt","cGauPhi3Explpt","run10 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Explpt","bGauPhi4Explpt","cGauPhi4Explpt","run10 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Explpt","bGauPhi5Explpt","cGauPhi5Explpt","run10 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run10firstStrPhiPol1hpt={
    {"aGauPhi1Pol1hpt","bGauPhi1Pol1hpt","cGauPhi1Pol1hpt","run16 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi2Pol1hpt","bGauPhi2Pol1hpt","cGauPhi2Pol1hpt","run16 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi3Pol1hpt","bGauPhi3Pol1hpt","cGauPhi3Pol1hpt","run16 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi4Pol1hpt","bGauPhi4Pol1hpt","cGauPhi4Pol1hpt","run16 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"},
    {"aGauPhi5Pol1hpt","bGauPhi5Pol1hpt","cGauPhi5Pol1hpt","run16 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1"}

};
vector < vector<TString> > run10firstStrPhiPol2hpt={
    {"aGauPhi1Pol2hpt","bGauPhi1Pol2hpt","cGauPhi1Pol2hpt","run10 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi2Pol2hpt","bGauPhi2Pol2hpt","cGauPhi2Pol2hpt","run10 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi3Pol2hpt","bGauPhi3Pol2hpt","cGauPhi3Pol2hpt","run10 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi4Pol2hpt","bGauPhi4Pol2hpt","cGauPhi4Pol2hpt","run10 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"},
    {"aGauPhi5Pol2hpt","bGauPhi5Pol2hpt","cGauPhi5Pol2hpt","run10 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run10firstStrPhiExphpt={
    {"aGauPhi1Exphpt","bGauPhi1Exphpt","cGauPhi1Exphpt","run10 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi2Exphpt","bGauPhi2Exphpt","cGauPhi2Exphpt","run10 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi3Exphpt","bGauPhi3Exphpt","cGauPhi3Exphpt","run10 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi4Exphpt","bGauPhi4Exphpt","cGauPhi4Exphpt","run10 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"},
    {"aGauPhi5Exphpt","bGauPhi5Exphpt","cGauPhi5Exphpt","run10 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, gaussian","mass GeV","counts",dir+"saves.pdf","area","mass","width","offset","exp amp","decay"}
};

vector < vector<double> > run12half_parBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e2,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2.5e2,1e4,0.97,1.00,15e-3,40e-3,1e6,1e8,-1e8,-5e7,1e6,5e7},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1e2,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,60e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12half_parBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

// shouldn't have to worry about these names too much, as the objects should be deleted after primary->Analysis() is completed
vector < vector<TString> > run12half_sBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","run12 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","run12 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","run12 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","run12 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","run12 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","run12 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run12half_sBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","run12 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","run12 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","run12 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","run12 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","run12 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","run12 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run12half_sBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","run12 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","run12 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","run12 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","run12 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","run12 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","run12 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run12half_sBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","run12 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","run12 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","run12 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","run12 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","run12 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run12half_sBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","run12 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","run12 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","run12 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","run12 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","run12 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run12half_sBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","run12 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","run12 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","run12 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","run12 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","run12 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run12half_sBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","run12 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","run12 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","run12 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","run12 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","run12 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run12half_sBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","run12 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","run12 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","run12 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","run12 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","run12 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run12half_sBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","run12 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","run12 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","run12 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","run12 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","run12 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};

vector < vector<double> > run12parBoundsBreitWigMassPtPol1={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e2,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPtPol2={
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e1,1e5,0.95,1.05,10e-3,40e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e5,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {2.5e2,1e4,0.97,1.00,15e-3,40e-3,1e6,1e8,-1e8,-5e7,1e6,5e7},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {50,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPtExp={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {1e2,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiPol1lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiPol2lpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {1e3,1e8,0.95,1.05,5e-3,60e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiExplpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiPol1hpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiPol2hpt={
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8},
    {0,1e8,0.95,1.05,25e-3,100e-3,-1e8,1e8,-1e8,1e8,-1e8,1e8}
};
vector < vector<double> > run12parBoundsBreitWigMassPhiExphpt={
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.98,1.00,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10},
    {0,1e8,0.95,1.05,5e-3,100e-3,-1e8,1e8,0,1e8,0,10}
};

// shouldn't have to worry about these names too much, as the objects should be deleted after primary->Analysis() is completed
vector < vector<TString> > run12sBreitWigMassPtPol1={
    {"aBreitWigPt1Pol1","bBreitWigPt1Pol1","cBreitWigPt1Pol1","run12 Mass distribution: 0.0<pT<0.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt2Pol1","bBreitWigPt2Pol1","cBreitWigPt2Pol1","run12 Mass distribution: 0.5<pT<1.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt3Pol1","bBreitWigPt3Pol1","cBreitWigPt3Pol1","run12 Mass distribution: 1.0<pT<1.5 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt4Pol1","bBreitWigPt4Pol1","cBreitWigPt4Pol1","run12 Mass distribution: 1.5<pT<2.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt5Pol1","bBreitWigPt5Pol1","cBreitWigPt5Pol1","run12 Mass distribution: 2.0<pT<3.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPt6Pol1","bBreitWigPt6Pol1","cBreitWigPt6Pol1","run12 Mass distribution: 3.0<pT<5.0 GeV (linear bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};

vector < vector<TString> > run12sBreitWigMassPtPol2={
    {"aBreitWigPt1Pol2","bBreitWigPt1Pol2","cBreitWigPt1Pol2","run12 Mass distribution: 0.0<pT<0.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt2Pol2","bBreitWigPt2Pol2","cBreitWigPt2Pol2","run12 Mass distribution: 0.5<pT<1.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt3Pol2","bBreitWigPt3Pol2","cBreitWigPt3Pol2","run12 Mass distribution: 1.0<pT<1.5 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt4Pol2","bBreitWigPt4Pol2","cBreitWigPt4Pol2","run12 Mass distribution: 1.5<pT<2.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt5Pol2","bBreitWigPt5Pol2","cBreitWigPt5Pol2","run12 Mass distribution: 2.0<pT<3.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"},
    {"aBreitWigPt6Pol2","bBreitWigPt6Pol2","cBreitWigPt6Pol2","run12 Mass distribution: 3.0<pT<5.0 GeV (quadratic bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly 2"}
};
vector < vector<TString> > run12sBreitWigMassPtExp={
    {"aBreitWigPt1exp","bBreitWigPt1exp","cBreitWigPt1exp","run12 Mass distribution: 0.0<pT<0.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt2exp","bBreitWigPt2exp","cBreitWigPt2exp","run12 Mass distribution: 0.5<pT<1.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt3exp","bBreitWigPt3exp","cBreitWigPt3exp","run12 Mass distribution: 1.0<pT<1.5 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt4exp","bBreitWigPt4exp","cBreitWigPt4exp","run12 Mass distribution: 1.5<pT<2.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt5exp","bBreitWigPt5exp","cBreitWigPt5exp","run12 Mass distribution: 2.0<pT<3.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPt6exp","bBreitWigPt6exp","cBreitWigPt6exp","run12 Mass distribution: 3.0<pT<5.0 GeV (exponential bkg), Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 0<pt<2 GeV
vector < vector<TString> > run12sBreitWigMassPhiPol1lpt={
    {"aBreitWigPhi1Pol1lpt","bBreitWigPhi1Pol1lpt","cBreitWigPhi1Pol1lpt","run12 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1lpt","bBreitWigPhi2Pol1lpt","cBreitWigPhi2Pol1lpt","run12 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1lpt","bBreitWigPhi3Pol1lpt","cBreitWigPhi3Pol1lpt","run12 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1lpt","bBreitWigPhi4Pol1lpt","cBreitWigPhi4Pol1lpt","run12 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1lpt","bBreitWigPhi5Pol1lpt","cBreitWigPhi5Pol1lpt","run12 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}
};
vector < vector<TString> > run12sBreitWigMassPhiPol2lpt={
    {"aBreitWigPhi1Pol2lpt","bBreitWigPhi1Pol2lpt","cBreitWigPhi1Pol2lpt","run12 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2lpt","bBreitWigPhi2Pol2lpt","cBreitWigPhi2Pol2lpt","run12 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2lpt","bBreitWigPhi3Pol2lpt","cBreitWigPhi3Pol2lpt","run12 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2lpt","bBreitWigPhi4Pol2lpt","cBreitWigPhi4Pol2lpt","run12 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2lpt","bBreitWigPhi5Pol2lpt","cBreitWigPhi5Pol2lpt","run12 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run12sBreitWigMassPhiExplpt={
    {"aBreitWigPhi1Explpt","bBreitWigPhi1Explpt","cBreitWigPhi1Explpt","run12 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Explpt","bBreitWigPhi2Explpt","cBreitWigPhi2Explpt","run12 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Explpt","bBreitWigPhi3Explpt","cBreitWigPhi3Explpt","run12 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Explpt","bBreitWigPhi4Explpt","cBreitWigPhi4Explpt","run12 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Explpt","bBreitWigPhi5Explpt","cBreitWigPhi5Explpt","run12 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 0<pT<2GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};
        // phi distribution 2<pt<5 GeV
vector < vector<TString> > run12sBreitWigMassPhiPol1hpt={
    {"aBreitWigPhi1Pol1hpt","bBreitWigPhi1Pol1hpt","cBreitWigPhi1Pol1hpt","run12 Mass distribution: 0.0pi<phi<0.1pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi2Pol1hpt","bBreitWigPhi2Pol1hpt","cBreitWigPhi2Pol1hpt","run12 Mass distribution: 0.1pi<phi<0.2pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi3Pol1hpt","bBreitWigPhi3Pol1hpt","cBreitWigPhi3Pol1hpt","run12 Mass distribution: 0.2pi<phi<0.3pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi4Pol1hpt","bBreitWigPhi4Pol1hpt","cBreitWigPhi4Pol1hpt","run12 Mass distribution: 0.3pi<phi<0.4pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"},
    {"aBreitWigPhi5Pol1hpt","bBreitWigPhi5Pol1hpt","cBreitWigPhi5Pol1hpt","run12 Mass distribution: 0.4pi<phi<0.5pi (linear bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1"}

};
vector < vector<TString> > run12sBreitWigMassPhiPol2hpt={
    {"aBreitWigPhi1Pol2hpt","bBreitWigPhi1Pol2hpt","cBreitWigPhi1Pol2hpt","run12 Mass distribution: 0.0pi<phi<0.1pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi2Pol2hpt","bBreitWigPhi2Pol2hpt","cBreitWigPhi2Pol2hpt","run12 Mass distribution: 0.1pi<phi<0.2pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi3Pol2hpt","bBreitWigPhi3Pol2hpt","cBreitWigPhi3Pol2hpt","run12 Mass distribution: 0.2pi<phi<0.3pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi4Pol2hpt","bBreitWigPhi4Pol2hpt","cBreitWigPhi4Pol2hpt","run12 Mass distribution: 0.3pi<phi<0.4pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"},
    {"aBreitWigPhi5Pol2hpt","bBreitWigPhi5Pol2hpt","cBreitWigPhi5Pol2hpt","run12 Mass distribution: 0.4pi<phi<0.5pi (quadratic bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","poly 1","poly2"}
};
vector < vector<TString> > run12sBreitWigMassPhiExphpt={
    {"aBreitWigPhi1Exphpt","bBreitWigPhi1Exphpt","cBreitWigPhi1Exphpt","run12 Mass distribution: 0.0pi<phi<0.1pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi2Exphpt","bBreitWigPhi2Exphpt","cBreitWigPhi2Exphpt","run12 Mass distribution: 0.1pi<phi<0.2pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi3Exphpt","bBreitWigPhi3Exphpt","cBreitWigPhi3Exphpt","run12 Mass distribution: 0.2pi<phi<0.3pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi4Exphpt","bBreitWigPhi4Exphpt","cBreitWigPhi4Exphpt","run12 Mass distribution: 0.3pi<phi<0.4pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"},
    {"aBreitWigPhi5Exphpt","bBreitWigPhi5Exphpt","cBreitWigPhi5Exphpt","run12 Mass distribution: 0.4pi<phi<0.5pi (exponential bkg) 2<pT<5GeV, Breit Wigner","mass GeV","counts",dir+"saves.pdf","Breit amp","mass","width","offset","exp amp","decay"}
};



#endif
