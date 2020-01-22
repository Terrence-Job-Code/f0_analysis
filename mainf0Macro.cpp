#include "rootheader.h"
#include "primary.h"
#include "primary.cpp"
#include "fits.h"
#include "fits.cpp"
#include "auxFun.cpp"

using namespace std;

void mainf0Macro()
{

    // some gStyle options for plotting

    //gStyle->SetOptFit(1111);

    gStyle->SetOptTitle(0);
	gStyle->SetPadTickX(0);
	gStyle->SetPadTickY(0);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptDate(0);
	gStyle->SetErrorX(0);
	gStyle->SetEndErrorSize(8);


    TString f11 = "DATA/OLDDATA/run11Data.root";
   
    TString f16_primary = "DATA/Run16_primary/run16data200_primary.root"; 
    TString f16_global = "DATA/Run16_global/run16data200_global.root"; 

    TString f14_primary = "DATA/Run14_nohft_primary/run14data200_primary.root"; 
    TString f14_global = "DATA/Run14_nohft_global/run14data200_global.root"; 

    TString f10_first = "DATA/Run10_set01/run10Data.root";

    //TString f12_half = "DATA/UURun12/halfData/run12UUhalf.root"; 
    TString f12_half = "DATA/UURun12/fullFirst/run12UU.root"; 
    //TString f16 = "DATA/OLDDATA/run16Data.root"; // this might not have event plane correction
    //TString f16 = "DATA/Run16_full_incorrect/run16Full.root"; 
    //TString f16 = "DATA/Run16_full_correct/run16Full.root"; 
 
    //TString f1611 = "DATA/OLDDATA/run16_and_11Data.root";
    //TString f11_part16 = "DATA/OLDDATA/run11full16third.root";

    //TString f16_notof = "DATA/Run16_notof_half/run16_50perData.root";

    //TString f11_test = "DATA/Run11_set01/run11Data_set01.root"; // ahh, this is a test run
    //TString f16_test = "DATA/Run16_set01/run16Data_set01.root";

    //TString f16_momCorrPart = "DATA/Run16_momCorr/run16_momCorrect.root";


    //TString f14_att01 = "DATA/Run14_noHFT_02/run14Data200.root";
    //TString f14_att01 = "DATA/Run14_noHFT_03/run14Data200.root"; // these are global tracks, somehow I missed that

    // RUN11 Gaussian analysis

    
    primaryFit *tempGau = new primaryFit(f11,6);

    tempGau->basicCurve(11);

    // creating the Gaussian struct
    fitContainer run11GaussCont;
    
    run11GaussCont.canvNames = "Gauss";
    run11GaussCont.runNumber = 11;
    run11GaussCont.fixMass=1;
    run11GaussCont.fixWidth=0;

    // I'm gonna try to eliminate these strings soon, make their changing more automated
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPtPol1,0); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPtPol2,1); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPtExp,2); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run11GaussCont,run11sGauMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPtPol1,0);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPtPol2,1);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPtExp,2);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run11GaussCont,run11parBoundsGauMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run11GaussCont,GaussFuncs);
    
    fillNumParamsContainerStrucLayer(&run11GaussCont,GaussNumParams);

    tempGau->Analysis(run11GaussCont);

    delete tempGau;

     

    // RUN11 Breit Wigner analysis

    primaryFit *run11BreitWig = new primaryFit(f11,6);
    
    fitContainer run11BreitWigCont;

    run11BreitWigCont.canvNames = "BreitWigner";
    run11BreitWigCont.runNumber = 11;
    run11BreitWigCont.functionCurve = curveBreitW;
    run11BreitWigCont.fixMass = 0;
    run11BreitWigCont.fixWidth = 0;

    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPtPol1,0); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPtPol2,1); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPtExp,2); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run11BreitWigCont,run11sBreitWigMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPtPol1,0);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPtPol2,1);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPtExp,2);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run11BreitWigCont,run11parBoundsBreitWigMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run11BreitWigCont,BreitWigFuncs);
    fillNumParamsContainerStrucLayer(&run11BreitWigCont,BreitWigNumParams);

    run11BreitWig->Analysis(run11BreitWigCont);

    delete run11BreitWig;

    // RUN11 BREIT SOD analysis
   
    primaryFit *run11BreitSod = new primaryFit(f11,6);
    fitContainer run11BreitSodCont;
    
    run11BreitSodCont.canvNames = "BreitSod";
    run11BreitSodCont.runNumber = 11;
    run11BreitSodCont.functionCurve = curveBreitSod;
    run11BreitSodCont.fixMass = 0;
    run11BreitSodCont.fixWidth = 0;

    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPtPol1,0); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPtPol2,1); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPtExp,2); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run11BreitSodCont,run11sBreitSodMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPtPol1,0);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPtPol2,1);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPtExp,2);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run11BreitSodCont,run11parBoundsBreitSodMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run11BreitSodCont,BreitSodFuncs);
    fillNumParamsContainerStrucLayer(&run11BreitSodCont,BreitSodNumParams);

    run11BreitSod->Analysis(run11BreitSodCont);

    delete run11BreitSod;

    // RUN11 MOD SOD analysis
    primaryFit *run11ModSod = new primaryFit(f11,6);
    fitContainer run11ModSodCont;
    
    run11ModSodCont.canvNames = "ModSod";
    run11ModSodCont.runNumber = 11;
    run11ModSodCont.functionCurve = curveModSod;
    run11ModSodCont.fixMass = 0;
    run11ModSodCont.fixWidth = 0;

    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPtPol1,0); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPtPol2,1); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPtExp,2); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run11ModSodCont,run11sModSodMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPtPol1,0);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPtPol2,1);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPtExp,2);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run11ModSodCont,run11parBoundsModSodMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run11ModSodCont,ModSodFuncs);
    fillNumParamsContainerStrucLayer(&run11ModSodCont,ModSodNumParams);

    run11ModSod->Analysis(run11ModSodCont);

    delete run11ModSod;

    // RUN11 ROSS STOD
    primaryFit *run11RossStod = new primaryFit(f11,6);
    fitContainer run11RossStodCont;
    
    run11RossStodCont.canvNames = "RossStod";
    run11RossStodCont.runNumber = 11;
    run11RossStodCont.functionCurve = curveRossStod;
    run11RossStodCont.fixMass = 0;
    run11RossStodCont.fixWidth = 0;

    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPtPol1,0); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPtPol2,1); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPtExp,2); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run11RossStodCont,run11sRossStodMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPtPol1,0);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPtPol2,1);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPtExp,2);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run11RossStodCont,run11parBoundsRossStodMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run11RossStodCont,RossStodFuncs);
    fillNumParamsContainerStrucLayer(&run11RossStodCont,RossStodNumParams);

    run11RossStod->Analysis(run11RossStodCont);

    delete run11RossStod;

    // RUN16 

    // Run 16 Gauss

    primaryFit *run16Gau = new primaryFit(f16_primary,6);
    fitContainer run16GauCont;
    
    run16GauCont.canvNames = "Gauss";
    run16GauCont.runNumber = 16;
  //run16GauCont.functionCurve = curveRossStod; // don't need this for gauss
    run16GauCont.fixMass = 1;
    run16GauCont.fixWidth = 0;

    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPtPol1,0); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPtPol2,1); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPtExp,2); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run16GauCont,run16sGauMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPtPol1,0);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPtPol2,1);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPtExp,2);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run16GauCont,run16parBoundsGauMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run16GauCont,GaussFuncs);
    fillNumParamsContainerStrucLayer(&run16GauCont,GaussNumParams);

    run16Gau->Analysis(run16GauCont);

    delete run16Gau;

    // run 16 Breit Wigner noTOF fit
    primaryFit *run16BreitW = new primaryFit(f16_primary,6);
    fitContainer run16BreitWigCont;
    
    run16BreitWigCont.canvNames = "BreitWig";
    run16BreitWigCont.runNumber = 16;
    run16BreitWigCont.functionCurve = curveBreitW;
    run16BreitWigCont.fixMass = 0;
    run16BreitWigCont.fixWidth = 0;

    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPtPol1,0); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPtPol2,1); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPtExp,2); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run16BreitWigCont,run16sBreitWigMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPtPol1,0);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPtPol2,1);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPtExp,2);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run16BreitWigCont,run16parBoundsBreitWigMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run16BreitWigCont,BreitWigFuncs);
    fillNumParamsContainerStrucLayer(&run16BreitWigCont,BreitWigNumParams);

    run16BreitW->Analysis(run16BreitWigCont); // maybe put the constructor fill in the analysis bit 

    //run16BreitW->f0pTtests();

    delete run16BreitW;

    // run 10 analysis
    primaryFit *run10firstFit = new primaryFit(f10_first,6); 
    fitContainer run10firstCont;

    run10firstCont.canvNames = "Gauss";
    run10firstCont.runNumber = 10;
//    run11testCont.functionCurve = curveBreitW;
    run10firstCont.fixMass = 1;
    run10firstCont.fixWidth = 0;

    fillStContainerStrucLayer(&run10firstCont,run10firstStrPtPol1,0); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPtPol2,1); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPtExp,2); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiExplpt,5); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run10firstCont,run10firstStrPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPtPol1,0);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPtPol2,1);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPtExp,2);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiExplpt,5);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run10firstCont,run10firstparBoundsPhiExphpt,8);

    fillFuncContainerStrucLayer(&run10firstCont,GaussFuncs);
    
    fillNumParamsContainerStrucLayer(&run10firstCont,GaussNumParams);

    run10firstFit->Analysis(run10firstCont);

    run10firstFit->basicCurve(10);

    delete run10firstFit;

    primaryFit *run14_primary = new primaryFit(f14_primary,6); 
    run14_primary->basicCurve(14);

    fitContainer run14BreitWigCont;

    run14BreitWigCont.canvNames = "BreitWigner";
    run14BreitWigCont.runNumber = 14;
    run14BreitWigCont.functionCurve = curveBreitW;
    run14BreitWigCont.fixMass = 1;
    run14BreitWigCont.fixWidth = 0;

    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPtPol1,0); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPtPol2,1); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPtExp,2); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run14BreitWigCont,run14sBreitWigMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPtPol1,0);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPtPol2,1);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPtExp,2);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run14BreitWigCont,run14parBoundsBreitWigMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run14BreitWigCont,BreitWigFuncs);
    fillNumParamsContainerStrucLayer(&run14BreitWigCont,BreitWigNumParams);

    run14_primary->Analysis(run14BreitWigCont);

    run14_primary->SinglePartDistrs();

    run14_primary->f0pTtests();

    delete run14_primary;

    primaryFit *run12_half = new primaryFit(f12_half,6); 
    run12_half->basicCurve(12);

    fitContainer run12half_BreitWigCont;

    run12half_BreitWigCont.canvNames = "BreitWigner";
    run12half_BreitWigCont.runNumber = 12;
    run12half_BreitWigCont.functionCurve = curveBreitW;
    run12half_BreitWigCont.fixMass = 1;
    run12half_BreitWigCont.fixWidth = 0;

    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPtPol1,0); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPtPol2,1); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPtExp,2); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiPol1lpt,3); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiPol2lpt,4); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiExplpt,5); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiPol1hpt,6); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiPol2hpt,7); 
    fillStContainerStrucLayer(&run12half_BreitWigCont,run12half_sBreitWigMassPhiExphpt,8); 
   
    // filling the parameters
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPtPol1,0);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPtPol2,1);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPtExp,2);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiPol1lpt,3);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiPol2lpt,4);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiExplpt,5);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiPol1hpt,6);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiPol2hpt,7);
    fillParamContainerStrucLayer(&run12half_BreitWigCont,run12half_parBoundsBreitWigMassPhiExphpt,8);

    fillFuncContainerStrucLayer(&run12half_BreitWigCont,BreitWigFuncs);
    fillNumParamsContainerStrucLayer(&run12half_BreitWigCont,BreitWigNumParams);

    run12_half->Analysis(run12half_BreitWigCont);

    //run12_primary->SinglePartDistrs();

    //run12_primary->f0pTtests();

    delete run12_half;



    // run 16 momentum correction
    // the momentum correction didn't work
    //primaryFit *run16momCorrPart = new primaryFit(f16_momCorrPart,6); 
    //fitContainer run16momCorrPartCont;

    //run16momCorrPartCont.canvNames = "Gauss";
    //run16momCorrPartCont.runNumber = 1602;
//  //  run11testCont.functionCurve = curveBreitW;
    //run16momCorrPartCont.fixMass = 1;
    //run16momCorrPartCont.fixWidth = 0;

    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPtPol1,0); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPtPol2,1); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPtExp,2); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiPol1lpt,3); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiPol2lpt,4); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiExplpt,5); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiPol1hpt,6); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiPol2hpt,7); 
    //fillStContainerStrucLayer(&run16momCorrPartCont,run16momCorrStrPhiExphpt,8); 
   
    //// filling the parameters
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPtPol1,0);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPtPol2,1);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPtExp,2);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiPol1lpt,3);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiPol2lpt,4);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiExplpt,5);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiPol1hpt,6);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiPol2hpt,7);
    //fillParamContainerStrucLayer(&run16momCorrPartCont,run16momCorrparBoundsPhiExphpt,8);

    //fillFuncContainerStrucLayer(&run16momCorrPartCont,GaussFuncs);
    //
    //fillNumParamsContainerStrucLayer(&run16momCorrPartCont,GaussNumParams);

    //run16momCorrPart->Analysis(run16momCorrPartCont);

    //run16momCorrPart->SinglePartDistrs();
    //run16momCorrPart->f0pTtests();
    //run16momCorrPart->correctionFit();

    // here are some auxillary functions not connected to primary.cpp
    plotDcaComp();
    MCdcaComp();
    MCpTeffic();

    cout << "END MAIN" << endl;
    
    return;
}//end of mainf0macro()
