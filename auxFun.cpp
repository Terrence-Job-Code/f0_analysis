#include "rootheader.h"
#include "fitvalues.h"
// for the same reason 

#ifndef AUXFUN_CPP_
#define AUXFUN_CPP_

//typedef double (*FitFuncPtr) (double *,double *);

const int extraStLength=20;
const int maxNumBins=6;
const int maxNumParam=7;
const int totalDistrs=9; // mass bckgrnds * distribution types
const int numFunParts=3;
const int numMassDistr=3;
const int numBckGrnds=3;
const int axisBounds=2;

// want to move to vectors if root lets me
struct fitContainer
{
    int parBoundsLength;
    int numBins[totalDistrs];

    // names for all of the fits
    TString inserts[totalDistrs][maxNumBins][extraStLength];
    
    double parameterBounds[totalDistrs][maxNumBins][2*maxNumParam]; // I want to delete this, almost updated the rest of my code
   
    // fitting functions 
    FitFuncPtr functionPieces[numMassDistr][numBckGrnds];
   
    FitFuncPtr functionCurve;

    // number of fitting parameters in each function 
    int numberOfParams[numMassDistr][numBckGrnds];
    
    double parameterValues[numMassDistr*numBckGrnds][maxNumBins][2*maxNumParam];
    
    vector < vector < vector<double> > > ParamVals;
    vector < vector<double> > FunctVals;
    vector < vector<double> > BGrndVals; 

    // the y axis range // I don't think this needs to be used now
    //double yAxisVals[numMassDistr*numBckGrnds][maxNumBins][axisBounds];

//    vector < vector<TString> > v2NameCalls;

    // a name to help keep canvases and saves from being deleted
    TString canvNames;

    int runNumber; // can I make this a string?
    int fixMass;
    int fixWidth;  
};

// may adapt later so that all previous 2d vectors become a single 3d vector
//
// I can also probably condense this to one function with different capabilities
// based on what it's filled with

void fillStContainerStrucLayer(struct fitContainer *thePass,vector < vector<TString> > &desiredWords,/* int *length,*/int pos)
{
    //if(length[0]<desiredWords.size()) cout << "array size is too small" << endl;

    for(int i=0;i<desiredWords.size();i++)
    {
        
        //if(length[1]<desiredWords[i].size()) cout << "ARRAY IS TOO SMALL" << endl;

        for(int j=0;j<desiredWords[i].size();j++)
        {
            thePass->inserts[pos][i][j]=desiredWords[i][j];
        }
    }

}

void fillParamContainerStrucLayer(struct fitContainer *thePass,vector < vector<double> > &params,/* int *length,*/int pos)
{
    //if(length[0]<params.size()) cout << "array size is too small" << endl;

    for(int i=0;i<params.size();i++)
    {
        
        //if(length[1]<params[i].size()) cout << "ARRAY IS TOO SMALL" << endl;

        for(int j=0;j<params[i].size();j++)
        {
            thePass->parameterBounds[pos][i][j]=params[i][j]; 
        }
    }

    return;
}

// if I have too many issues with the address, just manually enter the shit
void fillFuncContainerStrucLayer(struct fitContainer *thePass, vector < vector<FitFuncPtr> > &funcs)
{
    if(3<funcs.size()) cout << "array size is too small" << endl;

    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            thePass->functionPieces[i][j]=funcs[i][j];
        }
    }

    return;
}

void fillNumParamsContainerStrucLayer(struct fitContainer *thePass, vector < vector<int> > &numOfParams)
{
    for(int i=0;i<3;i++)
    {
        for(int j=0;j<3;j++)
        {
            thePass->numberOfParams[i][j]=numOfParams[i][j];
        }
    }  

    return;
}
//void fillAxisContainerStrucLayer(struct fitContainer *thePass, vector < vector<double> > &axis,int pos)
//{
//    for(int i=0;i<axis.size();i++)
//    {
//        for(int j=0;j<axis[i].size();j++)
//        {
//            thePass->yAxisVals[pos][i][j]=axis[i][j];
//        }
//    }
//}

// I'm tired of putting stuff in primary.cpp, it's becoming bloated
// I'll put the new plotting code pieces here
// I might even move over much of the primary.cpp code here
void plotDcaComp()
{
    TFile *fprim16 = new TFile("DATA/Run16_primary/run16data200_primary.root");
    //TFile *fglob16 = new TFile("DATA/Run16_global/run16data200_global.root");
    TFile *fglob16 = new TFile("DATA/Run16_full_correct/run16data200.root"); // remaking this data right now to be sure
    TFile *fprim14 = new TFile("DATA/Run14_nohft_primary/run14data200_primary.root");
    TFile *fglob14 = new TFile("DATA/Run14_nohft_global/run14data200_global.root");

    // getting the normalization numbers

    TH1F *primCent14;
    TH1F *globCent14;
    TH1F *primCent16;
    TH1F *globCent16;

    double numEventsPrim14;
    double numEventsGlob14;
    double numEventsPrim16;
    double numEventsGlob16;

    primCent14 = (TH1F*)fprim14->Get("centralityw");
    globCent14 = (TH1F*)fglob14->Get("centralityw");
    primCent16 = (TH1F*)fprim16->Get("centralityw");
    globCent16 = (TH1F*)fglob16->Get("centralityw");

    numEventsPrim14 = primCent14->GetBinContent(1);
    numEventsGlob14 = globCent14->GetBinContent(1);
    numEventsPrim16 = primCent16->GetBinContent(1);
    numEventsGlob16 = globCent16->GetBinContent(1);

    // finished the normalization numbers

    TH1F *ppDca14;
    TH1F *pgDca14;
    TH1F *ppDca16;
    TH1F *pgDca16;

    ppDca14 = (TH1F*)fprim14->Get("pgDca");
    pgDca14 = (TH1F*)fglob14->Get("pgDca");
    ppDca16 = (TH1F*)fprim16->Get("pgDca");
    pgDca16 = (TH1F*)fglob16->Get("pgDca");
    
    // normalize the dca plots

    double rawNum = 1;
    double scalNum = 1;

    for(int i = 1; i<ppDca14->GetNbinsX()+1;i++)
    {
        rawNum = ppDca14->GetBinContent(i);
        scalNum = rawNum/numEventsPrim14;
        ppDca14->SetBinContent(i,scalNum);

        rawNum = pgDca14->GetBinContent(i);
        scalNum = rawNum/numEventsGlob14;
        pgDca14->SetBinContent(i,scalNum);

        rawNum = ppDca16->GetBinContent(i);
        scalNum = rawNum/numEventsPrim16;
        ppDca16->SetBinContent(i,scalNum);

        rawNum = pgDca16->GetBinContent(i);
        scalNum = rawNum/numEventsGlob16;
        pgDca16->SetBinContent(i,scalNum);

    }


    TCanvas *c01 = new TCanvas();

    pgDca14->SetLineColor(2);
    pgDca14->SetMarkerColor(2);

    ppDca16->SetLineColor(3);
    ppDca16->SetMarkerColor(3);

    pgDca16->SetLineColor(4);
    pgDca16->SetMarkerColor(4);

    ppDca16->Draw();
    ppDca14->Draw("same");
    pgDca14->Draw("same");
    pgDca16->Draw("same");

    TLegend *leg = new TLegend(0.6,0.25,1.0,0.45);
    leg->AddEntry(ppDca14,"run 14 primary");
    leg->AddEntry(pgDca14,"run 14 global");
    leg->AddEntry(ppDca16,"run 16 primary");
    leg->AddEntry(pgDca16,"run 16 global");
    leg->Draw();

    c01->SaveAs("mainPlots/testPlots/DcaCompare.pdf");

    return;
}

// this code became a monster :(
void MCdcaComp()
{
    // compares monte carlo dca and data dca
    TFile *fpr16 = new TFile("DATA/Run16_primary/run16data200_primary.root");
    TFile *fpr14 = new TFile("DATA/Run14_nohft_primary/run14data200_primary.root");
    TFile *fgl16 = new TFile("DATA/Run16_full_correct/run16data200.root");
    TFile *fgl14 = new TFile("DATA/Run14_nohft_global/run14data200_global.root");
    //TFile *fpr11 = new TFile("DATA/OLDDATA/run11Data.root");
    TFile *fpr11 = new TFile("DATA/OLDDATA/run11test.root"); // using this because my full run11 data does not have ppDca

    //TFile *fmc14 = new TFile("pionEmbedding/run14/trial7/observe.root");
    //TFile *fmc14 = new TFile("pionEmbedding/run14/trial9/fullout.root");
    //TFile *fmc14 = new TFile("pionEmbedding/run14/trial10/partOut.root");
    TFile *fmc14 = new TFile("pionEmbedding/extraTesting01/run14embed.root");
    TFile *f14dat= new TFile("pionEmbedding/extraTesting01/run14data200_primaryTest01.root");

    // DCA plots

    TH1D *hMCPrDca;
    TH1D *hMCGlDca;
    TH1F *h14PrDca;
    TH1F *h14GlDca;
    TH1F *h16PrDca;
    TH1F *h16GlDca;
    TH1F *h11PrDca;
    TH2D *hMCdcavsPt;
    TH2D *hrun11dcavsPt;
    TH2D *hrun14dcavsPt;

    TH1D *hprojmc01;
    TH1D *hprojmc02;
    TH1D *hprojmc03;
    TH1D *hprojmc04;
    TH1D *hprojmc05;
    TH1D *hprojmc06;

    // can't use the run 11 comparisons yet 
    TH1D *hprojr11dat11;
    TH1D *hprojr11dat02;
    TH1D *hprojr11dat03;
    TH1D *hprojr11dat04;
    TH1D *hprojr11dat05;
    TH1D *hprojr11dat06;
 
    TH1D *hprojr14dat01;
    TH1D *hprojr14dat02;
    TH1D *hprojr14dat03;
    TH1D *hprojr14dat04;
    TH1D *hprojr14dat05;
    TH1D *hprojr14dat06;
 
    hMCPrDca = (TH1D*)fmc14->Get("PiPluPrimaryDCA");// this one doesn't seem important because global dca is used in the data
    hMCGlDca = (TH1D*)fmc14->Get("PiPluGlobalDCA");
    h14PrDca = (TH1F*)fpr14->Get("pgDca"); 
    h16PrDca = (TH1F*)fpr16->Get("pgDca");
    h11PrDca = (TH1F*)fpr11->Get("pgDca");
    hMCdcavsPt = (TH2D*)fmc14->Get("PipluDcavsPt");// x = pt; y = dca , for pt, 4 bins for each 1 GeV, starting at -1
    hrun14dcavsPt = (TH2D*)f14dat->Get("dcavsPt");
    
    // new bins coming up
    hprojmc01 = hMCdcavsPt->ProjectionY("pt01",1,80);//0-2 GeV
    hprojmc02 = hMCdcavsPt->ProjectionY("pt02",81,200);// 2-5 GeV
    hprojmc03 = hMCdcavsPt->ProjectionY("pt03",31,32);// 0.5 - 0.6 GeV
    hprojmc04 = hMCdcavsPt->ProjectionY("pt04",41,44);// 1.0 - 1.2 GeV
    hprojmc05 = hMCdcavsPt->ProjectionY("pt05",53,56); // 1.6 - 1.8 GeV
    hprojmc06 = hMCdcavsPt->ProjectionY("pt06",11,12);// 1.5 - 2.0 GeV

    hprojr14dat01 = hrun14dcavsPt->ProjectionY("r14pt01",10,12); // 0.5-0.6 GeV
    hprojr14dat02 = hrun14dcavsPt->ProjectionY("r14pt02",10,12); // 0.5-0.6 GeV
    hprojr14dat03 = hrun14dcavsPt->ProjectionY("r14pt03",11,12); // 0.5-0.6 GeV
    hprojr14dat04 = hrun14dcavsPt->ProjectionY("r14pt04",21,24); // 1.0-1.2 GeV
    hprojr14dat05 = hrun14dcavsPt->ProjectionY("r14pt05",33,36); // 1.6-1.8 GeV
    hprojr14dat06 = hrun14dcavsPt->ProjectionY("r14pt06",10,12); // 0.5-0.6 GeV

    cout << hMCGlDca->Integral() << endl;
    cout << hprojmc01->Integral() << endl;

    // the scale function applies the action to all bins

    hMCGlDca->Scale(1./hMCGlDca->GetBinWidth(1));
    hMCGlDca->Scale(1./hMCGlDca->Integral());

    h11PrDca->Scale(1./h11PrDca->GetBinWidth(1));
    h11PrDca->Scale(1./h11PrDca->Integral());

    h14PrDca->Scale(1./h14PrDca->GetBinWidth(1));
    h14PrDca->Scale(1./h14PrDca->Integral());

    hprojmc01->Scale(1./hprojmc01->GetBinWidth(1));
    hprojmc01->Scale(1./hprojmc01->Integral());

    hprojmc02->Scale(1./hprojmc02->GetBinWidth(1));
    hprojmc02->Scale(1./hprojmc02->Integral());

    hprojmc03->Scale(1./hprojmc03->GetBinWidth(1));
    hprojmc03->Scale(1./hprojmc03->Integral());

    hprojmc04->Scale(1./hprojmc04->GetBinWidth(1));
    hprojmc04->Scale(1./hprojmc04->Integral());

    hprojmc05->Scale(1./hprojmc05->GetBinWidth(1));
    hprojmc05->Scale(1./hprojmc05->Integral());

    hprojr14dat01->Scale(1./hprojr14dat01->GetBinWidth(1));
    hprojr14dat01->Scale(1./hprojr14dat01->Integral());

    hprojr14dat02->Scale(1./hprojr14dat02->GetBinWidth(1));
    hprojr14dat02->Scale(1./hprojr14dat02->Integral());

    hprojr14dat03->Scale(1./hprojr14dat03->GetBinWidth(1));
    hprojr14dat03->Scale(1./hprojr14dat03->Integral());

    hprojr14dat04->Scale(1./hprojr14dat04->GetBinWidth(1));
    hprojr14dat04->Scale(1./hprojr14dat04->Integral());

    hprojr14dat05->Scale(1./hprojr14dat05->GetBinWidth(1));
    hprojr14dat05->Scale(1./hprojr14dat05->Integral());

    TCanvas *c01 = new TCanvas();
    //c01->SetLogy();

    hMCGlDca->GetXaxis()->SetTitle("DCA cm");
    hMCGlDca->GetYaxis()->SetTitle("entries/total entries");
    hMCGlDca->GetXaxis()->SetRangeUser(0,2);

    h11PrDca->GetXaxis()->SetTitle("DCA cm");
    h11PrDca->GetYaxis()->SetTitle("entries/total entries");
    h11PrDca->GetXaxis()->SetRangeUser(0,2);
    
    h14PrDca->GetXaxis()->SetRangeUser(0,2);

    hMCGlDca->SetLineColor(1);
    h14PrDca->SetLineColor(2);
    h16PrDca->SetLineColor(3);  
    h11PrDca->SetLineColor(4);
    hprojmc01->SetLineColor(6); 

    hMCGlDca->Draw("hist");
    h11PrDca->Draw("same hist");
    h14PrDca->Draw("same hist");
    hprojmc01->Draw("same hist");
    //hMCGlDca->Draw("same hist");
    //h16PrDca->Draw("same");
    //h11PrDca->Draw("same");

    TLegend *leg = new TLegend(0.6,0.25,1.0,0.45);
    //leg->AddEntry(hMCPrDca,"run 14 monte carlo primary dca");
    leg->AddEntry(hMCGlDca,"run 14 monte carlo global dca");
    leg->AddEntry(h14PrDca,"run 14 data dca");
    //leg->AddEntry(h16PrDca,"run 16 data dca");
    leg->AddEntry(h11PrDca,"run 11 data dca");
    leg->AddEntry(hprojmc01,"run 14 embedding, 0-2 GeV");
    leg->Draw();

    c01->SaveAs("mainPlots/testPlots/DcaMonteCarloCompare.pdf");
    //c01->Close();

    TCanvas *c02 = new TCanvas();

    hMCGlDca->GetXaxis()->SetTitle("DCA cm");
    hMCGlDca->GetYaxis()->SetTitle("entries/total entries");
    hMCGlDca->GetXaxis()->SetRangeUser(0,2);

    h11PrDca->GetXaxis()->SetTitle("DCA cm");
    h11PrDca->GetYaxis()->SetTitle("entries/total entries");
    h11PrDca->GetXaxis()->SetRangeUser(0,2);
    
    h14PrDca->GetXaxis()->SetRangeUser(0,2);

    hMCGlDca->SetLineColor(1);
    h14PrDca->SetLineColor(2);
    h16PrDca->SetLineColor(3);  
    h11PrDca->SetLineColor(4);
    hprojmc02->SetLineColor(6); 

    hMCGlDca->Draw("hist");
    h11PrDca->Draw("same hist");
    h14PrDca->Draw("same hist");
    hprojmc02->Draw("same hist");
    //hMCGlDca->Draw("same hist");
    //h16PrDca->Draw("same");
    //h11PrDca->Draw("same");

    TLegend *leg02 = new TLegend(0.6,0.25,1.0,0.45);
    //leg->AddEntry(hMCPrDca,"run 14 monte carlo primary dca");
    leg02->AddEntry(hMCGlDca,"run 14 monte carlo global dca");
    leg02->AddEntry(h14PrDca,"run 14 data dca");
    //leg->AddEntry(h16PrDca,"run 16 data dca");
    leg02->AddEntry(h11PrDca,"run 11 data dca");
    leg02->AddEntry(hprojmc02,"run 14 embedding, 2-5 GeV");
    leg02->Draw();

    c02->SaveAs("mainPlots/embedPlots/dcaComp02.pdf");

    // all histograms will be binned projections after this point
    TCanvas *c03 = new TCanvas();

    hMCGlDca->GetXaxis()->SetTitle("DCA cm");
    hMCGlDca->GetYaxis()->SetTitle("entries/total entries");
    hMCGlDca->GetXaxis()->SetRangeUser(0,2);

    h11PrDca->GetXaxis()->SetTitle("DCA cm");
    h11PrDca->GetYaxis()->SetTitle("entries/total entries");
    h11PrDca->GetXaxis()->SetRangeUser(0,2);
    
    hprojr14dat03->GetXaxis()->SetRangeUser(0,2);

    hMCGlDca->SetLineColor(1);
    hprojr14dat03->SetLineColor(2);
    h16PrDca->SetLineColor(3);  
    h11PrDca->SetLineColor(4);
    hprojmc03->SetLineColor(6); 

    hMCGlDca->Draw("hist");
    h11PrDca->Draw("same hist");
    hprojr14dat03->Draw("same hist");
    hprojmc03->Draw("same hist");
    //hMCGlDca->Draw("same hist");
    //h16PrDca->Draw("same");
    //h11PrDca->Draw("same");

    TLegend *leg03 = new TLegend(0.6,0.25,1.0,0.45);
    leg03->AddEntry(hMCGlDca,"run 14 monte carlo global dca");
    leg03->AddEntry(hprojr14dat03,"run 14 data dca, 0.5-0.6 GeV");
    leg03->AddEntry(h11PrDca,"run 11 data dca");
    leg03->AddEntry(hprojmc03,"run 14 embedding, 0.5-0.6 GeV");
    leg03->Draw();

    c03->SaveAs("mainPlots/embedPlots/dcaComp03.pdf");

    TCanvas *c04 = new TCanvas();

    hMCGlDca->GetXaxis()->SetTitle("DCA cm");
    hMCGlDca->GetYaxis()->SetTitle("entries/total entries");
    hMCGlDca->GetXaxis()->SetRangeUser(0,2);

    h11PrDca->GetXaxis()->SetTitle("DCA cm");
    h11PrDca->GetYaxis()->SetTitle("entries/total entries");
    h11PrDca->GetXaxis()->SetRangeUser(0,2);
    
    hprojr14dat03->GetXaxis()->SetRangeUser(0,2);

    hMCGlDca->SetLineColor(1);
    hprojr14dat04->SetLineColor(2);
    h16PrDca->SetLineColor(3);  
    h11PrDca->SetLineColor(4);
    hprojmc04->SetLineColor(6); 

    hMCGlDca->Draw("hist");
    h11PrDca->Draw("same hist");
    hprojr14dat04->Draw("same hist");
    hprojmc04->Draw("same hist");
    //hMCGlDca->Draw("same hist");
    //h16PrDca->Draw("same");
    //h11PrDca->Draw("same");

    TLegend *leg04 = new TLegend(0.6,0.25,1.0,0.45);
    leg04->AddEntry(hMCGlDca,"run 14 monte carlo global dca");
    leg04->AddEntry(hprojr14dat04,"run 14 data dca, 1.0-1.2 GeV");
    leg04->AddEntry(h11PrDca,"run 11 data dca");
    leg04->AddEntry(hprojmc04,"run 14 embedding, 1.0-1.2 GeV");
    leg04->Draw();

    c04->SaveAs("mainPlots/embedPlots/dcaComp04.pdf");

    TCanvas *c05 = new TCanvas();

    hMCGlDca->GetXaxis()->SetTitle("DCA cm");
    hMCGlDca->GetYaxis()->SetTitle("entries/total entries");
    hMCGlDca->GetXaxis()->SetRangeUser(0,2);

    h11PrDca->GetXaxis()->SetTitle("DCA cm");
    h11PrDca->GetYaxis()->SetTitle("entries/total entries");
    h11PrDca->GetXaxis()->SetRangeUser(0,2);
    
    hprojr14dat05->GetXaxis()->SetRangeUser(0,2);

    hMCGlDca->SetLineColor(1);
    hprojr14dat05->SetLineColor(2);
    h16PrDca->SetLineColor(3);  
    h11PrDca->SetLineColor(4);
    hprojmc05->SetLineColor(6); 

    hMCGlDca->Draw("hist");
    h11PrDca->Draw("same hist");
    hprojr14dat05->Draw("same hist");
    hprojmc05->Draw("same hist");
    //hMCGlDca->Draw("same hist");
    //h16PrDca->Draw("same");
    //h11PrDca->Draw("same");

    TLegend *leg05 = new TLegend(0.6,0.25,1.0,0.45);
    leg05->AddEntry(hMCGlDca,"run 14 monte carlo global dca");
    leg05->AddEntry(hprojr14dat05,"run 14 data dca, 1.6-1.8 GeV");
    leg05->AddEntry(h11PrDca,"run 11 data dca");
    leg05->AddEntry(hprojmc05,"run 14 embedding, 1.6-1.8 GeV");
    leg05->Draw();

    c05->SaveAs("mainPlots/embedPlots/dcaComp05.pdf");

    return;
}

// for now I'll just plot the efficiencies
void MCpTeffic()
{
    TFile *fpr16 = new TFile("DATA/Run16_primary/run16data200_primary.root");
    TFile *fpr14 = new TFile("DATA/Run14_nohft_primary/run14data200_primary.root");
    TFile *fgl16 = new TFile("DATA/Run16_full_correct/run16data200.root");
    TFile *fgl14 = new TFile("DATA/Run14_nohft_global/run14data200_global.root");

    //TFile *fmc14 = new TFile("pionEmbedding/run14/trial7/observe.root");
    //TFile *fmc14 = new TFile("pionEmbedding/run14/trial9/fullout.root");
    TFile *fmc14 = new TFile("pionEmbedding/run14/trial10/partOut.root");

    // number of events
    TH1D *McEvents;
    TH1D *run14Cent;
    TH1D *run16Cent;

    double numMcEv;
    double num14Ev;
    double num16Ev;

    McEvents = (TH1D*)fmc14->Get("DeltaVertexX");
    run14Cent = (TH1D*)fpr14->Get("centralityw");
    run16Cent = (TH1D*)fpr16->Get("centralityw");

    numMcEv = McEvents->Integral();
    num14Ev = run14Cent->GetBinContent(1); 
    num16Ev = run16Cent->GetBinContent(1); 

    // pt distributions
    TH1D *hMCPt;
    TH1D *hRcPrPt;
    TH1D *hRcGlPt;
    TH1D *hpr16Pt;
    TH1D *hpr14Pt;

    hMCPt = (TH1D*)fmc14->Get("PiPluMcPt");     
    hRcPrPt = (TH1D*)fmc14->Get("PiPluPrimaryPt");
    hRcGlPt = (TH1D*)fmc14->Get("PiPluGlobalPt");
    hpr16Pt = (TH1D*)fpr16->Get("ppt");
    hpr14Pt = (TH1D*)fpr14->Get("ppt");

    // can't normalize by unity
    hMCPt->Scale(1./numMcEv);
    hRcPrPt->Scale(1./numMcEv);
    hRcGlPt->Scale(1./numMcEv);
    hpr16Pt->Scale(1./num16Ev);
    hpr14Pt->Scale(1./num14Ev);

    // I should probably include the data pT distributions as well for comparison
    // will also normalize to unity
    TCanvas *c00 = new TCanvas();
    c00->SetLogy();
    hMCPt->SetLineColor(3);
    hRcGlPt->SetLineColor(2);
    hpr16Pt->SetLineColor(4);
    hpr14Pt->SetLineColor(6);
   
    hpr16Pt->SetTitle("#pi + pT distributions");
    hpr16Pt->GetXaxis()->SetTitle("pT GeV");
    hpr16Pt->GetYaxis()->SetTitle("entries/events");

    hpr16Pt->Draw("hist");
    hMCPt->Draw("hist same");
    hRcPrPt->Draw("hist same");
    hRcGlPt->Draw("hist same");
    hpr16Pt->Draw("hist same");
    hpr14Pt->Draw("hist same");

    TLegend *leg01 = new TLegend(0.7,0.1,0.9,0.3);
    leg01->AddEntry(hMCPt,"run 14 mc");
    leg01->AddEntry(hRcPrPt,"run 14 reconstructed pr");
    leg01->AddEntry(hRcGlPt,"run 14 reconstructed gl");
    leg01->AddEntry(hpr16Pt,"run 16 primary data");
    leg01->AddEntry(hpr14Pt,"run 14 primary data");
    leg01->Draw();


    c00->SaveAs("mainPlots/testPlots/montecarlopTDists.pdf");

    double ratioPrMc = 0;
    double ratioGlMc = 0;

    for(int i = 0; i < hMCPt->GetNbinsX() + 1; i++)
    {
        if(hMCPt->GetBinContent(i)!=0)
        { 
            ratioPrMc = hRcPrPt->GetBinContent(i)/hMCPt->GetBinContent(i);
            ratioGlMc = hRcGlPt->GetBinContent(i)/hMCPt->GetBinContent(i);

            hRcPrPt->SetBinContent(i,ratioPrMc);   
            hRcGlPt->SetBinContent(i,ratioGlMc);
        }
        else
        {
            cout << i << "\t Mc is zero" << endl;
        }
    }

    TCanvas *c01 = new TCanvas();
    //c01->SetLogy();

    hRcPrPt->SetTitle("Efficiency vs Pt");
    hRcPrPt->GetXaxis()->SetTitle("pT");
    hRcPrPt->GetYaxis()->SetTitle("reconstructed/montecarlo");

    hRcGlPt->SetLineColor(2);
    hRcPrPt->Draw("hist");
    hRcGlPt->Draw("hist same");

    TLegend *leg = new TLegend(0.1,0.5,0.5,0.7);
    leg->AddEntry(hRcPrPt,"run 14 pr/mc");
    leg->AddEntry(hRcGlPt,"run 14 gl/mc");
    leg->Draw();

    c01->SaveAs("mainPlots/testPlots/montecarloefficiency.pdf");
    //c01->Close();

    return;
}

// because these are profiles of 2d histograms, I shouldn't have to normalize them
void nHitsvsPt()
{
    TFile *fdat14 = new TFile("pionEmbedding/extraTesting01/run14data200_primaryTest01.root");
    TFile *fmc14 = new TFile("pionEmbedding/extraTesting01/run14embed.root");

    TH2D *hmcnHits;
    TH2D *hdatnHits;

    hmcnHits = (TH2D*)fmc14->Get("PiplunHitsvsPt");
    hdatnHits = (TH2D*)fdat14->Get("nHitsvsPt");

    TCanvas *c0 = new TCanvas();
    
    hmcnHits->ProfileX()->GetYaxis()->SetTitle("nHits");
    hmcnHits->ProfileX()->GetXaxis()->SetTitle("pt GeV");
    hmcnHits->ProfileX()->Draw();

    hdatnHits->ProfileX()->SetLineColor(2);
    hdatnHits->ProfileX()->SetMarkerColor(2);
    hdatnHits->ProfileX()->Draw("same");

    c0->SaveAs("mainPlots/embedPlots/nHitsvsptpfx.pdf");

    TCanvas *c1 = new TCanvas();

    hmcnHits->ProfileY()->GetXaxis()->SetTitle("nHits");
    hmcnHits->ProfileY()->GetYaxis()->SetTitle("pt GeV");
    hmcnHits->ProfileY()->Draw();

    hdatnHits->ProfileY()->SetLineColor(2);
    hdatnHits->ProfileY()->SetMarkerColor(2);
    hdatnHits->ProfileY()->Draw("same");

    c1->SaveAs("mainPlots/embedPlots/nHitsvsptpfy.pdf");

    // projections

    TH1D *hprojmc01;
    TH1D *hprojmc02;
    TH1D *hprojmc03;
    TH1D *hprojmc04;
    TH1D *hprojmc05;
    TH1D *hprojmc06;

    // can't use the run 11 comparisons yet 
    TH1D *hprojr11dat11;
    TH1D *hprojr11dat02;
    TH1D *hprojr11dat03;
    TH1D *hprojr11dat04;
    TH1D *hprojr11dat05;
    TH1D *hprojr11dat06;
 
    TH1D *hprojr14dat01;
    TH1D *hprojr14dat02;
    TH1D *hprojr14dat03;
    TH1D *hprojr14dat04;
    TH1D *hprojr14dat05;
    TH1D *hprojr14dat06;
 
    hprojmc01 = hmcnHits->ProjectionY("pt01",1,80);//0-2 GeV
    hprojmc02 = hmcnHits->ProjectionY("pt02",81,200);// 2-5 GeV
    hprojmc03 = hmcnHits->ProjectionY("pt03",31,32);// 0.5 - 0.6 GeV
    hprojmc04 = hmcnHits->ProjectionY("pt04",41,44);// 1.0 - 1.2 GeV
    hprojmc05 = hmcnHits->ProjectionY("pt05",53,56); // 1.6 - 1.8 GeV
    hprojmc06 = hmcnHits->ProjectionY("pt06",11,12);// 1.5 - 2.0 GeV

    hprojr14dat01 = hdatnHits->ProjectionY("r14pt01",10,12); // 0.5-0.6 GeV
    hprojr14dat02 = hdatnHits->ProjectionY("r14pt02",10,12); // 0.5-0.6 GeV
    hprojr14dat03 = hdatnHits->ProjectionY("r14pt03",11,12); // 0.5-0.6 GeV
    hprojr14dat04 = hdatnHits->ProjectionY("r14pt04",21,24); // 1.0-1.2 GeV
    hprojr14dat05 = hdatnHits->ProjectionY("r14pt05",33,36); // 1.6-1.8 GeV
    hprojr14dat06 = hdatnHits->ProjectionY("r14pt06",10,12); // 0.5-0.6 GeV

    hprojmc01->Scale(1./hprojmc01->GetBinWidth(1));
    hprojmc01->Scale(1./hprojmc01->Integral());

    hprojmc02->Scale(1./hprojmc02->GetBinWidth(1));
    hprojmc02->Scale(1./hprojmc02->Integral());

    hprojmc03->Scale(1./hprojmc03->GetBinWidth(1));
    hprojmc03->Scale(1./hprojmc03->Integral());

    hprojmc04->Scale(1./hprojmc04->GetBinWidth(1));
    hprojmc04->Scale(1./hprojmc04->Integral());

    hprojmc05->Scale(1./hprojmc05->GetBinWidth(1));
    hprojmc05->Scale(1./hprojmc05->Integral());

    hprojr14dat01->Scale(1./hprojr14dat01->GetBinWidth(1));
    hprojr14dat01->Scale(1./hprojr14dat01->Integral());

    hprojr14dat02->Scale(1./hprojr14dat02->GetBinWidth(1));
    hprojr14dat02->Scale(1./hprojr14dat02->Integral());

    hprojr14dat03->Scale(1./hprojr14dat03->GetBinWidth(1));
    hprojr14dat03->Scale(1./hprojr14dat03->Integral());

    hprojr14dat04->Scale(1./hprojr14dat04->GetBinWidth(1));
    hprojr14dat04->Scale(1./hprojr14dat04->Integral());

    hprojr14dat05->Scale(1./hprojr14dat05->GetBinWidth(1));
    hprojr14dat05->Scale(1./hprojr14dat05->Integral());


    TCanvas *c03 = new TCanvas();

    hprojmc03->GetXaxis()->SetTitle("nHits");
    hprojmc03->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat03->GetXaxis()->SetTitle("nHits");
    hprojr14dat03->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat03->SetLineColor(2);
    hprojmc03->SetLineColor(6); 

    hprojr14dat03->Draw("hist");
    hprojmc03->Draw("same hist");

    TLegend *leg03 = new TLegend(0.6,0.25,1.0,0.45);
    leg03->AddEntry(hprojr14dat03,"run 14 data nHits, 0.5-0.6 GeV");
    leg03->AddEntry(hprojmc03,"run 14 embedding, 0.5-0.6 GeV");
    leg03->Draw();

    c03->SaveAs("mainPlots/embedPlots/nHitsComp03.pdf");

    TCanvas *c04 = new TCanvas();

    hprojmc04->GetXaxis()->SetTitle("nHits");
    hprojmc04->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat04->GetXaxis()->SetTitle("nHits");
    hprojr14dat04->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat04->SetLineColor(2);
    hprojmc04->SetLineColor(6); 

    hprojr14dat04->Draw("hist");
    hprojmc04->Draw("same hist");

    TLegend *leg04 = new TLegend(0.6,0.25,1.0,0.45);
    leg04->AddEntry(hprojr14dat04,"run 14 data nHits, 1.0-1.2 GeV");
    leg04->AddEntry(hprojmc04,"run 14 embedding, 1.0-1.2 GeV");
    leg04->Draw();

    c04->SaveAs("mainPlots/embedPlots/nHitsComp04.pdf");

    TCanvas *c05 = new TCanvas();

    hprojmc05->GetXaxis()->SetTitle("nHits");
    hprojmc05->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat05->GetXaxis()->SetTitle("nHits");
    hprojr14dat05->GetYaxis()->SetTitle("entries/total entries");

    hprojr14dat05->SetLineColor(2);
    hprojmc05->SetLineColor(6); 

    hprojr14dat05->Draw("hist");
    hprojmc05->Draw("same hist");

    TLegend *leg05 = new TLegend(0.6,0.25,1.0,0.45);
    leg05->AddEntry(hprojr14dat05,"run 14 data nHits, 1.6-1.8 GeV");
    leg05->AddEntry(hprojmc05,"run 14 embedding, 1.6-1.8 GeV");
    leg05->Draw();

    c05->SaveAs("mainPlots/embedPlots/nHitsComp05.pdf");



    return;
}

void gdcavsPt()
{
    TFile *fdat14 = new TFile("pionEmbedding/extraTesting01/run14data200_primaryTest01.root");
    TFile *fmc14 = new TFile("pionEmbedding/extraTesting01/run14embed.root");

    TH2D *hmcdca;
    TH2D *hdatdca;

    hmcdca = (TH2D*)fmc14->Get("PipluDcavsPt");
    hdatdca = (TH2D*)fdat14->Get("dcavsPt");

    //PmcdcaX = hmcdca->ProfileX();
    //PdatdcaX = hdatdca->ProfileX();

    TCanvas *c0 = new TCanvas();
   
    //PmcdcaX->Draw();
    //PmcdcaX->GetYaxis()->SetTitle("dca cm");
    //PmcdcaX->GetXaxis()->SetTitle("pt GeV");

    //PdatdcaX->SetMarkerStyle(2);
    //PdatdcaX->SetLineColor(2);
    //PdatdcaX->Draw("same");
    //PdatdcaX->SetMarkerColor(2);

    //PmcdcaX->GetYaxis()->SetTitle("dca cm");
    //PmcdcaX->GetXaxis()->SetTitle("pt GeV");


    //hmcdca->ProfileX()->Draw();
    //hmcdca->ProfileX()->GetYaxis()->SetTitle("dca cm");
    //hmcdca->ProfileX()->GetXaxis()->SetTitle("pt GeV");

    //hdatdca->ProfileX()->SetMarkerStyle(2);
    //hdatdca->ProfileX()->SetLineColor(2);
    //hdatdca->ProfileX()->Draw("same");
    //hdatdca->ProfileX()->SetMarkerColor(2);

    TLegend *leg0 = new TLegend(0.7,0.1,0.9,0.3);
    leg0->AddEntry(hmcdca->ProfileX(),"run 14 embedding");
    leg0->AddEntry(hdatdca->ProfileX(),"run 14 data");
    //leg0->Draw();

    c0->SaveAs("mainPlots/embedPlots/dcavsptpfx.pdf");

    TCanvas *c1 = new TCanvas();

    hmcdca->ProfileY()->Draw();
    hmcdca->ProfileY()->GetYaxis()->SetTitle("pt GeV");
    hmcdca->ProfileY()->GetXaxis()->SetTitle("dca cm");

    hdatdca->ProfileY()->SetMarkerStyle(2);
    hdatdca->ProfileY()->SetLineColor(2);
    hdatdca->ProfileY()->Draw("same");
    hdatdca->ProfileY()->SetMarkerColor(2);

    TLegend *leg1 = new TLegend(0.7,0.1,0.9,0.3);
    leg1->AddEntry(hmcdca->ProfileY(),"run 14 embedding");
    leg1->AddEntry(hdatdca->ProfileY(),"run 14 data");
    //leg1->Draw();

    c1->SaveAs("mainPlots/embedPlots/dcavsptpfy.pdf");


    return;
}

// I want to make a test to see how small an amount of data I can use for run11

#endif
