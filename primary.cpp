// #include "rootheader.h"
#include "primary.h"
#include "fitvalues.h"
#include "fits.h"
//#include "container.h"

// reference
// phys Rev Lett 89 27302 Adler

using namespace std;

// USE "/TEMPORARY" TO FIND TEMPORARY EDITS!!!!

primaryFit::primaryFit(TString fileName,int centrality)
{
    file = new TFile(fileName);
    cent=centrality;
    char ULvar[50];
    char LSvar[50]; 
    for(int i=0;i<numCents;i++)
    {
        sprintf(ULvar,"TMdphiUL0c%i",i+1);
        sprintf(LSvar,"TMdphiLS0c%i",i+1);
        m_ULhists[i] = (TH3D*)file->Get(ULvar);
        m_LShists[i] = (TH3D*)file->Get(LSvar); 
    }
    for(int i=0;i<cent;i++)
    {
        m_ULhists[0]->Add(m_ULhists[i+1],1);
        m_LShists[0]->Add(m_LShists[i+1],1);
    }

}
primaryFit::~primaryFit()
{
    cout << "destroying object" << endl;
    //delete ULhists[numCents];
    //delete LShists[numCents];
    delete file;
}

void primaryFit::fillConstr(fitContainer intermediateBox)
{
    
    // this function only needs to be concerned about initial values

    m_infoBox.canvNames = intermediateBox.canvNames;

    m_infoBox.runNumber = intermediateBox.runNumber;
    m_infoBox.fixMass = intermediateBox.fixMass;
    m_infoBox.fixWidth = intermediateBox.fixWidth;

    m_infoBox.parBoundsLength = intermediateBox.parBoundsLength;

    for(int i = 0; i < totalDistrs; i++)
    {

        m_infoBox.numBins[i] = intermediateBox.numBins[i];

        for(int j = 0; j<maxNumBins; j++)
        {

            for(int k = 0; k<2*maxNumParam; k++)
            {
                m_infoBox.inserts[i][j][k] = intermediateBox.inserts[i][j][k];
            }

        }

    }

    for(int i = 0; i < numMassDistr; i++)
    {
        for(int j = 0; j < numBckGrnds;j++)
        {
            m_infoBox.functionPieces[i][j] = intermediateBox.functionPieces[i][j];
            m_infoBox.numberOfParams[i][j] = intermediateBox.numberOfParams[i][j];
        }
    }

    for(int i = 0; i < totalDistrs; i++ )
    {
        for(int j = 0; j < maxNumBins; j++)
        {
            for(int k = 0; k < (2*maxNumParam);k++ )
            {
                m_infoBox.parameterBounds[i][j][k] = intermediateBox.parameterBounds[i][j][k];
            }
        }
    }

    return;
}

// just for an easy plot call
void primaryFit::basicCurve(int runNum)
{
    // normalize these distributions
    // full x axis range
    TH1D *h1u0[3];
    TH1D *h1l0[3];

    // x axis limited to f0 peak
    TH1D *h1u1[3];
    TH1D *h1l1[3];

    // some names
    char ULname[256];
    char LSname[256];
    char ULgraph[256];
    char LSgraph[256];
    char fullName[256];
    char f0Name[256];
    char fullGraph[256];

    sprintf(ULname,"%iULfull.pdf",runNum);
    sprintf(LSname,"%iLSfull.pdf",runNum);
    sprintf(ULgraph,"run%i UL distribution",runNum);
    sprintf(LSgraph,"run%i LS distribution",runNum);
    sprintf(fullName,"%iFullDistr.pdf",runNum);
    sprintf(f0Name,"%if0PeakDistr.pdf",runNum);
    sprintf(fullGraph,"run%i Mass Distr UL-LS",runNum);

    h1u0[0]=m_ULhists[0]->ProjectionX("UL0");
    h1l0[0]=m_LShists[0]->ProjectionX("LS0");
    
    TCanvas *c00 = new TCanvas("UL");

    h1u0[0]->SetTitle(ULgraph);
    h1u0[0]->GetXaxis()->SetTitle("Mass GeV");
    h1u0[0]->Draw();
    h1l0[0]->SetLineColor(2);
    h1l0[0]->SetMarkerColor(2);
    h1l0[0]->SetMarkerSize(1.5);
    h1l0[0]->Draw("same");

    TLegend *leg = new TLegend(0.71,0.25,1.0,0.45);
    leg->AddEntry(h1u0[0],"unlike sign");
    leg->AddEntry(h1l0[0],"like sign");
    leg->Draw();

    TString subDir = "ULDistrs/";    

    c00->SaveAs(dir+subDir+ULname);
    c00->Close();

    //TCanvas *c01 = new TCanvas("LS");

    //h1l0[0]->SetTitle(LSgraph);
    //h1l0[0]->GetXaxis()->SetTitle("Mass GeV");
    //h1l0[0]->Draw();

    //c01->SaveAs(dir+LSname);
    //c01->Close();
    
    // h1u0[1]=ULhists[0]->ProjectionY("");
    // h1l0[1]=LShists[0]->ProjectionY("");
    h1u0[0]->Add(h1l0[0],-1);
    
    TCanvas *c0 = new TCanvas("full");
    
    h1u0[0]->SetTitle(fullGraph);
    h1u0[0]->GetYaxis()->SetTitle("Counts");
    h1u0[0]->GetXaxis()->SetTitle("Mass GeV");
    h1u0[0]->Draw();
    
    c0->SaveAs(dir+subDir+fullName);
    c0->Close();

    h1u1[0]=m_ULhists[0]->ProjectionX("UL1");
    h1l1[0]=m_LShists[0]->ProjectionX("LS1");
    // h1u1[1]=ULhists[0]->ProjectionY("");
    // h1l1[1]=LShists[0]->ProjectionY("");
    h1u1[0]->Add(h1l1[0],-1);
    
    TCanvas *c1 = new TCanvas("full");
    
    h1u1[0]->SetTitle(fullGraph);
    h1u1[0]->GetYaxis()->SetTitle("Counts");
    h1u1[0]->GetXaxis()->SetTitle("Mass GeV");
    h1u1[0]->SetAxisRange(0.8,1.3,"X");
    h1u1[0]->Draw();
    
    c1->SaveAs(dir+subDir+f0Name);
    c1->Close();

}

TH1D *primaryFit::funFit(int *bins, int count1, int count2, int count3, int countPass)
{
    // count1 = i
    // count2 = j
    // count3 = k
    // countPass = the total number of times through the loop so far

    // switch these names over to to be dependent on run number and fit type
    TH1D *h1u;
    TH1D *h1l;
    // can make bins a vector so that I don't need to pass it
    h1u=m_ULhists[0]->ProjectionX(m_infoBox.inserts[3*count1+count2][count3][0],bins[0],bins[1],bins[2],bins[3]);
    h1l=m_LShists[0]->ProjectionX(m_infoBox.inserts[3*count1+count2][count3][1],bins[0],bins[1],bins[2],bins[3]);
    h1u->Add(h1l,-1);
    h1u->SetAxisRange(0.85,1.3,"X");
    //h1u->SetAxisRange(axisVal[0],axisVal[1],"Y");
   
    char canvTitle[256];

    TCanvas *c = new TCanvas(m_infoBox.inserts[3*count1+count2][count3][2]);
    
    h1u->SetTitle(m_infoBox.inserts[3*count1+count2][count3][3]);
    h1u->GetXaxis()->SetTitle(m_infoBox.inserts[3*count1+count2][count3][4]);
    // commented out just to look nicer
    //h1u->GetYaxis()->SetTitle(m_infoBox.inserts[3*count1+count2][count3][5]);
    h1u->Draw();
    
    TF1 *fit = new TF1("fit",m_infoBox.functionPieces[count2][0],0,3,m_infoBox.numberOfParams[count2][0]);

    // will try to get rid of the necessity of these with sprintf soon
    if(m_infoBox.numberOfParams[count2][0]==5) fit->SetParNames(m_infoBox.inserts[3*count1+count2][count3][7],
                                        m_infoBox.inserts[3*count1+count2][count3][8],
                                        m_infoBox.inserts[3*count1+count2][count3][9],
                                        m_infoBox.inserts[3*count1+count2][count3][10],
                                        m_infoBox.inserts[3*count1+count2][count3][11]);

    if(m_infoBox.numberOfParams[count2][0]==6) fit->SetParNames(m_infoBox.inserts[3*count1+count2][count3][7],
                                        m_infoBox.inserts[3*count1+count2][count3][8],
                                        m_infoBox.inserts[3*count1+count2][count3][9],
                                        m_infoBox.inserts[3*count1+count2][count3][10],
                                        m_infoBox.inserts[3*count1+count2][count3][11],
                                        m_infoBox.inserts[3*count1+count2][count3][12]);

    if(m_infoBox.numberOfParams[count2][0]==7) fit->SetParNames(m_infoBox.inserts[3*count1+count2][count3][7],
                                        m_infoBox.inserts[3*count1+count2][count3][8],
                                        m_infoBox.inserts[3*count1+count2][count3][9],
                                        m_infoBox.inserts[3*count1+count2][count3][10],
                                        m_infoBox.inserts[3*count1+count2][count3][11],
                                        m_infoBox.inserts[3*count1+count2][count3][12],
                                        m_infoBox.inserts[3*count1+count2][count3][13]);
  
    // control the fit->SetRange() with an if statement

    if(count1==0) fit->SetRange(0.8,1.2); // originally set for (0.9,1.1) 
    if(count1>0) fit->SetRange(0.8,1.2);

    for(int i = 0; i < m_infoBox.numberOfParams[count2][0]; i++)
    {
        fit->SetParLimits(i, m_infoBox.parameterBounds[3*count1+count2][count3][(2*i)],
                             m_infoBox.parameterBounds[3*count1+count2][count3][(2*i)+1]);
    }

    if(m_infoBox.fixMass==1)
    {
        fit->FixParameter(1,f0mass);
    }
    if(m_infoBox.fixWidth==1)
    {
        fit->FixParameter(2,f0width_KK); // will want to make this more variable later
    }

    fit->SetLineColor(46);
    h1u->Fit(fit,"R");
    //fit->Draw("same");
    for(int i=0;i<m_infoBox.numberOfParams[count2][0];i++)
    {
        m_infoBox.ParamVals[3*count1+count2][count3].push_back( fit->GetParameter(i) );
    }
    for(int i=0;i<m_infoBox.numberOfParams[count2][0];i++)
    {

        m_infoBox.ParamVals[3*count1+count2][count3].push_back(fit->GetParError(i));
    }
    TF1 *bkg = new TF1("bkg",m_infoBox.functionPieces[count2][1],0,3,m_infoBox.numberOfParams[count2][1]);
    TF1 *funn= new TF1("funn",m_infoBox.functionPieces[count2][2],0,3,m_infoBox.numberOfParams[count2][2]);
    for(int i = m_infoBox.numberOfParams[count2][2]-1;i<m_infoBox.numberOfParams[count2][0];i++)
    {
        bkg->SetParameter(i,m_infoBox.ParamVals[3*count1+count2][count3][i]);
        m_infoBox.BGrndVals[countPass-1].push_back(m_infoBox.ParamVals[3*count1+count2][count3][i]);
    }
    for(int i = 0;i<m_infoBox.numberOfParams[count2][2]-1;i++)
    {
        funn->SetParameter(i,m_infoBox.ParamVals[3*count1+count2][count3][i]);
        m_infoBox.FunctVals[countPass-1].push_back(m_infoBox.ParamVals[3*count1+count2][count3][i]);
    }
    funn->SetParameter(m_infoBox.numberOfParams[count2][2]-1,h1u->GetBinContent(200)/3);
    m_infoBox.FunctVals[countPass-1].push_back((h1u->GetBinContent(200))/3);
    if(count1==0) bkg->SetRange(0.8,1.2);
    if(count1>0) bkg->SetRange(0.8,1.2);
    if(count1==0) funn->SetRange(0.8,1.2); // originally set for (0.9,1.1)
    if(count1>0) funn->SetRange(0.8,1.2);
    bkg->SetLineColor(9);
    funn->SetLineColor(8);
    
    //bkg->Draw("same");
    //funn->Draw("same");

    //A legend function to be written
    //

    //c->SaveAs(names[6]);
    //c->Close();



    return h1u;
}
void primaryFit::Analysis(struct fitContainer infoBox)
{ 
    // gStyle->SetOptFit(1111);
   
    fillConstr(infoBox); 

    // BEGIN ANALYSIS SECTION //

    int numBackgrounds=3;
    int distrTypes=3;
    int numHistograms=0;

    int thisCount=0;

    // Filling vector and elements for pushing
    vector<TH1D*> vHistPush; 
    vector < vector<TH1D*> > indivHists;

    vector<double> doubPush;
    vector < vector<double> > doubMatr;
    
    // declaring fitting function vectors
    // these all need same structure as indivHists
    vector < vector <TF1*> > tfFitFuncs;
    vector < vector <TF1*> > tfBckGrnds;

    vector <TF1*> tf1Push;

    vector<FitFuncPtr> funcAddVec;
    vector<FitFuncPtr> bckgAddVec;

    // creating 1d histograms for different fit types
    
    // because these for loops already exist, I will just let i j and k be inputs for the fnc, may change later
    for(int i=0;i<distrTypes;i++)
    {
        if(i==0) numHistograms=6; // pt distribution
        if(i>=1) numHistograms=5; // phi distributions

        for(int j=0;j<numBackgrounds;j++)
        {
   
            indivHists.push_back(vHistPush);         
            tfFitFuncs.push_back(tf1Push);
            tfBckGrnds.push_back(tf1Push);

            funcAddVec.push_back(m_infoBox.functionPieces[j][2]);
            bckgAddVec.push_back(m_infoBox.functionPieces[j][1]);

            m_infoBox.ParamVals.push_back(doubMatr); // This isn't working for some reason

            for(int k=0;k<numHistograms;k++)
            {
                //cout << "run number\t" << thisCount << endl;
                thisCount++; 
                //cout << "identify\t" << infoBox.inserts[3*i+j][k][0] << endl;
             
                m_infoBox.ParamVals[3*i+j].push_back(doubPush);
                
                m_infoBox.FunctVals.push_back(doubPush);
                m_infoBox.BGrndVals.push_back(doubPush);

                if(i==0) indivHists[3*i+j].push_back(funFit( ptDisBin[k], i, j, k, thisCount) );
                
                if(i>=1) indivHists[3*i+j].push_back(funFit( phiDisBin[i-1][k], i, j, k, thisCount));

            }

        }
    }

    // filling functions and background
    
    int CountTf1Fill=0;
    
    for(int i=0;i<indivHists.size();i++)
    {
        for(int j=0;j<indivHists[i].size();j++)
        {

            tfFitFuncs[i].push_back(funcFit(m_infoBox.FunctVals[CountTf1Fill],funcAddVec[i],CountTf1Fill));
            tfBckGrnds[i].push_back(bkgFit(m_infoBox.BGrndVals[CountTf1Fill],bckgAddVec[i],CountTf1Fill));
            
            CountTf1Fill++;
        }
    } 


    // making the multiplots
    multHistPlot(indivHists,tfFitFuncs,tfBckGrnds,m_infoBox.canvNames,m_infoBox);

    // yield variables for extracting v2 
    vector < vector<double> > YieldForv2;
    vector < vector<double> > YieldForv2Err;
    vector < vector<double> > xAxisforv2;
    vector < vector<double> > xAxisforv2Err; 

    // extracted v2 and corresponding pT
    vector < vector<double> > v2pTvals;
    vector < vector<double> > v2pTerrs;
    vector < vector<double> > thev2;
    vector < vector<double> > thev2Err; 

    // for gaussian 
    if(m_infoBox.canvNames=="Gauss")
    {

        for(int i=0;i<m_infoBox.ParamVals.size();i++)
        {
            
            YieldForv2.push_back(doubPush);
            YieldForv2Err.push_back(doubPush);
            xAxisforv2.push_back(doubPush);
            xAxisforv2Err.push_back(doubPush);

            //////////// CHANGE ASAP!!! ///////// 
            // this needs to be temporary. I need to calculate the pt and phi 
            // from the histograms
            if(m_infoBox.ParamVals[i].size()==5)
            {
                ////////////////// PLEASE CHANGE THIS TO THE DATA /////////////////////
                xAxisforv2[i].push_back(0.05*PI); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(0.15*PI); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(0.25*PI); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(0.35*PI); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(0.45*PI); xAxisforv2Err[i].push_back(0);
            }
            if(m_infoBox.ParamVals[i].size()==6)
            {
                xAxisforv2[i].push_back(0.25); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(0.75); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(1.25); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(1.75); xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(2.5);  xAxisforv2Err[i].push_back(0);
                xAxisforv2[i].push_back(4);    xAxisforv2Err[i].push_back(0);
            }

            for(int j=0;j<m_infoBox.ParamVals[i].size();j++)
            {
                YieldForv2[i].push_back( m_infoBox.ParamVals[i][j][0]); // the area
                YieldForv2Err[i].push_back( m_infoBox.ParamVals[i][j][ m_infoBox.ParamVals[i][j].size()/2  ] ); // the area error    
            }
        }

    }// end of Gaussian if
    else
    {
        // fill the yields 
        int bgrndCount = 0;
        for(int i = 0; i <indivHists.size() ; i++)
        {
            xAxisforv2.push_back(doubPush);
            xAxisforv2Err.push_back(doubPush);
            YieldForv2.push_back(doubPush);
            YieldForv2Err.push_back(doubPush); 

            // pt yield
            if(indivHists[i].size()==6)
            {
                xAxisforv2[i].push_back(0.25);
                xAxisforv2[i].push_back(0.75);
                xAxisforv2[i].push_back(1.25);
                xAxisforv2[i].push_back(1.75);
                xAxisforv2[i].push_back(2.5);
                xAxisforv2[i].push_back(4);

                // the error is zero
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
            }
            // phi yield
            if(indivHists[i].size()==5)
            {
                xAxisforv2[i].push_back(0.05*PI);
                xAxisforv2[i].push_back(0.15*PI);
                xAxisforv2[i].push_back(0.25*PI);
                xAxisforv2[i].push_back(0.35*PI);
                xAxisforv2[i].push_back(0.45*PI);
                
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
                xAxisforv2Err[i].push_back(0);
            }

            if(indivHists[i].size()<5 || indivHists[i].size()>6)
            {
                printf("number of histogram errors\nnumHists=%lu\n",indivHists[i].size());
            }

            // filling the area
            for(int j = 0; j < indivHists[i].size(); j++)
            {
                YieldForv2[i].push_back( YieldCalc(indivHists[i][j],m_infoBox.functionPieces[i%3][1],m_infoBox.BGrndVals[bgrndCount]) );
                YieldForv2Err[i].push_back( YieldErrorCalc(indivHists[i][j],m_infoBox.functionPieces[i%3][1],m_infoBox.BGrndVals[bgrndCount]) );
                bgrndCount++; 
            }

        } 

    } // end of amplitude fits/else
   
    // plot yields and extract v2
    YieldPlot(YieldForv2,YieldForv2Err,xAxisforv2,xAxisforv2Err,infoBox); 

    // get weighted pT for v2
    weightv2PT();

    // write v2 to file
    char tmpFitName[256];
    for(int i = 0; i < m_infoBox.canvNames.Sizeof();i++)
    {
        tmpFitName[i]=m_infoBox.canvNames[i];
    }

    char fixedMassName[256] = "fixedMass";
    char fixedWidthName[256] = "fixedWidth";

    if(infoBox.fixMass==1)
    {
        strcat(tmpFitName,fixedMassName);    
    }
    if(infoBox.fixWidth==1)
    {
        strcat(tmpFitName,fixedWidthName);
    }

    char tmpRunName[256]; // include this in sprintf name
   
    // should have a different file for each background fit

    char v2FileName[256];
    //TString v2FileName="STAR_results/f0_pt_v2.txt";

    // GET DEVIATIONS TO WRITE!!!!! 
    for(int i = 0; i<m_v2.size(); i++)
    {
        sprintf(v2FileName,"STAR_results/f0_pt_v2_run%i%s",m_infoBox.runNumber,tmpFitName);
        
        if(i==0) strcat(v2FileName,"pol1.txt");
        if(i==1) strcat(v2FileName,"pol2.txt");
        if(i==2) strcat(v2FileName,"exp.txt");
        if(i>2) cout << "m_v2 is too large\tsize="<<m_v2.size() << endl;

        ofstream v2file;
        v2file.open(v2FileName);
 
        v2file << m_weightedv2pT[0]<< " " << m_weightedv2pT[1] << endl; 
        v2file << m_v2[i][0] <<" " << m_v2[i][1] << endl;
        v2file << m_v2Err[i][0] <<" " << m_v2Err[i][1] << endl;

        v2file.close();
    } 
    
    
    // deleting pointers
    for(int i = 0;i<vHistPush.size();i++)
    {
        delete vHistPush[i];
    }
    for(int i = 0;i<indivHists.size();i++)
    {

        for(int j = 0;j<indivHists[i].size();j++)
        {
            delete indivHists[i][j];   
        }
    }
    

    return;
}
// various plotting functions
void primaryFit::multHistPlot(vector < vector<TH1D*> > &hists,vector < vector<TF1*> > &TheFit,vector < vector<TF1*> > &TheBck,TString type,fitContainer infoBox)
{
    char canv[256];
    char save[256]; 

    char TMPnames[256];

    for(int i=0;i<type.Sizeof();i++)
    {
        TMPnames[i]=type[i];
    }

    char fixedMassName[256] = "fixedMass";
    char fixedWidthName[256] = "fixedWidth";
    
    if(infoBox.fixMass==1)
    {
        strcat(TMPnames,fixedMassName);
    }
    if(infoBox.fixWidth==1)
    {
        strcat(TMPnames,fixedWidthName);
    }

    for(int i=0;i<hists.size();i++)
    {
    
        // Will NEED to add more options for the names later
        // maybe just a string argument at the top
        sprintf(canv,"%scm%i",TMPnames,i+1);  
        sprintf(save,"mainPlots/massDistrs/run%i%smulPlot%i.pdf",infoBox.runNumber,TMPnames,i+1); // include run later

        TCanvas *cm = new TCanvas(canv);
        
        if(hists[i].size()==4 || hists[i].size()==3) cm->Divide(2,2); 
        if(hists[i].size()==5 || hists[i].size()==6) cm->Divide(3,2);
        
        for(int j=0;j<hists[i].size();j++)
        {
            cm->cd(j+1);
            hists[i][j]->Draw();

            TheFit[i][j]->Draw("same");
            TheBck[i][j]->Draw("same");

            TheFit[i][j]->SetLineColor(8);
            TheBck[i][j]->SetLineColor(9);
        }

        cm->SaveAs(save);
        cm->Close();
    }

    return;
}
void primaryFit::YieldPlot(vector < vector<double> > &yields, vector < vector<double> > &yieldErr, vector < vector <double> > &xAxis, 
        vector < vector<double> > &xAxisErr, fitContainer infoBox)
{

    // this will create 3 multigraphs per call
    int numYieldPlots=5;
    char IndivCanvNames[numYieldPlots][256];
    char multiCanvNames[256];

    char IndivCanvSaves[256];// probably won't use this one
    char multiCanvSaves[256];

    char funcName[256];

    for(int i = 0; i < infoBox.canvNames.Sizeof() ; i++)
    {
        funcName[i]=infoBox.canvNames[i];
    }

    char fixedMassName[256] = "fixedMass";
    char fixedWidthName[256] = "fixedWidth";
    
    if(infoBox.fixMass==1)
    {
        strcat(funcName,fixedMassName);
    }
    if(infoBox.fixWidth==1)
    {
        strcat(funcName,fixedWidthName);
    }

    if( yields.size()/3!=3 ) cout << "MORE YIELDS THAN EXPECTED" << endl;
    
    string bkg1 = "pol1";
    string bkg2 = "pol2";
    string bkg3 = "exp";

    char loader[5];

    vector<double> v2push;

    // the order for the yields is based on indivHists order
    for(int i = 0; i < (yields.size()/3) ; i++)
    {
        if(i==0) for(int j = 0; j<4; j++) loader[j]=bkg1[j];
        if(i==1) for(int j = 0; j<4; j++) loader[j]=bkg2[j];
        if(i==2) for(int j = 0; j<3; j++) loader[j]=bkg3[j];
   
        double tmpYieldDec[6];
        double tmpYieldDecErr[6];
        double tmpxAxisDec[6];
        double tmpxAxisDecErr[6];

        double tmpYieldlpt[5];
        double tmpYieldlptErr[5];
        double tmpYieldhpt[5];
        double tmpYieldhptErr[5];
        double tmpPhixAxis[5];
        double tmpPhixAxisErr[5];

        for(int j = 0;j<yields[i].size();j++)
        {
            tmpYieldDec[j]=yields[i][j];
            tmpYieldDecErr[j]=yieldErr[i][j];
            tmpxAxisDec[j]=xAxis[i][j];
            tmpxAxisDecErr[j]=xAxisErr[i][j]; 
        }
        for(int j = 0;j<yields[i+3].size();j++)
        {
            tmpYieldlpt[j]=yields[i+3][j];
            tmpYieldlptErr[j]=yieldErr[i+3][j];
            tmpYieldhpt[j]=yields[i+6][j];
            tmpYieldhptErr[j]=yieldErr[i+6][j];
            tmpPhixAxis[j]=xAxis[i+3][j];
            tmpPhixAxisErr[j]=xAxisErr[i+3][j];
        }

        sprintf(multiCanvNames,"run%i%syields%i",infoBox.runNumber,funcName,i);
        sprintf(multiCanvSaves,"mainPlots/yields/run%i%syields%s.pdf",infoBox.runNumber,funcName,loader);

        sprintf(IndivCanvNames[0],"dN/dpT vs pT run%i %s fit %s bkg",infoBox.runNumber,funcName,loader);       
        sprintf(IndivCanvNames[1],"dN/dpT vs pT (log) run%i %s fit %s bkg",infoBox.runNumber,funcName,loader);       
        sprintf(IndivCanvNames[2],"dN/pT dpT vs pT run%i %s fit %s bkg",infoBox.runNumber,funcName,loader);       
        sprintf(IndivCanvNames[3],"dN/d#phi 0<pT<2GeV run%i %s fit %s bkg",infoBox.runNumber,funcName,loader);       
        sprintf(IndivCanvNames[4],"dN/d#phi 2<pT<5GeV run%i %s fit %s bkg",infoBox.runNumber,funcName,loader);       
        
        TGraphErrors *yieldGraphs[5];
       
        TCanvas *cmg = new TCanvas(multiCanvNames);
        //cmg->Divide(3,2);
        //
        // TEMPORARY EDITS
        //
        // cmg->Divide(2,1);

        // variable for (dN/pT)/dpT
        double yieldOverPt[6];
        for(int smallCounter = 0; smallCounter < yields[i].size(); smallCounter++)
        {
            yieldOverPt[smallCounter]=(yields[i][smallCounter]/xAxis[i][smallCounter]);   
        }
       
        cmg->cd(1);
        yieldGraphs[3] = new TGraphErrors(yields[i+3].size(),tmpPhixAxis,tmpYieldlpt,tmpPhixAxisErr,tmpYieldlptErr);// v2 fit
        yieldGraphs[3]->SetMarkerStyle(21);
        yieldGraphs[3]->SetTitle(IndivCanvNames[3]);
        yieldGraphs[3]->GetXaxis()->SetTitle("#phi");
        yieldGraphs[3]->GetYaxis()->SetTitle("dN/d#phi");
        yieldGraphs[3]->Draw();
        TF1 *fv2lpt = new TF1("v2lptfit",v2cosfit,0,3,2);
        fv2lpt->SetParNames("Amp","v2");
        //need to get the v2 parameters from this function still
        yieldGraphs[3]->Fit(fv2lpt,"R");
        //fv2lpt->Draw("same");       
        
        m_v2.push_back(v2push);
        m_v2Err.push_back(v2push);
        m_v2[i].push_back( fv2lpt->GetParameter(1) );
        m_v2Err[i].push_back( fv2lpt->GetParError(1) );
        
        //cmg->cd(2);
        yieldGraphs[4] = new TGraphErrors(yields[i+6].size(),tmpPhixAxis,tmpYieldhpt,tmpPhixAxisErr,tmpYieldhptErr);// v2 fit
        yieldGraphs[4]->SetMarkerStyle(21);
        yieldGraphs[4]->SetTitle(IndivCanvNames[4]);
        yieldGraphs[4]->GetXaxis()->SetTitle("#phi");
        yieldGraphs[4]->GetYaxis()->SetTitle("dN/d#phi");
        yieldGraphs[4]->Draw();
        TF1 *fv2hpt = new TF1("v2hptfit",v2cosfit,0,3,2);
        fv2hpt->SetParNames("Amp","v2");
        yieldGraphs[4]->Fit(fv2hpt,"R");
        fv2hpt->Draw("same");
        
        m_v2[i].push_back( fv2hpt->GetParameter(1) );
        m_v2Err[i].push_back( fv2hpt->GetParError(1) );
 
        cmg->SaveAs(multiCanvSaves);
        cmg->Close();

    }

    return;
}
// This function creates the Pt values for the v2 vs pT  
void primaryFit::weightv2PT() 
{
    // this needs to be redone. Do the weighting by calling each bin and adding entries in the same region
    // ->GetBinCenter(),, ->GetBinContent
    // Y projection has pt

    TH1D *ulHigh;
    TH1D *ulLow;
    TH1D *lsHigh;
    TH1D *lsLow;

    double valueTopHigh=0;
    double valueTopLow=0;
    double valueBotHigh=0;
    double valueBotLow=0;
      
    ulLow = m_ULhists[0]->ProjectionY(); 
    ulHigh = m_ULhists[0]->ProjectionY();
    lsLow = m_LShists[0]->ProjectionY();
    lsHigh = m_LShists[0]->ProjectionY();

    ulLow->Add(lsLow,-1);
    ulHigh->Add(lsHigh,-1);

    // the for loop numbers are predetermined by binning and ptCut choices
    for(int i=1; i<9;i++)
    {
        if(ulLow->GetBinCenter(i)>=2 ) cout << "v2 pT calculation error low" << endl;
        valueTopLow = valueTopLow + ( ulLow->GetBinContent(i) )*( ulLow->GetBinCenter(i) );
        valueBotLow = valueBotLow + ( ulLow->GetBinContent(i) );
    }

    for(int i=9; i<21;i++)
    {
        if(ulHigh->GetBinCenter(i)<2 ) cout << "v2 pT calc error high" << endl; 
        valueTopHigh = valueTopHigh + ( ulHigh->GetBinContent(i) )*( ulHigh->GetBinCenter(i) );
        valueBotHigh = valueBotHigh + ( ulHigh->GetBinContent(i) );
    }

    double ptLow = valueTopLow/valueBotLow;
    double ptHigh = valueTopHigh/valueBotHigh;

    m_weightedv2pT.push_back(ptLow);
    m_weightedv2pT.push_back(ptHigh);

    // CALCULATE THE ERROR 

    return;
}

double primaryFit::YieldCalc(TH1D *histo,FitFuncPtr bkgFunc,vector<double> &bkgVals)
{
    // bin 160 : mass = 0.7975
    // bin 260 : mass = 1.2975
    // bin 198 : mass = 0.9875

    // bin 181 : mass = 0.90250
    // bin 220 : mass = 1.09750

    int binStart = 181;
    int binEnd = 220;

    // the bkg is only defined from 0.9 to 1.1

    TF1 *bkg = new TF1("bkg",bkgFunc,0,3,bkgVals.size());
    
    for(int i = 0; i<bkgVals.size();i++)
    {
        bkg->SetParameter(i,bkgVals[i]);
    }

    // width and peak aren't important based on how the functions are designed so far
    double width = ( histo->GetBinCenter(161) - histo->GetBinCenter(160) );
    double peak = histo->GetBinCenter(198);

    double xPoint = 0.9;

    double yield = 0;

    double yieldErr = 0; // worth looking at this again

    for(int i = binStart; i <  binEnd+1 ;i++)
    {
        xPoint = histo->GetBinCenter(i);
        yield = yield + ( histo->GetBinContent(i) ) - ( bkg->Eval(xPoint) );
        yieldErr = yieldErr + ( histo->GetBinError(i) )*( histo->GetBinError(i) )*width*width;
    }
 
    yield = yield*width;

    yieldErr = TMath::Sqrt(yieldErr);

    return yield;
}
// there might be a more elegant way, but this works for now
double primaryFit::YieldErrorCalc(TH1D *histo,FitFuncPtr bkgFunc,vector<double> &bkgVals)
{
    // bin 160 : mass = 0.7975
    // bin 260 : mass = 1.2975
    // bin 198 : mass = 0.9875

    // bin 181 : mass = 0.90250
    // bin 220 : mass = 1.09750

    int binStart = 181;
    int binEnd = 220;

    // the bkg is only defined from 0.9 to 1.1

    TF1 *bkg = new TF1("bkg",bkgFunc,0,3,bkgVals.size());
    
    for(int i = 0; i<bkgVals.size();i++)
    {
        bkg->SetParameter(i,bkgVals[i]);
    }

    // width and peak aren't important based on how the functions are designed so far
    double width = ( histo->GetBinCenter(161) - histo->GetBinCenter(160) );
    double peak = histo->GetBinCenter(198);

    double xPoint = 0.9;

    double yield = 0;

    double yieldErr = 0; // worth looking at this again

    for(int i = binStart; i <  binEnd+1 ;i++)
    {
        xPoint = histo->GetBinCenter(i);
        yield = yield + ( histo->GetBinContent(i) ) - ( bkg->Eval(xPoint) );
        yieldErr = yieldErr + ( histo->GetBinError(i) )*( histo->GetBinError(i) )*width*width;
    }
 
    yield = yield*width;

    yieldErr = TMath::Sqrt(yieldErr);

    return yieldErr;
}

void primaryFit::f0splits()
{
    // creates and saves split histograms like funFit, but without the fit 

    return;
}


void primaryFit::singleParticlepT()
{
    // creates the single particle pT distribution    

    return;
}

// Below are HARD CODED functions that won't change by variable choice
// choosing to make them part of the class, may change that later
void primaryFit::compareRuns_01()
{
    // compares UL-LS Mass distributions for run11 and run16

    TFile *f01 = new TFile("DATA/OLDDATA/run11Data.root");
    //TFile *f02 = new TFile("DATA/OLDDATA/run16Data.root");
    TFile *f02 = new TFile("DATA/OLDDATA/run16test.root");

    int centBins = 9;
    int cent = 6; // choosing my centralities to add

    char ULvar01[50];
    char LSvar01[50];

    TH3D *ULhists01[9];
    TH3D *LShists01[9];
    TH3D *ULhists02[9];
    TH3D *LShists02[9];
    TH1F *vert01;
    TH1F *vert02;

    for(int i=0;i<centBins;i++){
        sprintf(ULvar01,"TMdphiUL0c%i",i+1);
        sprintf(LSvar01,"TMdphiLS0c%i",i+1);
        ULhists01[i]=(TH3D*)f01->Get(ULvar01);
        LShists01[i]=(TH3D*)f01->Get(LSvar01); 
    }
    for(int i=0;i<cent;i++){
        ULhists01[0]->Add(ULhists01[i+1],1);
        LShists01[0]->Add(LShists01[i+1],1);
    }

    char ULvar02[50];
    char LSvar02[50]; 
    for(int i=0;i<centBins;i++){
        sprintf(ULvar02,"TMdphiUL0c%i",i+1);
        sprintf(LSvar02,"TMdphiLS0c%i",i+1);
        ULhists02[i]=(TH3D*)f02->Get(ULvar02);
        LShists02[i]=(TH3D*)f02->Get(LSvar02); 
    }
    for(int i=0;i<cent;i++){
        ULhists02[0]->Add(ULhists02[i+1],1);
        LShists02[0]->Add(LShists02[i+1],1);
    }

    // getting the number of events
    vert01 = (TH1F*)f01->Get("centralityw");
    vert02 = (TH1F*)f02->Get("centralityw");

    // may get type conversion error
    double numEventsRun11 = vert01->Integral(2,7);
    double numEventsRun16 = vert02->Integral(2,7);

    //cout << numEventsRun11 << endl;
    //cout << numEventsRun16 << endl;

    //return;

    TH1D *h1u01;
    TH1D *h1l01;
    TH1D *h1u02;
    TH1D *h1l02;

    //h1u01 = ULhists01[0]->ProjectionX("UL01");
    //h1l01 = LShists01[0]->ProjectionX("LS01");
    
    //h1u02 = ULhists02[0]->ProjectionX("UL02");
    //h1l02 = LShists02[0]->ProjectionX("LS02");

    h1u01 = ULhists01[0]->ProjectionX("UL01",1,8,1,20);
    h1l01 = LShists01[0]->ProjectionX("LS01",1,8,1,20);
    
    h1u02 = ULhists02[0]->ProjectionX("UL02",1,8,1,20);
    h1l02 = LShists02[0]->ProjectionX("LS02",1,8,1,20);


    // need the phi and pt bins

    double oldEntry = 0;
    double updateEntry = 0;

    double sumCheck = 0;

    h1u01->Scale(1/numEventsRun11);
    h1l01->Scale(1/numEventsRun11);
    h1u02->Scale(1/numEventsRun16);
    h1l02->Scale(1/numEventsRun16);

    //for(int i=0;i<h1u01->GetNbinsX();i++)
    //{
    //    oldEntry = h1u01->GetBinContent(i);
    //    sumCheck = sumCheck + oldEntry*h1u01->GetBinWidth(i);
    //    updateEntry = oldEntry/numEventsRun11;
    //    h1u01->SetBinContent(i,updateEntry);
    //}

    //for(int i=0;i<h1l01->GetNbinsX();i++)
    //{
    //    oldEntry = h1l01->GetBinContent(i);
    //    sumCheck = sumCheck + oldEntry*h1l01->GetBinWidth(i);
    //    updateEntry = oldEntry/numEventsRun11;
    //    h1l01->SetBinContent(i,updateEntry);
    //}

    //for(int i=0;i<h1u02->GetNbinsX();i++)
    //{
    //    oldEntry = h1u02->GetBinContent(i);
    //    sumCheck = sumCheck + oldEntry*h1u02->GetBinWidth(i);
    //    updateEntry = oldEntry/numEventsRun16;
    //    h1u02->SetBinContent(i,updateEntry);
    //}

    //for(int i=0;i<h1l02->GetNbinsX();i++)
    //{
    //    oldEntry = h1l02->GetBinContent(i);
    //    //sumCheck = sumCheck + oldEntry*h1l02->GetBinWidth(i);
    //    updateEntry = oldEntry/numEventsRun16;
    //    h1l02->SetBinContent(i,updateEntry);
    //}

    TCanvas *c00 = new TCanvas("UL+LSdistr");

    h1u02->SetAxisRange(-10,20,"Y");

    h1u02->Draw("hist");
    h1l01->SetMarkerColor(2);
    h1l01->SetLineColor(2);
    h1u02->SetMarkerColor(3);
    h1u02->SetLineColor(3);
    h1l02->SetMarkerColor(4);
    h1l02->SetLineColor(4);
    h1l01->Draw("hist same");
    h1u01->Draw("hist same");
    h1l02->Draw("hist same");

    c00->SaveAs(dir+"Run11andRun16ULandLS.pdf");
    c00->Close();

    h1u01->Add(h1l01,-1);
    h1u02->Add(h1l02,-1);

    // normalize run11 entries by events
    // ALREADY NORMALIZED THE RUNS

    sumCheck = 0;
    
    oldEntry = 0;
    updateEntry = 0;

    TCanvas *cLarge = new TCanvas("01");

    h1u01->SetAxisRange(-0.05,0.1,"Y");
    
    h1u01->Draw("HIST");
    h1u02->Draw("HIST same");
    h1u02->SetLineColor(8);  

    cLarge->SaveAs("mainPlots/Run11Run16CompareWhole.pdf");
    cLarge->Close(); 

    // for an arbitrary scale up
    oldEntry=0;
    updateEntry=0;
    for(int i = 0;i<h1u02->GetNbinsX();i++)
    {
        oldEntry = h1u02->GetBinContent(i);
        updateEntry=oldEntry*2;
        h1u02->SetBinContent(i,updateEntry);
    }

    TCanvas *cSmall = new TCanvas("02");

    h1u01->SetAxisRange(-1e-3,.05,"Y");

    h1u01->Draw("HIST");
    h1u02->Draw("HIST same");
    h1u01->SetAxisRange(0.8,1.3,"X");
    h1u02->SetLineColor(8);

    cSmall->SaveAs("mainPlots/Run11Run16Comparef0.pdf");
    cSmall->Close();

    if(f01) delete f01;
    if(f02) delete f02;

    return;
}

// creates plots examining UL, LS, and UL-LS for different centralities for a given root file
void primaryFit::compareRuns_02()
{
    TFile *f0 = new TFile("DATA/OLDDATA/run16Data.root");

    TFile *f1 = new TFile("DATA/OLDDATA/run11Data.root");

    int centBins = 9;

    char ULvar[256];
    char LSvar[256];

    TH3D *ULhists[centBins];
    TH3D *LShists[centBins];

    TH1D *ul1d[centBins];
    TH1D *ls1d[centBins];

    TH3D *run11ULhists[centBins];
    TH3D *run11LShists[centBins];

    TH1D *run11ul1d[centBins];
    TH1D *run11ls1d[centBins];

    // need to do integral for centralityw from the correct centralities I'm examining
    TH1F *multNum;
    multNum = (TH1F*)f0->Get("RefMult");
    double numEvents = multNum->GetEntries();

    for(int i=0;i<centBins;i++){
        sprintf(ULvar,"TMdphiUL0c%i",i+1);
        sprintf(LSvar,"TMdphiLS0c%i",i+1);

        ULhists[i]=(TH3D*)f0->Get(ULvar);
        LShists[i]=(TH3D*)f0->Get(LSvar); 
    
        ul1d[i] = ULhists[i]->ProjectionX();
        ls1d[i] = LShists[i]->ProjectionX();

        run11ULhists[i]=(TH3D*)f1->Get(ULvar);
        run11LShists[i]=(TH3D*)f1->Get(LSvar); 
    
        run11ul1d[i] = ULhists[i]->ProjectionX();
        run11ls1d[i] = LShists[i]->ProjectionX();
    }

    // currently the code does not do anything with the run 11 distributions.
    // Add new files and canvases to store them in

    char savePlaceULandLS[256];
    char savePlaceULminusLS[256];
    char canvasName[256];
    char ULminLScanvname[256];

    char run11SavePlaceULLS[256];
    char run11SavePlaceULminLS[256];
    char run11CanvName[256];
    char run11ULminLScanv[256];

    // creating the plots
    for(int i = 0;i<centBins;i++)
    {
        sprintf(canvasName,"c%i",i);
        sprintf(savePlaceULandLS,"mainPlots/testPlots/ULLScent%i.pdf",i);
        sprintf(savePlaceULminusLS,"mainPlots/testPlots/UL-LScent%i.pdf",i);
        sprintf(ULminLScanvname,"c_%i",i);

        sprintf(run11SavePlaceULLS,"mainPlots/testPlots/run11ULLScent%i",i);
        //sprintf();

        TCanvas *c0 = new TCanvas(canvasName);

        double oldEntry = 0;
        double updEntry = 0;

        double upperBound = 0.1;
        double lowerBound = -0.01;

        // for some reason, the bins are filled up to N+1

        // normalizing the mass distributions
        for(int j = 0;j<ul1d[i]->GetNbinsX()+1;j++)
        {
            oldEntry = ul1d[i]->GetBinContent(j);
            updEntry = oldEntry/numEvents;
            ul1d[i]->SetBinContent(j,updEntry);
        }

        oldEntry = 0;
        updEntry = 0;

        for(int j = 0;j<ls1d[i]->GetNbinsX()+1;j++)
        {
            oldEntry = ls1d[i]->GetBinContent(j);
            updEntry = oldEntry/numEvents;
            ls1d[i]->SetBinContent(j,updEntry);    
        }

        upperBound=(ul1d[i]->GetBinContent(ul1d[i]->GetMaximumBin() ) );
        lowerBound=(ul1d[i]->GetBinContent(ul1d[i]->GetMinimumBin() ) );

        ul1d[i]->SetAxisRange(lowerBound,upperBound,"Y");
        ul1d[i]->Draw("hist");
        ls1d[i]->Draw("hist same");
        ls1d[i]->SetMarkerColor(2);
        ls1d[i]->SetLineColor(2);
        
        c0->SaveAs(savePlaceULandLS);
        c0->Close();

        upperBound = 0.1;
        lowerBound = -0.1;

        TCanvas *c1 = new TCanvas(ULminLScanvname);

        ul1d[i]->Add(ls1d[i],-1);
       
        upperBound=(ul1d[i]->GetBinContent(ul1d[i]->GetMaximumBin() ) );
        lowerBound=(ul1d[i]->GetBinContent(ul1d[i]->GetMinimumBin() ) );

        ul1d[i]->SetAxisRange(lowerBound,upperBound,"Y");
        ul1d[i]->Draw("hist");
        
        c1->SaveAs(savePlaceULminusLS);
        c1->Close();

        delete c0; 
        delete c1;   
    }

    /*
    delete f0;
    delete ULhists[centBins];
    delete LShists[centBins];
    delete ul1d[centBins];
    delete ls1d[centBins];
    */

    return;
} 
// made some single charged particle distributions, this function will plot them
void primaryFit::SinglePartDistrs()
{
    // may want to turn things into functions if this has to get any longer
    // I may make it so that this program performs the functions of the previous two compare runs

    // files
    TFile *f11 = new TFile("DATA/OLDDATA/run11test.root");
    //TFile *f16 = new TFile("DATA/OLDDATA/run16eventtest_03.root");
    //TFile *f11 = new TFile("DATA/OLDDATA/run11pidTest.root");
    //TFile *f11 = new TFile("DATA/Rroot");
    //TFile *f16 = new TFile("DATA/OLDDATA/run16pidTest.root");
    TFile *f16 = new TFile("DATA/Run16_full_correct/run16data200.root");
    //TFile *f16 = new TFile("DATA/OLDDATA/run16noTOFglob99job.root");
    //TFile *f11 = new TFile("DATA/OLDDATA/run11Data.root");
    // TFile *f16 = new TFile("DATA/OLDDATA/run16Data.root");
    //TFile *fprim16 = new TFile("DATA/OLDDATA/run16primary.root"); // primary tracks for run16

    TFile *f14 = new TFile("DATA/Run14_test02/run14data200.root");
    //TFile *f14 = new TFile("DATA/Run14_noHFT_03/run14data200.root");

    int centBins = 9;

    double run11TotEvents;
    double run16TotEvents;
    double run14TotEvents;


    double run11centEvents[centBins];
    double run16centEvents[centBins];
    double run14centEvents[centBins];

    TH3D *run11ULhists3d[centBins];
    TH3D *run11LShists3d[centBins];
        
    TH3D *run16ULhists3d[centBins];
    TH3D *run16LShists3d[centBins];
    
    TH3D *run14ULhists3d[centBins];
    TH3D *run14LShists3d[centBins];

    // getting the histograms
    
    char ULName[256];
    char LSName[256];
    
    for(int i = 0; i <centBins ; i++)
    {
        sprintf(ULName,"TMdphiUL0c%i",i+1);
        sprintf(LSName,"TMdphiLS0c%i",i+1);

        run11ULhists3d[i] = (TH3D*)f11->Get(ULName);
        run11LShists3d[i] = (TH3D*)f11->Get(LSName);
        run16ULhists3d[i] = (TH3D*)f16->Get(ULName);
        run16LShists3d[i] = (TH3D*)f16->Get(LSName);
        run14ULhists3d[i] = (TH3D*)f14->Get(ULName);
        run14LShists3d[i] = (TH3D*)f14->Get(LSName);
    }

    // adding the centralities from 20 to 80%
    for(int i = 1; i<6;i++)
    {
        run11ULhists3d[0]->Add(run11ULhists3d[i],1);
        run11LShists3d[0]->Add(run11LShists3d[i],1);
        run16ULhists3d[0]->Add(run16ULhists3d[i],1);
        run16LShists3d[0]->Add(run16LShists3d[i],1);
        run14ULhists3d[0]->Add(run14ULhists3d[i],1);
        run14LShists3d[0]->Add(run14LShists3d[i],1);
    }

    // bin=2 is the first centrality region (the most peripheral)
    TH1F *run11Centralities = (TH1F*)f11->Get("centralityw");
    TH1F *run16Centralities = (TH1F*)f16->Get("centralityw");
    TH1F *run14Centralities = (TH1F*)f14->Get("centralityw");

    run11TotEvents = run11Centralities->GetBinContent(1);
    run16TotEvents = run16Centralities->GetBinContent(1);
    run14TotEvents = run14Centralities->GetBinContent(1);


    for(int i = 2;i<10;i++) // numbers because of structure of centrality events
    {
        run11centEvents[i-2] = run11Centralities->GetBinContent(i);
        run16centEvents[i-2] = run16Centralities->GetBinContent(i);
        run14centEvents[i-2] = run14Centralities->GetBinContent(i);
    }
   
    double eventNormrun11 = run11Centralities->Integral(2,7);
    double eventNormrun16 = run16Centralities->Integral(2,7);
    double eventNormrun14 = run14Centralities->Integral(2,7);

    // other plots

    // mass distribution ratios
    
    run11ULhists3d[0]->Add(run11LShists3d[0],-1);
    run16ULhists3d[0]->Add(run16LShists3d[0],-1);
    run14ULhists3d[0]->Add(run14LShists3d[0],-1);

    TH1D *run11ullsphiDistr[2][5];
    TH1D *run16ullsphiDistr[2][5];
    TH1D *run14ullsphiDistr[2][5];

    char run11projNames[256];
    char run16projNames[256];
    char run14projNames[256];
    char canvProjNames[256];
    char canvProjSaves[256];

    char canvDivNames[256];
    char canvDivSaves[256];

    for(int i = 0; i < 2; i++)
    {

            
        sprintf(canvProjNames,"canvProj%i",i);
        TCanvas *cp = new TCanvas(canvProjNames);
        cp->Divide(3,2);

        // comparing f0 directly (no ratio)
        for(int j = 0;j<5;j++)
        {
            sprintf(run11projNames,"run11proj%iname%i",i,j);
            sprintf(run16projNames,"run16proj%iname%i",i,j);
            sprintf(run14projNames,"run14proj%iname%i",i,j);
    
            run11ullsphiDistr[i][j] = run11ULhists3d[0]->ProjectionX(run11projNames,phiDisBin[i][j][0],phiDisBin[i][j][1],phiDisBin[i][j][2],phiDisBin[i][j][3]);
            run16ullsphiDistr[i][j] = run16ULhists3d[0]->ProjectionX(run16projNames,phiDisBin[i][j][0],phiDisBin[i][j][1],phiDisBin[i][j][2],phiDisBin[i][j][3]);
            run14ullsphiDistr[i][j] = run14ULhists3d[0]->ProjectionX(run14projNames,phiDisBin[i][j][0],phiDisBin[i][j][1],phiDisBin[i][j][2],phiDisBin[i][j][3]);

            // set entries for normalization
            double entryBE=0;
            double entryAF=0;
            
            for(int k = 0; k < ( run11ullsphiDistr[i][j]->GetNbinsX() +1);k++)
            {
                entryBE = run11ullsphiDistr[i][j]->GetBinContent(k);
                entryAF = entryBE/eventNormrun11;
                run11ullsphiDistr[i][j]->SetBinContent(k,entryAF);

                entryBE=0;
                entryAF=0;

                entryBE = run16ullsphiDistr[i][j]->GetBinContent(k);
                entryAF = entryBE/eventNormrun16;
                run16ullsphiDistr[i][j]->SetBinContent(k,entryAF);

                entryBE=0;
                entryAF=0;

                entryBE = run14ullsphiDistr[i][j]->GetBinContent(k);
                entryAF = entryBE/eventNormrun14;
                run14ullsphiDistr[i][j]->SetBinContent(k,entryAF);

                entryBE=0;
                entryAF=0;
            }            
            
            // draw
            cp->cd(j+1); // Add A Legend!!!
            
            run11ullsphiDistr[i][j]->Draw("hist");
            run11ullsphiDistr[i][j]->SetAxisRange(0.8,1.1,"X");

            run16ullsphiDistr[i][j]->SetLineColor(2);
            run16ullsphiDistr[i][j]->SetMarkerColor(2);
            run16ullsphiDistr[i][j]->Draw("hist same");

            run14ullsphiDistr[i][j]->SetLineColor(3);
            run14ullsphiDistr[i][j]->SetMarkerColor(3);
            run14ullsphiDistr[i][j]->Draw("hist same");
        }
     
        sprintf(canvProjSaves,"mainPlots/massDistrs/f0comps%i.pdf",i);
        cp->SaveAs(canvProjSaves);
        cp->Close();
        
         
        sprintf(canvDivNames,"cRatio%i",i);
        sprintf(canvDivSaves,"mainPlots/massDistrs/f0ratio%i.pdf",i);

        TCanvas *cr = new TCanvas(canvDivNames);
        cr->Divide(3,2);

        // This is a ratio
        for(int j = 0; j < 5;j++ )
        {
            double entryBE11 = 0;
            double entryAF11 = 0;
            double entryBE16 = 0;
            double entryAF16 = 0;
            double entryBE14 = 0;
            double entryAF14 = 0;
            double entryFin = 0;

            for(int k = 0; k < (run11ULhists3d[i]->GetNbinsX() + 1 );k++)
            {
                entryBE11 = run11ullsphiDistr[i][j]->GetBinContent(k);
                entryAF11 = entryBE11/eventNormrun11;

                entryBE16 = run16ullsphiDistr[i][j]->GetBinContent(k);
                entryAF16 = entryBE16/eventNormrun16;

                entryBE14 = run14ullsphiDistr[i][j]->GetBinContent(k);
                entryAF14 = entryBE14/eventNormrun14;
 
                //entryFin = entryAF11/entryAF16;
                entryFin = entryAF11/entryAF14;

                if(entryAF16==0) entryFin=-999;

                run11ullsphiDistr[i][j]->SetBinContent(k,entryFin);
            } 

            cr->cd(j+1);
            run11ullsphiDistr[i][j]->SetAxisRange(0.8,1.3,"X");
            run11ullsphiDistr[i][j]->Draw("hist");
        }

        cr->SaveAs(canvDivSaves);
        cr->Close();

    }

    // getting pt distribution of f0
        // UL and LS have already been subtracted in 3d plots
    TH1D *run11pTf0;
    TH1D *run16pTf0;

    //TFile *f16TOFTPC = new TFile("DATA/OLDDATA/run16eventtest_03.root");
    TFile *f16TOFTPC = new TFile("DATA/Run16_set01/run16Data_set01.root");
    TH1D *run16pTf0TOFTPC;
    TH3D *run16ULTOFTPC[centBins];
    TH3D *run16LSTOFTPC[centBins];

    for(int i = 0;i<centBins;i++)
    {
        sprintf(ULName,"TMdphiUL0c%i",i+1);
        sprintf(LSName,"TMdphiLS0c%i",i+1);

        run16ULTOFTPC[i] = (TH3D*)f16TOFTPC->Get(ULName);
        run16LSTOFTPC[i] = (TH3D*)f16TOFTPC->Get(LSName);
    }

    for(int i = 1;i<6;i++)
    {
        run16ULTOFTPC[0]->Add(run16ULTOFTPC[i],1);
        run16LSTOFTPC[0]->Add(run16LSTOFTPC[i],1);
    }

    run16ULTOFTPC[0]->Add(run16LSTOFTPC[0],-1);
    run16pTf0TOFTPC = run16ULTOFTPC[0]->ProjectionY("run16pTf0TofTpc",180,220,1,20);

    TH1F *run16hTofTpcNumEvents;
    run16hTofTpcNumEvents = (TH1F*)f16TOFTPC->Get("centralityw");
    double run16TOFTPCnumEvents=run16hTofTpcNumEvents->Integral(2,7);

    run11pTf0 = run11ULhists3d[0]->ProjectionY("run11pTf0",180,220,1,20);
    run16pTf0 = run16ULhists3d[0]->ProjectionY("run16pTf0",180,220,1,20);

    // scaling entries
  
    double entryBE11=0;
    double entryAF11=0;
    double entryBE16=0;
    double entryAF16=0;

    //printf("numEvents11=%f\tnumEvents16=%f\n",eventNormrun11,eventNormrun16);

    run11pTf0->Scale(1/eventNormrun11);
    run16pTf0->Scale(1/eventNormrun16);
    run16pTf0TOFTPC->Scale(1/run16TOFTPCnumEvents);

    // the run11pTf0->GetNbins isn't working....
    //for(int i = 1; i< 21; i++)
    //{
    //    entryBE11 = run11pTf0->GetBinContent(i);   
    //    entryAF11 = entryBE11/eventNormrun11;
    //    run11pTf0->SetBinContent(i,entryAF11);

    //    entryBE16 = run16pTf0->GetBinContent(i);
    //    entryAF16 = entryBE16/eventNormrun16;
    //    run16pTf0->SetBinContent(i,entryAF16);

    //    entryBE16 = run16pTf0TOFTPC->GetBinContent(i);
    //    entryAF16 = entryBE16/run16TOFTPCnumEvents;
    //    run16pTf0TOFTPC->SetBinContent(i,5*entryAF16);

    //    //printf("%f\t%f\t%f\t%f\t%i\n%f\t%f\n",entryBE11,entryAF11,entryBE16,entryAF16,i,eventNormrun11,eventNormrun16);
    //    //printf("run11bins=%i\n",run11pTf0->GetNbinsX()); 
    //} 

    // this is an important one

    TCanvas *cf0pT = new TCanvas("cf0pT");
    cf0pT->SetLogy();

    run11pTf0->Draw("hist");
    run16pTf0->SetLineColor(2);
    run16pTf0->SetMarkerColor(2);
    run16pTf0->Draw("hist same");
    run16pTf0TOFTPC->SetLineColor(3);
    run16pTf0TOFTPC->SetMarkerColor(3);
    run16pTf0TOFTPC->Draw("hist same");

    cf0pT->SaveAs("mainPlots/kinematics/f0pTcompare.pdf");
    cf0pT->Close();
    
    // end of pt distr for f0


    // charged particle kinematics
    TH1F *run11hppt;
    TH1F *run11hnpt;
    TH1F *run11hpeta;
    TH1F *run11hneta;
    TH1F *run11hpphi;
    TH1F *run11hnphi;
    TH1F *run11hpgDca;
    TH1F *run11hngDca;
    TH1F *run11hpnHitsdEdx;
    TH1F *run11hnnHitsdEdx;
    TH1F *run11hpnHitsFit;
    TH1F *run11hnnHitsFit;
    TH1F *run11hpnHitsRatio;
    TH1F *run11hnnHitsRatio;

    TH1F *run14hppt;
    TH1F *run14hnpt;
    TH1F *run14hpeta;
    TH1F *run14hneta;
    TH1F *run14hpphi;
    TH1F *run14hnphi;
    TH1F *run14hpgDca;
    TH1F *run14hngDca;
    TH1F *run14hpnHitsdEdx;
    TH1F *run14hnnHitsdEdx;
    TH1F *run14hpnHitsFit;
    TH1F *run14hnnHitsFit;
    TH1F *run14hpnHitsRatio;
    TH1F *run14hnnHitsRatio;

    TH1F *run16hppt;
    TH1F *run16hnpt;
    TH1F *run16hpeta;
    TH1F *run16hneta;
    TH1F *run16hpphi;
    TH1F *run16hnphi;
    TH1F *run16hpgDca;
    TH1F *run16hngDca;
    TH1F *run16hpnHitsdEdx;
    TH1F *run16hnnHitsdEdx;
    TH1F *run16hpnHitsFit;
    TH1F *run16hnnHitsFit;
    TH1F *run16hpnHitsRatio;
    TH1F *run16hnnHitsRatio;

    // comparing primary tracks bf and after TOF and PID
    // separate group of cuts, trying to keep isolated
    TH1F *run16primhppt;
    TH1F *run16primhbfppt;
    TH1F *run16primhmidppt;

    TH1F *run16glhppt;
    TH1F *run16glhbfppt;
    TH1F *run16glhmidppt;

    // also, just check ppt at the moment, quite similar to npt already
    TFile *fprimTOF16 = new TFile("DATA/OLDDATA/run16primaryTOF.root");
    TFile *fglobTOF16 = new TFile("DATA/OLDDATA/run16globalTOF.root");

    run16primhppt = (TH1F*)fprimTOF16->Get("ppt");
    run16primhbfppt = (TH1F*)fprimTOF16->Get("bfppt");
    run16primhmidppt = (TH1F*)fprimTOF16->Get("midppt");

    run16glhppt = (TH1F*)fglobTOF16->Get("ppt");
    run16glhbfppt = (TH1F*)fglobTOF16->Get("bfppt");
    run16glhmidppt = (TH1F*)fglobTOF16->Get("midppt");

    // keep an eye out to see if global ppt changes
    
    TCanvas *cprimbfppt = new TCanvas("cprimbfppt");
    cprimbfppt->SetLogy();

    run16glhbfppt->SetLineColor(1);
    run16glhbfppt->Draw();
    run16primhbfppt->SetLineColor(2);
    run16primhbfppt->Draw("same");

    cprimbfppt->SaveAs("mainPlots/kinematics/primvsglobbfppt.pdf");
    cprimbfppt->Close();

    TCanvas *cprimmidppt = new TCanvas("cprimmidppt");
    cprimmidppt->SetLogy();

    run16glhmidppt->SetLineColor(1);
    run16glhmidppt->Draw();
    run16primhmidppt->SetLineColor(2);
    run16primhmidppt->Draw("same");

    cprimmidppt->SaveAs("mainPlots/kinematics/primvsglobmidppt.pdf");
    cprimmidppt->Close();

    TCanvas *cprimafppt = new TCanvas("cprimafppt");
    cprimafppt->SetLogy();

    run16glhppt->SetLineColor(1);
    run16glhppt->Draw();
    run16primhppt->SetLineColor(2);
    run16primhppt->Draw("same");

    cprimafppt->SaveAs("mainPlots/kinematics/primvsglobafppt.pdf");
    cprimafppt->Close();

   
    // RUN COMPARISON CHARGED PARTICLE KINEMATICS

    // end of isolated cuts

    // positive variables
    run11hppt = (TH1F*)f11->Get("ppt");
    run11hpeta = (TH1F*)f11->Get("peta");
    run11hpphi = (TH1F*)f11->Get("pphi"); 
    run11hpnHitsdEdx = (TH1F*)f11->Get("pnHitsdEdx");
    run11hpnHitsFit = (TH1F*)f11->Get("pnHitsFit");
    run11hpnHitsRatio = (TH1F*)f11->Get("pnHitsRatio");
    run11hpgDca = (TH1F*)f11->Get("pgDca");   

    run14hppt = (TH1F*)f14->Get("ppt");
    run14hpeta = (TH1F*)f14->Get("peta");
    run14hpphi = (TH1F*)f14->Get("pphi"); 
    run14hpnHitsdEdx = (TH1F*)f14->Get("pnHitsdEdx");
    run14hpnHitsFit = (TH1F*)f14->Get("pnHitsFit");
    run14hpnHitsRatio = (TH1F*)f14->Get("pnHitsRatio");
    run14hpgDca = (TH1F*)f14->Get("pgDca");   

    run16hppt = (TH1F*)f16->Get("ppt");
    run16hpeta = (TH1F*)f16->Get("peta");
    run16hpphi = (TH1F*)f16->Get("pphi"); 
    run16hpnHitsdEdx = (TH1F*)f16->Get("pnHitsdEdx");
    run16hpnHitsFit = (TH1F*)f16->Get("pnHitsFit");
    run16hpnHitsRatio = (TH1F*)f16->Get("pnHitsRatio");
    run16hpgDca = (TH1F*)f16->Get("pgDca");   

    // negative variables
    run11hnpt = (TH1F*)f11->Get("npt");
    run11hneta = (TH1F*)f11->Get("neta");
    run11hnphi = (TH1F*)f11->Get("nphi"); 
    run11hnnHitsdEdx = (TH1F*)f11->Get("nnHitsdEdx");
    run11hnnHitsFit = (TH1F*)f11->Get("nnHitsFit");
    run11hnnHitsRatio = (TH1F*)f11->Get("nnHitsRatio");
    run11hngDca = (TH1F*)f11->Get("ngDca");   

    run14hnpt = (TH1F*)f14->Get("npt");
    run14hneta = (TH1F*)f14->Get("neta");
    run14hnphi = (TH1F*)f14->Get("nphi"); 
    run14hnnHitsdEdx = (TH1F*)f14->Get("nnHitsdEdx");
    run14hnnHitsFit = (TH1F*)f14->Get("nnHitsFit");
    run14hnnHitsRatio = (TH1F*)f14->Get("nnHitsRatio");
    run14hngDca = (TH1F*)f14->Get("ngDca");   

    run16hnpt = (TH1F*)f16->Get("npt");
    run16hneta = (TH1F*)f16->Get("neta");
    run16hnphi = (TH1F*)f16->Get("nphi"); 
    run16hnnHitsdEdx = (TH1F*)f16->Get("nnHitsdEdx");
    run16hnnHitsFit = (TH1F*)f16->Get("nnHitsFit");
    run16hnnHitsRatio = (TH1F*)f16->Get("nnHitsRatio");
    run16hngDca = (TH1F*)f16->Get("ngDca");   


    double rawBins = 0;
    double scaledBins = 0;
      
    run11hppt->Scale(1/run11TotEvents);
    run11hnpt->Scale(1/run11TotEvents);
    run14hppt->Scale(1/run14TotEvents);
    run14hnpt->Scale(1/run14TotEvents);
    run16hppt->Scale(1/run16TotEvents);
    run16hnpt->Scale(1/run16TotEvents); 

    //for(int i = 0; i < (run11hppt->GetNbinsX()) ; i++)
    //{
    //    rawBins = run11hppt->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hppt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hnpt->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hnpt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hppt->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hppt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hnpt->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hnpt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hppt->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hppt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hnpt->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hnpt->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpeta->Scale(1/run11TotEvents);
    run11hneta->Scale(1/run11TotEvents);
    run14hpeta->Scale(1/run14TotEvents);
    run14hneta->Scale(1/run14TotEvents);
    run16hpeta->Scale(1/run16TotEvents);
    run16hneta->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpeta->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpeta->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpeta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hneta->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hneta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpeta->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpeta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hneta->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hneta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpeta->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpeta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hneta->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hneta->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpphi->Scale(1/run11TotEvents);
    run11hnphi->Scale(1/run11TotEvents);
    run14hpphi->Scale(1/run14TotEvents);
    run14hnphi->Scale(1/run14TotEvents);
    run16hpphi->Scale(1/run16TotEvents);
    run16hnphi->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpphi->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hnphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hnphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hnphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hnphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hnphi->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hnphi->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpgDca->Scale(1/run11TotEvents);
    run11hngDca->Scale(1/run11TotEvents);
    run14hpgDca->Scale(1/run14TotEvents);
    run14hngDca->Scale(1/run14TotEvents);
    run16hpgDca->Scale(1/run16TotEvents);
    run16hngDca->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpgDca->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpgDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpgDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hngDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hngDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpgDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpgDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hngDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hngDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpgDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpgDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hngDca->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hngDca->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpnHitsdEdx->Scale(1/run11TotEvents);
    run11hnnHitsdEdx->Scale(1/run11TotEvents);
    run14hpnHitsdEdx->Scale(1/run14TotEvents);
    run14hnnHitsdEdx->Scale(1/run14TotEvents);
    run16hpnHitsdEdx->Scale(1/run16TotEvents);
    run16hnnHitsdEdx->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpnHitsdEdx->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hnnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hnnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hnnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hnnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hnnHitsdEdx->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hnnHitsdEdx->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpnHitsFit->Scale(1/run11TotEvents);
    run11hnnHitsFit->Scale(1/run11TotEvents);
    run14hpnHitsFit->Scale(1/run14TotEvents);
    run14hnnHitsFit->Scale(1/run14TotEvents);
    run16hpnHitsFit->Scale(1/run16TotEvents);
    run16hnnHitsFit->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpnHitsFit->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hnnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hnnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hnnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hnnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hnnHitsFit->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hnnHitsFit->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    run11hpnHitsRatio->Scale(1/run11TotEvents);
    run11hnnHitsRatio->Scale(1/run11TotEvents);
    run14hpnHitsRatio->Scale(1/run14TotEvents);
    run14hnnHitsRatio->Scale(1/run14TotEvents);
    run16hpnHitsRatio->Scale(1/run16TotEvents);
    run16hnnHitsRatio->Scale(1/run16TotEvents);

    //for(int i = 0; i < ( run11hpnHitsRatio->GetNbinsX() ) ; i++)
    //{
    //    rawBins = run11hpnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hpnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run11hnnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run11TotEvents;
    //    run11hnnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hpnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hpnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run14hnnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run14TotEvents;
    //    run14hnnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hpnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hpnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;

    //    rawBins = run16hnnHitsRatio->GetBinContent(i+1);
    //    scaledBins = rawBins/run16TotEvents;
    //    run16hnnHitsRatio->SetBinContent(i+1,scaledBins);

    //    rawBins = 0;
    //    scaledBins = 0;
    //}

    // plots, will compare just the runs
        
    TCanvas *cppt = new TCanvas("cppt");
    cppt->Divide(2,3);
    cppt->cd(1)->SetLogy();

    run11hppt->Draw();

    run16hppt->SetMarkerColor(2);
    run16hppt->SetLineColor(2);
    run16hppt->Draw("same");

    run14hppt->SetMarkerColor(3);
    run14hppt->SetLineColor(3);
    run14hppt->Draw("same");
   
    TLegend *lppt = new TLegend(0.6,0.25,1.0,0.45);
    lppt->AddEntry(run11hppt,"run11");
    lppt->AddEntry(run14hppt,"run14");
    lppt->AddEntry(run16hppt,"run16");
    lppt->Draw();

    //cppt->SaveAs("mainPlots/kinematics/ppt.pdf");
    //cppt->Close();
    
    //TCanvas *cnpt = new TCanvas("cnpt");
    //cnpt->SetLogy();
    cppt->cd(2)->SetLogy();

    run11hnpt->Draw();

    run16hnpt->SetLineColor(2);
    run16hnpt->Draw("same");

    run14hnpt->SetLineColor(3);
    run14hnpt->Draw("same");

    //cnpt->SaveAs("mainPlots/kinematics/npt.pdf");
    //cnpt->Close();
    //cppt->SaveAs("mainPlots/kinematics/test01.pdf");

    //TCanvas *cpeta = new TCanvas("cpeta");
    cppt->cd(3);

    run11hpeta->Draw();
    run16hpeta->SetLineColor(2);
    run16hpeta->Draw("same");
    run14hpeta->SetLineColor(3);
    run14hpeta->Draw("same");

    //cpeta->SaveAs("mainPlots/kinematics/peta.pdf");
    //cpeta->Close();

    //TCanvas *cneta = new TCanvas("cneta");
    cppt->cd(4);

    run11hneta->Draw();
    run16hneta->SetLineColor(2);
    run16hneta->Draw("same");

    run14hneta->SetLineColor(3);
    run14hneta->Draw("same");
   
    //cneta->SaveAs("mainPlots/kinematics/neta.pdf");
    //cneta->Close();
    
    //TCanvas *cpphi = new TCanvas("cpphi");
    cppt->cd(5);

    run11hpphi->Draw();
    run16hpphi->SetLineColor(2);
    run16hpphi->Draw("same");

    run14hpphi->SetLineColor(3);
    run14hpphi->Draw("same");

    //cpphi->SaveAs("mainPlots/kinematics/pphi.pdf");
    //cpphi->Close();

    //TCanvas *cnphi = new TCanvas("cnphi");
    cppt->cd(6);

    run11hnphi->Draw();

    run14hnphi->SetLineColor(3);
    run14hnphi->Draw("same");

    run16hnphi->SetLineColor(2);
    run16hnphi->Draw("same");

    cppt->SaveAs("mainPlots/kinematics/singPartKinemat.pdf");

    //cnphi->SaveAs("mainPlots/kinematics/nphi.pdf");
    //cnphi->Close();

    TCanvas *cpnHitsdEdx = new TCanvas("cpnHitsdEdx");
    
    run11hpnHitsdEdx->Draw();
    run16hpnHitsdEdx->SetLineColor(2);
    run16hpnHitsdEdx->Draw("same");

    run14hpnHitsdEdx->SetLineColor(3);
    run14hpnHitsdEdx->Draw("same");

    cpnHitsdEdx->SaveAs("mainPlots/kinematics/pnHitsdEdx.pdf");
    cpnHitsdEdx->Close();

    TCanvas *cnnHitsdEdx = new TCanvas("cnnHitsdEdx");
    
    run11hnnHitsdEdx->Draw();
    run16hnnHitsdEdx->SetLineColor(2);
    run16hnnHitsdEdx->Draw("same");

    run14hnnHitsdEdx->SetLineColor(3);
    run14hnnHitsdEdx->Draw("same");

    cnnHitsdEdx->SaveAs("mainPlots/kinematics/nnHitsdEdx.pdf");
    cnnHitsdEdx->Close();

    TCanvas *cpnHitsFit = new TCanvas("cpnHitsFit");
    
    run11hpnHitsFit->Draw();
    run16hpnHitsFit->SetLineColor(2);
    run16hpnHitsFit->Draw("same");

    run14hpnHitsFit->SetLineColor(3);
    run14hpnHitsFit->Draw("same");

    cpnHitsFit->SaveAs("mainPlots/kinematics/pnHitsFit.pdf");
    cpnHitsFit->Close();

    TCanvas *cnnHitsFit = new TCanvas("cnnHitsFit");
    
    run11hnnHitsFit->Draw();
    run16hnnHitsFit->SetLineColor(2);
    run16hnnHitsFit->Draw("same");

    run14hnnHitsFit->SetLineColor(3);
    run14hnnHitsFit->Draw("same");

    cnnHitsFit->SaveAs("mainPlots/kinematics/nnHitsFit.pdf");
    cnnHitsFit->Close();

    TCanvas *cpnHitsRatio = new TCanvas("cpnHitsRatio");
    
    run11hpnHitsRatio->Draw();
    run16hpnHitsRatio->SetLineColor(2);
    run16hpnHitsRatio->Draw("same");

    run14hpnHitsRatio->SetLineColor(3);
    run14hpnHitsRatio->Draw("same");

    cpnHitsRatio->SaveAs("mainPlots/kinematics/pnHitsRatio.pdf");
    cpnHitsRatio->Close();

    TCanvas *cnnHitsRatio = new TCanvas("cnnHitsRatio");
    
    run11hnnHitsRatio->Draw();
    run16hnnHitsRatio->SetLineColor(2);
    run16hnnHitsRatio->Draw("same");

    run14hnnHitsRatio->SetLineColor(3);
    run14hnnHitsRatio->Draw("same");

    cnnHitsRatio->SaveAs("mainPlots/kinematics/nnHitsRatio.pdf");
    cnnHitsRatio->Close();

    TCanvas *cpgDca = new TCanvas("cpgDca");
    
    run16hpgDca->SetLineColor(2);
    run16hpgDca->Draw();
    run11hpgDca->Draw("same");

    run14hpgDca->SetLineColor(3);
    run14hpgDca->Draw("same");

    cpgDca->SaveAs("mainPlots/kinematics/pgDca.pdf");
    cpgDca->Close();

    TCanvas *cngDca = new TCanvas("cngDCa");

    run16hngDca->SetLineColor(2);
    run16hngDca->Draw();
    run11hngDca->Draw("same");

    run14hngDca->SetLineColor(3);
    run14hngDca->Draw("same");

    cngDca->SaveAs("mainPlots/kinematics/ngDca.pdf");
    cngDca->Close();

    return;
}

// this function examines the f0 peak for different rootfiles
void primaryFit::f0pTtests()
{
    int centBins=9;
    int myCents=6;

    vector<TH1D*> histPush;

    // these names for f11 are misleading. Fix Them!!
    // next data set
    TFile *f11DataOld = new TFile("DATA/OLDDATA/run11Data.root"); // full Data set
    //TFile *f11DataOld = new TFile("DATA/Run11_set01/run11Data_set01.root"); // 100mil events 
    TH3D *run11ULOld3d[centBins];
    TH3D *run11LSOld3d[centBins];

    vector < vector<TH1D*> >run11OldX;
    TH1D *run11OldY;

    double run11OldCentEvents;
    double run11OldTotEvents;

    // next data set
    //TFile *f16DataOld = new TFile("DATA/OLDDATA/run16Data.root"); // full Data set, does not have flat event plane
    //TH3D *run16ULOld3d[centBins];
    //TH3D *run16LSOld3d[centBins];
   
    //vector < vector<TH1D*> > run16OldX;
    //TH1D *run16OldY;

    //double run16OldCentEvents;
    //double run16OldTotEvents;

    // next data set
    TFile *f16noTOFglob = new TFile("DATA/OLDDATA/run16noTOFglob99job.root"); // only 99 jobs out of 5000
    TH3D *run16ULnoTOFglob3d[centBins];
    TH3D *run16LSnoTOFglob3d[centBins];

    vector < vector<TH1D*> > run16noTOFglobX;
    TH1D *run16noTOFglobY;

    double run16noTOFglobCentEvents;
    double run16noTOFglobTotEvents;
    
    // next data set
    //TFile *f16flatevent = new TFile("DATA/OLDDATA/run16eventtest_03.root"); // not sure how many data files in this, but I can always get the events
    TFile *f16flatevent = new TFile("DATA/Run16_full_correct/Run16data200.root"); // not sure how many data files in this, but I can always get the events
    TH3D *run16ULflatevent3d[centBins];
    TH3D *run16LSflatevent3d[centBins];
   
    vector < vector<TH1D*> > run16flatX;
    TH1D *run16flatY;

    double run16flatCentEvents;
    double run16flatTotEvents;

    // next data set
    TFile *f16_35per = new TFile("DATA/OLDDATA/run16_35per.root");
    TH3D *run16UL35per3d[centBins];
    TH3D *run16LS35per3d[centBins];

    vector < vector<TH1D*> > run16_35perX;
    TH1D *run16_35perY;

    double run16_35perCentEvents;
    double run16_35perTotEvents;
    
    // next data set
    //TFile *f16_50per = new TFile("DATA/Run16_notof_half/run16_50perData.root");
    TFile *f16_50per = new TFile("DATA/Run16_set01/run16Data_set01.root");
    TH3D *run16UL50per3d[centBins];
    TH3D *run16LS50per3d[centBins];

    vector < vector<TH1D*> > run16_50perX;
    TH1D *run16_50perY;

    double run16_50perCentEvents;
    double run16_50perTotEvents;

    // next data set
    TFile *f16_momCorrPart = new TFile("DATA/Run16_momCorr/run16_momCorrect.root");
    TH3D *run16ULmomCorrPart3d[centBins];
    TH3D *run16LSmomCorrPart3d[centBins];
    
    vector < vector<TH1D*> > run16_momCorrPartX;
    TH1D *run16_momCorrPartY;

    double run16_momCorrPartCentEvents;
    double run16_momCorrPartTotEvents;

    // next data set
    TFile *f16_momCorr02Part = new TFile("DATA/Run16_momCorr_02/run16_momCorrect_02.root");
    TH3D *run16ULmomCorr02Part3d[centBins];
    TH3D *run16LSmomCorr02Part3d[centBins];

    vector < vector<TH1D*> > run16_momCorr02PartX;
    TH1D *run16_momCorr02PartY;
    TH1F *run16_momCorr02PartCent; // changing where I place this
    
    double run16_momCorr02PartCentEvents;
    double run16_momCorr02PartTotEvents;

    // next data set
    //TFile *f10_first = new TFile("DATA/Run10_set01/run10Data.root");
    //TFile *f10_first = new TFile("DATA/Run14_noHFT_03/run14data200.root");
    TFile *f10_first = new TFile("DATA/Run14_test02/run14data200.root");
    TH3D *run10ULfirst3d[centBins];
    TH3D *run10LSfirst3d[centBins];

    vector < vector<TH1D*> > run10firstX;
    TH1D *run10firstY;
    TH1F *run10firstCent;

    double run10firstCentEvents;
    double run10firstTotEvents;

    TFile *f14_sec = new TFile("DATA/Run14_noHFT_03/run14data200.root");
    
    TH3D *run14UL3d[centBins];
    TH3D *run14LS3d[centBins];

    vector < vector<TH1D*> > run14X;
    TH1D *run14Y;
    TH1F *run14Cent;

    double run14CentEvents;
    double run14TotEvents;


    // END OF DATA SET DECLARATIONS

    char ULname[256];
    char LSname[256];

    for(int i = 0; i < centBins; i++)
    {
        sprintf(ULname,"TMdphiUL0c%i",i+1);         
        sprintf(LSname,"TMdphiLS0c%i",i+1);

        run11ULOld3d[i] = (TH3D*) f11DataOld->Get(ULname);    
        run11LSOld3d[i] = (TH3D*) f11DataOld->Get(LSname);   

        //run16ULOld3d[i] = (TH3D*) f16DataOld->Get(ULname);
        //run16LSOld3d[i] = (TH3D*) f16DataOld->Get(LSname);

        run16ULnoTOFglob3d[i] = (TH3D*) f16noTOFglob->Get(ULname); 
        run16LSnoTOFglob3d[i] = (TH3D*) f16noTOFglob->Get(LSname);

        run16ULflatevent3d[i] = (TH3D*) f16flatevent->Get(ULname);
        run16LSflatevent3d[i] = (TH3D*) f16flatevent->Get(LSname);

        run16UL35per3d[i] = (TH3D*) f16_35per->Get(ULname);
        run16LS35per3d[i] = (TH3D*) f16_35per->Get(LSname);

        run16UL50per3d[i] = (TH3D*) f16_50per->Get(ULname);
        run16LS50per3d[i] = (TH3D*) f16_50per->Get(LSname);

        run16ULmomCorrPart3d[i] = (TH3D*) f16_momCorrPart->Get(ULname);
        run16LSmomCorrPart3d[i] = (TH3D*) f16_momCorrPart->Get(LSname);

        run16ULmomCorr02Part3d[i] = (TH3D*) f16_momCorr02Part->Get(ULname);
        run16LSmomCorr02Part3d[i] = (TH3D*) f16_momCorr02Part->Get(LSname);

        run10ULfirst3d[i] = (TH3D*) f10_first->Get(ULname);
        run10LSfirst3d[i] = (TH3D*) f10_first->Get(LSname);
    }

    // adding centralities
    for(int i = 1; i< myCents; i++)
    {
        run11ULOld3d[0]->Add(run11ULOld3d[i],1);
        run11LSOld3d[0]->Add(run11LSOld3d[i],1);
        
        //run16ULOld3d[0]->Add(run16ULOld3d[i],1);
        //run16LSOld3d[0]->Add(run16LSOld3d[i],1);
        
        run16ULnoTOFglob3d[0]->Add(run16ULnoTOFglob3d[i],1);
        run16LSnoTOFglob3d[0]->Add(run16LSnoTOFglob3d[i],1);
        
        run16ULflatevent3d[0]->Add(run16ULflatevent3d[i],1);
        run16LSflatevent3d[0]->Add(run16LSflatevent3d[i],1);

        run16UL35per3d[0]->Add(run16UL35per3d[i],1);
        run16LS35per3d[0]->Add(run16LS35per3d[i],1);

        run16UL50per3d[0]->Add(run16UL50per3d[i],1);
        run16LS50per3d[0]->Add(run16LS50per3d[i],1);

        run16ULmomCorrPart3d[0]->Add(run16ULmomCorrPart3d[i],1);
        run16LSmomCorrPart3d[0]->Add(run16LSmomCorrPart3d[i],1);
    
        run16ULmomCorr02Part3d[0]->Add(run16ULmomCorr02Part3d[i],1);
        run16LSmomCorr02Part3d[0]->Add(run16LSmomCorr02Part3d[i],1);

        run10ULfirst3d[0]->Add(run10ULfirst3d[i],1);
        run10LSfirst3d[0]->Add(run10LSfirst3d[i],1);
    }

    // subtracting ul and ls
    run11ULOld3d[0]->Add(run11LSOld3d[0],-1);
    //run16ULOld3d[0]->Add(run16LSOld3d[0],-1);
    run16ULnoTOFglob3d[0]->Add(run16LSnoTOFglob3d[0],-1);
    run16ULflatevent3d[0]->Add(run16LSflatevent3d[0],-1);
    run16UL35per3d[0]->Add(run16LS35per3d[0],-1);
    run16UL50per3d[0]->Add(run16LS50per3d[0],-1);
    run16ULmomCorrPart3d[0]->Add(run16LSmomCorrPart3d[0],-1);
    run16ULmomCorr02Part3d[0]->Add(run16LSmomCorr02Part3d[0],-1);
    run10ULfirst3d[0]->Add(run10LSfirst3d[0],-1);

    TCanvas *tmpC = new TCanvas();

    run16UL50per3d[0]->ProjectionX()->Draw();

    tmpC->SaveAs("tmp.pdf");
    tmpC->Close();

    // I should move these up with the rest of the data sets 
    // getting the number of events
    TH1F *run11OldCent;
    //TH1F *run16OldCent;
    TH1F *run16noTOFglobCent;
    TH1F *run16flatCent;
    TH1F *run16_35perCent;
    TH1F *run16_50perCent;
    TH1F *run16_momCorrPartCent;
    
    
    run11OldCent = (TH1F*) f11DataOld->Get("centralityw");
    //run16OldCent = (TH1F*) f16DataOld->Get("centralityw");
    run16noTOFglobCent = (TH1F*) f16noTOFglob->Get("centralityw");
    run16flatCent = (TH1F*) f16flatevent->Get("centralityw");
    run16_35perCent = (TH1F*) f16_35per->Get("centralityw");
    run16_50perCent = (TH1F*) f16_50per->Get("centralityw");
    run16_momCorrPartCent = (TH1F*) f16_momCorrPart->Get("centralityw");
    run16_momCorr02PartCent = (TH1F*) f16_momCorr02Part->Get("centralityw");
    run10firstCent = (TH1F*) f10_first->Get("centralityw");

    run11OldTotEvents = run11OldCent->GetBinContent(1);
    run11OldCentEvents = run11OldCent->Integral(2,7);

    //run16OldTotEvents = run16OldCent->GetBinContent(1);
    //run16OldCentEvents = run16OldCent->Integral(2,7);

    run16noTOFglobTotEvents = run16noTOFglobCent->GetBinContent(1);
    run16noTOFglobCentEvents = run16noTOFglobCent->Integral(2,7); 

    run16flatTotEvents = run16flatCent->GetBinContent(1);
    run16flatCentEvents = run16flatCent->Integral(2,7);

    run16_35perTotEvents = run16_35perCent->GetBinContent(1);
    run16_35perCentEvents = run16_35perCent->Integral(2,7);

    run16_50perTotEvents = run16_50perCent->GetBinContent(1);
    run16_50perCentEvents = run16_50perCent->Integral(2,7);

    run16_momCorrPartTotEvents = run16_momCorrPartCent->GetBinContent(1);
    run16_momCorrPartCentEvents = run16_momCorrPartCent->Integral(2,7);

    run16_momCorr02PartTotEvents = run16_momCorr02PartCent->GetBinContent(1);
    run16_momCorr02PartCentEvents = run16_momCorr02PartCent->Integral(2,7);

    run10firstTotEvents = run10firstCent->GetBinContent(1);
    run10firstCentEvents = run10firstCent->Integral(2,7);

    // getting pT distribution with specific mass cuts

    run11OldY = run11ULOld3d[0]->ProjectionY("run11OldpT",180,220,1,20);
    //run16OldY = run16ULOld3d[0]->ProjectionY("run16OldpT",180,220,1,20);
    run16noTOFglobY = run16ULnoTOFglob3d[0]->ProjectionY("run16noTOFglobpT",180,220,1,20);
    run16flatY = run16ULflatevent3d[0]->ProjectionY("run16flateventpT",180,220,1,20);
    run16_35perY = run16UL35per3d[0]->ProjectionY("run16_35percent",180,220,1,20);
    run16_50perY = run16UL50per3d[0]->ProjectionY("run16_50percent",180,220,1,20);
    run16_momCorrPartY = run16ULmomCorrPart3d[0]->ProjectionY("run16_momCorrectionPartial",180,220,1,20);
    run16_momCorr02PartY = run16ULmomCorr02Part3d[0]->ProjectionY("run16 with TOF and MomCorr",180,220,1,20);
    run10firstY = run10ULfirst3d[0]->ProjectionY("run10 with TOF",180,220,1,20);

    char PTprojNames[256];
    // pt projections, over all phi
    
    run11OldX.push_back(histPush);
    //run16OldX.push_back(histPush);
    run16noTOFglobX.push_back(histPush);
    run16flatX.push_back(histPush);
    run16_35perX.push_back(histPush);
    run16_50perX.push_back(histPush);
    run16_momCorrPartX.push_back(histPush);
    run16_momCorr02PartX.push_back(histPush);
    run10firstX.push_back(histPush);
   
    for(int i = 0; i < 6;i++)
    {
        sprintf(PTprojNames,"11old%i",i);
        run11OldX[0].push_back( run11ULOld3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );    
       
        //sprintf(PTprojNames,"16old%i",i); 
        //run16OldX[0].push_back( run16ULOld3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );    
        
        sprintf(PTprojNames,"16noTOFglob%i",i);
        run16noTOFglobX[0].push_back( run16ULnoTOFglob3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );    
        
        sprintf(PTprojNames,"16flat%i",i);
        run16flatX[0].push_back( run16ULflatevent3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );   

        sprintf(PTprojNames,"16_35per%i",i);
        run16_35perX[0].push_back( run16UL35per3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );

        sprintf(PTprojNames,"16_50per%i",i);
        run16_50perX[0].push_back( run16UL50per3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );

        sprintf(PTprojNames,"16_momCorrPart%i",i);
        run16_momCorrPartX[0].push_back( run16ULmomCorrPart3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );

        sprintf(PTprojNames,"16_momCorr02Part%i",i);
        run16_momCorr02PartX[0].push_back( run16ULmomCorr02Part3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );


        sprintf(PTprojNames,"10TOF%i",i);
        run10firstX[0].push_back( run10ULfirst3d[0]->ProjectionX(PTprojNames,ptDisBin[i][0],ptDisBin[i][1],ptDisBin[i][2],ptDisBin[i][3]) );

    }

    // scaling mass distributions
    
    //cout << run11ULOld3d[0]->GetNbinsX() << endl;

    double oldEntry = 0;
    double aftEntry = 0;

    // can use
    // hist->Scale(numEvents);

    for(int j = 0; j < run11OldX[0].size(); j++)
    {
        for(int i = 1; i < (run11ULOld3d[0]->GetNbinsX() + 1 ) ; i++)
        {
            oldEntry = run11OldX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run11OldCentEvents;
            run11OldX[0][j]->SetBinContent(i,aftEntry);

            //oldEntry = run16OldX[0][j]->GetBinContent(i);
            //aftEntry = oldEntry/run16OldCentEvents;
            //run16OldX[0][j]->SetBinContent(i,aftEntry);

            //cout << aftEntry << "\t" << i << endl;
            //cout << run16OldCentEvents << endl;;
            
            oldEntry = run16noTOFglobX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16noTOFglobCentEvents;
            run16noTOFglobX[0][j]->SetBinContent(i,aftEntry);

            oldEntry = run16flatX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16flatCentEvents;
            run16flatX[0][j]->SetBinContent(i,aftEntry);

            oldEntry = run16_35perX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16_35perCentEvents;
            run16_35perX[0][j]->SetBinContent(i,aftEntry);

            oldEntry = run16_50perX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16_50perCentEvents;
            run16_50perX[0][j]->SetBinContent(i,aftEntry);
 
            oldEntry = run16_momCorrPartX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16_momCorrPartCentEvents;
            run16_momCorrPartX[0][j]->SetBinContent(i,aftEntry);
 
            oldEntry = run16_momCorr02PartX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run16_momCorr02PartCentEvents;
            run16_momCorr02PartX[0][j]->SetBinContent(i,aftEntry);
 
            oldEntry = run10firstX[0][j]->GetBinContent(i);
            aftEntry = oldEntry/run10firstCentEvents;
            run10firstX[0][j]->SetBinContent(i,aftEntry);
            
        } 
    }

    // scaling pT distribution 
    
    //cout << run11ULOld3d[0]->GetNbinsY() << endl;

    run11OldY->Scale(1/run11OldCentEvents);
    run16noTOFglobY->Scale(1/run16noTOFglobCentEvents);
    run16flatY->Scale(1/run16flatCentEvents);
    run16_35perY->Scale(1/run16_35perCentEvents);
    run16_50perY->Scale(1/run16_50perCentEvents);
    run16_momCorrPartY->Scale(1/run16_momCorrPartCentEvents);
    run16_momCorr02PartY->Scale(1/run16_momCorr02PartCentEvents);
    run10firstY->Scale(1/run10firstCentEvents);

    //for(int i = 1; i < run11ULOld3d[0]->GetNbinsY() +1;i++)
    //{
    //    oldEntry = run11OldY->GetBinContent(i);
    //    aftEntry = oldEntry/run11OldCentEvents;
    //    run11OldY->SetBinContent(i,aftEntry);

    //    //oldEntry = run16OldY->GetBinContent(i);
    //    //aftEntry = oldEntry/run16OldCentEvents;
    //    //run16OldY->SetBinContent(i,aftEntry);

    //    oldEntry = run16noTOFglobY->GetBinContent(i);
    //    aftEntry = oldEntry/run16noTOFglobCentEvents;
    //    run16noTOFglobY->SetBinContent(i,aftEntry);

    //    oldEntry = run16flatY->GetBinContent(i);
    //    aftEntry = oldEntry/run16flatCentEvents;
    //    run16flatY->SetBinContent(i,aftEntry);

    //    oldEntry = run16_35perY->GetBinContent(i);
    //    aftEntry = oldEntry/run16_35perCentEvents;
    //    run16_35perY->SetBinContent(i,aftEntry);

    //    oldEntry = run16_50perY->GetBinContent(i);
    //    aftEntry = oldEntry/run16_50perCentEvents;
    //    run16_50perY->SetBinContent(i,aftEntry);

    //    oldEntry = run16_momCorrPartY->GetBinContent(i);
    //    aftEntry = oldEntry/run16_momCorrPartCentEvents;
    //    run16_momCorrPartY->SetBinContent(i,aftEntry);

    //    oldEntry = run16_momCorr02PartY->GetBinContent(i);
    //    aftEntry = oldEntry/run16_momCorr02PartCentEvents;
    //    run16_momCorr02PartY->SetBinContent(i,aftEntry);

    //    oldEntry = run10firstY->GetBinContent(i);
    //    aftEntry = oldEntry/run10firstCentEvents;
    //    run10firstY->SetBinContent(i,aftEntry);

    //}


    // Drawing the projections


    // f0 normalized mass distribution
    TCanvas *cmasspT = new TCanvas();
    cmasspT->Divide(3,2);
    TLegend *pTlegend[6];

    char smallTitlesPt[256];
    double localMax;
    double localMin;
   
    // kshort 
    //double xMax=0.6;
    //double xMin=0.3;

    //f0
    double xMax = 1.2;
    double xMin = 0.8;

    for(int i = 0; i < run11OldX[0].size() ; i++)
    {
        cmasspT->cd(i+1);
        run11OldX[0][i]->Draw("hist");
        run11OldX[0][i]->SetAxisRange(xMin,xMax,"X");
        
        localMax = run11OldX[0][i]->GetBinContent( run11OldX[0][i]->GetMaximumBin() );
        localMin = run11OldX[0][i]->GetBinContent( run11OldX[0][i]->GetMinimumBin() );

        sprintf(smallTitlesPt,"Mass Distribution with %.1f < pT < %.1f GeV ", (ptDisBin[i][0]*0.25 - 0.25) , (ptDisBin[i][1]*0.25 ) );

        run11OldX[0][i]->SetTitle(smallTitlesPt);
        run11OldX[0][i]->GetXaxis()->SetTitle("Mass GeV");
        run11OldX[0][i]->GetYaxis()->SetTitle("Counts/Events"); 
        // print out the number of events in each data set

        //run16OldX[0][i]->Draw("same hist");
        //run16OldX[0][i]->SetMarkerColor(2);
        //run16OldX[0][i]->SetLineColor(2);
        //run16OldX[0][i]->SetAxisRange(xMin,xMax,"X");

        //run16noTOFglobX[0][i]->Draw("same hist");
        //run16noTOFglobX[0][i]->SetMarkerColor(3);
        //run16noTOFglobX[0][i]->SetLineColor(3);
        //run16noTOFglobX[0][i]->SetAxisRange(xMin,xMax,"X");

        run16flatX[0][i]->Draw("same hist");
        run16flatX[0][i]->SetMarkerColor(3);
        run16flatX[0][i]->SetLineColor(3);
        run16flatX[0][i]->SetAxisRange(xMin,xMax,"X");

        //run16_35perX[0][i]->Draw("same hist");
        run16_35perX[0][i]->SetMarkerColor(2);
        run16_35perX[0][i]->SetLineColor(2);
        run16_35perX[0][i]->SetAxisRange(xMin,xMax,"X");

        //run16_50perX[0][i]->Draw("same hist");
        run16_50perX[0][i]->SetMarkerColor(3);
        run16_50perX[0][i]->SetLineColor(3);
        run16_50perX[0][i]->SetAxisRange(xMin,xMax,"X");

        //run16_momCorrPartX[0][i]->Draw("same hist");
        run16_momCorrPartX[0][i]->SetMarkerColor(5);
        run16_momCorrPartX[0][i]->SetLineColor(5);
        run16_momCorrPartX[0][i]->SetAxisRange(xMin,xMax,"X");

        //run16_momCorr02PartX[0][i]->Draw("same hist");
        run16_momCorr02PartX[0][i]->SetMarkerColor(2);
        run16_momCorr02PartX[0][i]->SetLineColor(2);
        run16_momCorr02PartX[0][i]->SetAxisRange(xMin,xMax,"X");

        run10firstX[0][i]->Draw("same hist");
        run10firstX[0][i]->SetMarkerColor(4);
        run10firstX[0][i]->SetLineColor(4);
        run10firstX[0][i]->SetAxisRange(xMin,xMax,"X");


        //if(localMax < run16OldX[0][i]->GetBinContent( run16OldX[0][i]->GetMaximumBin() ) ) 
        //{
        //    localMax = run16OldX[0][i]->GetBinContent( run16OldX[0][i]->GetMaximumBin() );
        //}

        //if(localMax < run16noTOFglobX[0][i]->GetBinContent( run16noTOFglobX[0][i]->GetMaximumBin() ) ) 
        //{
        //    localMax = run16noTOFglobX[0][i]->GetBinContent( run16noTOFglobX[0][i]->GetMaximumBin() );
        //}

        if(localMax < run16flatX[0][i]->GetBinContent( run16flatX[0][i]->GetMaximumBin() ) ) 
        {
            localMax = run16flatX[0][i]->GetBinContent( run16flatX[0][i]->GetMaximumBin() );
        }

        if(localMax < run16_35perX[0][i]->GetBinContent( run16_35perX[0][i]->GetMaximumBin() ) ) 
        {
            localMax = run16_35perX[0][i]->GetBinContent( run16_35perX[0][i]->GetMaximumBin() );
        }

        if(localMax < run16_momCorr02PartX[0][i]->GetBinContent( run16_momCorr02PartX[0][i]->GetMaximumBin() ) ) 
        {
            localMax = run16_momCorr02PartX[0][i]->GetBinContent( run16_momCorr02PartX[0][i]->GetMaximumBin() );
        }

        if(localMax < run10firstX[0][i]->GetBinContent( run10firstX[0][i]->GetMaximumBin() ) ) 
        {
            localMax = run10firstX[0][i]->GetBinContent( run10firstX[0][i]->GetMaximumBin() );
        }


        //if(localMin > run16OldX[0][i]->GetBinContent( run16OldX[0][i]->GetMinimumBin() ) ) 
        //{
        //    localMin = run16OldX[0][i]->GetBinContent( run16OldX[0][i]->GetMinimumBin() );
        //}

        //if(localMin > run16noTOFglobX[0][i]->GetBinContent( run16noTOFglobX[0][i]->GetMinimumBin() ) && i!=5) 
        //{
        //    localMax = run16noTOFglobX[0][i]->GetBinContent( run16noTOFglobX[0][i]->GetMinimumBin() );
        //}

        if(localMin > run16flatX[0][i]->GetBinContent( run16flatX[0][i]->GetMinimumBin() ) ) 
        {
            localMin = run16flatX[0][i]->GetBinContent( run16flatX[0][i]->GetMinimumBin() );
        }

        if(localMin > run16_35perX[0][i]->GetBinContent( run16_35perX[0][i]->GetMinimumBin() ) ) 
        {
            localMin = run16_35perX[0][i]->GetBinContent( run16_35perX[0][i]->GetMinimumBin() );
        }

        if(localMax < run16_momCorr02PartX[0][i]->GetBinContent( run16_momCorr02PartX[0][i]->GetMinimumBin() ) ) 
        {
            localMax = run16_momCorr02PartX[0][i]->GetBinContent( run16_momCorr02PartX[0][i]->GetMinimumBin() );
        }

        if(localMax < run10firstX[0][i]->GetBinContent( run10firstX[0][i]->GetMinimumBin() ) ) 
        {
            localMax = run10firstX[0][i]->GetBinContent( run10firstX[0][i]->GetMinimumBin() );
        }



        // set the y axis range to be based on the highest and lowest values coming from any of plots
        run11OldX[0][i]->SetAxisRange(localMin,localMax,"Y");


        // The legend
        pTlegend[i] = new TLegend(0.6,0.25,1.0,0.45);
        pTlegend[i]->AddEntry(run11OldX[0][i],"run 11");
        //pTlegend[i]->AddEntry(run16OldX[0][i],"run 16 TOF cut, not flat event plane");
        pTlegend[i]->AddEntry(run16flatX[0][i],"run 16");
        //pTlegend[i]->AddEntry(run16noTOFglobX[0][i],"run 16 no TOF cut");
        //pTlegend[i]->AddEntry(run16_35perX[0][i],"run 16 no TOF cut, 35 percent of data");
        //pTlegend[i]->AddEntry(run16_50perX[0][i],"run 16 no TOF cut, 50 percent of data");
        //pTlegend[i]->AddEntry(run16_momCorrPartX[0][i],"run 16 no TOF cut momentum correction");
        //pTlegend[i]->AddEntry(run16_momCorr02PartX[0][i],"run 16 TOF cut momentum correction");
        //pTlegend[i]->AddEntry(run10firstX[0][i],"run 10 TOF");
        pTlegend[i]->AddEntry(run10firstX[0][i],"run 14");

        pTlegend[i]->Draw();
        
    }
    
    cmasspT->SaveAs("mainPlots/f0testPlots/f0MassDistrs.pdf");
    cmasspT->Close();

    // CAN CHANGE TO RATIO HERE!!

    // f0pT distribution here
    TCanvas *cpTdistf0Reg = new TCanvas();
    cpTdistf0Reg->SetLogy();
    TLegend *f0pTleg;

    run11OldY->SetTitle("Number of entries in f0 peak region vs pT");
    run11OldY->GetXaxis()->SetTitle("pT GeV");
    run11OldY->GetYaxis()->SetTitle("Counts/Events");

    run11OldY->Draw("hist");

    //run16OldY->Draw("same hist");
    //run16OldY->SetMarkerColor(2);
    //run16OldY->SetLineColor(2);

    //run16noTOFglobY->Draw("same hist");
    run16noTOFglobY->SetMarkerColor(3);
    run16noTOFglobY->SetLineColor(3);

    run16flatY->Draw("same hist");
    run16flatY->SetMarkerColor(3);
    run16flatY->SetLineColor(3);
 
    //run16_35perY->Draw("same hist");
    run16_35perY->SetMarkerColor(2);
    run16_35perY->SetLineColor(2);
 
    //run16_50perY->Draw("same hist");
    run16_50perY->SetMarkerColor(3);
    run16_50perY->SetLineColor(3);
 
    //run16_momCorrPartY->Draw("same hist");
    run16_momCorrPartY->SetMarkerColor(5);
    run16_momCorrPartY->SetLineColor(5);
 
    //run16_momCorr02PartY->Draw("same hist");
    run16_momCorr02PartY->SetMarkerColor(2);
    run16_momCorr02PartY->SetLineColor(2);
 
    run10firstY->Draw("same hist");
    run10firstY->SetMarkerColor(4);
    run10firstY->SetLineColor(4);
   
    // the legend 
    f0pTleg = new TLegend(0.71,0.25,1.0,0.45);
    f0pTleg->AddEntry(run11OldY,"run 11");
    //f0pTleg->AddEntry(run16OldY,"run 16 TOF cut, uncorrected event plane");
    //f0pTleg->AddEntry(run16noTOFglobY,"run 16 no TOF cut");
    f0pTleg->AddEntry(run16flatY,"run 16");
    //f0pTleg->AddEntry(run16_35perY,"run 16 no TOF cut, 35 percent data");
    //f0pTleg->AddEntry(run16_50perY,"run 16 no TOF cut, 50 percent data");
    //f0pTleg->AddEntry(run16_momCorrPartY,"run 16 no TOF cut momentum correction");
    //f0pTleg->AddEntry(run16_momCorr02PartY,"run 16 TOF cut momentum correction");
    //f0pTleg->AddEntry(run10firstY,"run 10 TOF");
    f0pTleg->AddEntry(run10firstY,"run 14");

    f0pTleg->Draw();

    cpTdistf0Reg->SaveAs("mainPlots/f0testPlots/f0pTdist.pdf");
    cpTdistf0Reg->Close();

    // f0pTdistr ratios
    
    double ratNum = 0;
    
    for(int i = 1; i < run11ULOld3d[0]->GetNbinsY() +1;i++)
    {
        //oldEntry = run16flatY->GetBinContent(i);
        oldEntry = run11OldY->GetBinContent(i); // run10 is the run 14 data at the moment
        ratNum = run10firstY->GetBinContent(i)/run11OldY->GetBinContent(i);
        cout << ratNum << endl;
        if(oldEntry!=0) aftEntry = run10firstY->GetBinContent(i)/oldEntry;
        else aftEntry = 0;
        run11OldY->SetBinContent(i,aftEntry);
    }

    TCanvas *c_f0_11_14Rat = new TCanvas();
    run11OldY->Draw("hist");
    run11OldY->SetTitle("run14 counts divided by run11 counts vs p_{T}");
    c_f0_11_14Rat->SaveAs("mainPlots/f0testPlots/f011and14rat.pdf");
    c_f0_11_14Rat->Close();

    return;
}

// this function is for a momentum correction for TPC by TOF
// due to the suspicion that the HFT is lowering our momentum
// this function is now depracated and serves no purpose
void primaryFit::correctionFit()
{
    TFile *f16_corr = new TFile("DATA/Run16_momCorr_03/run16_momCorrect04.root");
    TH2D *toftpcP;
    TProfile *avgPplt;
    TF1 *fnegExp = new TF1("expFit","[0] - [1]*TMath::Exp(-[2]*TMath::Power(x,[3]))",0,10);

    TF1 *finvPoly1 = new TF1("invPol1","[3]-[0]/([1]+[2]*pow(x,[4]))",0,10);

    TF1 *flog = new TF1("logFit","[0] + [1]*TMath::Log([2]*x + [3])",0,10);


    toftpcP = (TH2D*) f16_corr->Get("pcompareBF");
    avgPplt = toftpcP->ProfileX();

    fnegExp->SetParLimits(0,1.15,1.25);
    fnegExp->SetParLimits(1,0.8,1.1);
    fnegExp->SetParLimits(2,1.2,1.7);
    fnegExp->SetParLimits(3,1.0,2.0);

    finvPoly1->SetParLimits(0,1,50);
    finvPoly1->SetParLimits(1,1e-5,50);
    finvPoly1->SetParLimits(2,1e-5,50);
    finvPoly1->SetParLimits(3,0.9,1.3);
    finvPoly1->SetParLimits(4,3.0,5);

    flog->SetParLimits(2,1e-5,500);
    flog->SetParLimits(3,1e-5,500);

    avgPplt->Fit(fnegExp,"R");    
    avgPplt->Fit(finvPoly1,"R");
    //avgPplt->Fit(flog,"R");

    fnegExp->SetLineColor(2);
    finvPoly1->SetLineColor(3);
    //flog->SetLineColor(4);

    TCanvas *cCorr = new TCanvas();

    avgPplt->Draw();
    fnegExp->Draw("same");
    finvPoly1->Draw("same");
    //flog->Draw("same");

    TLegend *leg01 = new TLegend(0.71,0.25,1.0,0.45);

    cCorr->SaveAs("mainPlots/correctionTests/run16before.pdf");
    cCorr->Close();

    TF1 *fline = new TF1("line","x",0,10);

    double oldEntry=0;
    double newEntry=0;

    for(int i=1;i<avgPplt->GetNbinsX()+1;i++)
    {
        oldEntry = avgPplt->GetBinContent(i);
        if(oldEntry>0)
        {
        newEntry = oldEntry/(fnegExp->Eval(avgPplt->GetBinCenter(i)));
        avgPplt->SetBinContent(i,newEntry);
        //printf("oldvalue=%f\tnewvalue=%f\nbinValue=%i\tbinCenter=%f\n",oldEntry,newEntry,i,avgPplt->GetBinCenter(i));
        }
    }

    TCanvas *cFin = new TCanvas();

    avgPplt->SetAxisRange(-0.5,10,"Y");
    avgPplt->SetMarkerSize(2);

    avgPplt->Draw();
    fline->SetLineColor(2);
    fnegExp->SetLineColor(3);
    fline->Draw("same");
    fnegExp->Draw("same");

    cFin->SaveAs("mainPlots/correctionTests/run16after.pdf");
    cFin->Close();

    return;
}
