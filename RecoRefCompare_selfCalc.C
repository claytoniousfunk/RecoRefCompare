









double JES_fit(double *x, double *par){
    double xx = x[0];
    double fitval = par[0]*TMath::Exp(-1.0*par[1]*x[0]) + par[2]*TMath::Exp(-1.0*par[3]*x[0]) + par[4]*TMath::Exp(-1.0*par[5]*x[0])+par[6];
    return fitval;
}

double fitfxn(double *x, double *par){
    double pt = x[0];
    //int bin = r->FindBin(pt);
    //double fitval = par[0] + par[1]*r->GetBinContent(bin) + par[2]*TMath::Exp(-1*par[3]*r->GetBinContent(bin));
    //double fitval = par[0] + par[1]*pt + par[2]*TMath::Exp(par[3]*pt);
    //double fitval = par[0]*TMath::Exp(par[1]*(pt-100)) + par[2]*TMath::Exp(par[3]*(pt-100));
    //double fitval = par[0] + par[1]*(pt-10) + par[2]*(pt-10)*(pt-10) + par[3]*(pt-10)*(pt-10)*(pt-10) + par[4]*(pt-10)*(pt-10)*(pt-10)*(pt-10);
    //double fitval = par[0]
    double fitval = par[0]*TMath::Exp(-0.5*TMath::Power(((x[0] - par[1])/par[2]),2));
    return fitval;
}




void RecoRefCompare_selfCalc(){

    gStyle->SetOptFit(0111);

    bool ispp = 0;
    bool isPeriph = 0;
    bool isNom = 0;
    
    bool isC1 = 0;
    bool isC2 = 0;
    bool isC3 = 0;
    bool isC4 = 1;


    bool isCent = 0;

    TFile *f1, *f2;



    if(ispp) f1 = TFile::Open("~/Analysis/code/skimming/PYTHIA_scan/rootFiles/PYTHIA_response_noHLT_pthat15_30June22.root");
    else f1 = TFile::Open("/home/clayton/Analysis/code/skimming/PYTHIAHYDJET_scan/rootFiles/response/PYTHIAHYDJET_response_10Aug22.root");

    TH2D *h1, *h2, *h3, *h4, *h5, *h6;

    

    if(ispp){
        f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets",h1);
        f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets",h2);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets",h3);
	f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets",h4);
        f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets",h5);
        f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets",h6);
    }
    else{
        if(isC1){
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C1",h1);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C1",h2);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C1",h3);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C1",h4);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C1",h5);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C1",h6);
        }
        else if(isC2){
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C2",h1);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C2",h2);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C2",h3);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C2",h4);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C2",h5);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C2",h6);
        }
        else if(isC3){
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C3",h1);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C3",h2);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C3",h3);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C3",h4);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C3",h5);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C3",h6);
        }
        else if(isC4){
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_allJets_C4",h1);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_udJets_C4",h2);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_sJets_C4",h3);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_gJets_C4",h4);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_cJets_C4",h5);
            f1->GetObject("h_matchedRecoJetPtOverGenJetPt_genJetPt_bJets_C4",h6);
        }
	else ;
    }


   
    
    
    

    TH1D *hx_test=h2->ProjectionY();
    TAxis *xaxis = hx_test->GetXaxis();
    
    
     const int Nbins = 16;
     double pt_axis[Nbins] = {50,60,70,80,90,100,110,120,140,160,180,210,250,300,380,500};
    //const int Nbins = 15;
    //double pt_axis[Nbins] = {100,110,120,140,160,180,200,220,250,280,320,360,400,450,500};
    double x[Nbins-1];
    double ept_axis[Nbins-1];

    
    int firstxbin = 0;
    int lastxbin = 0;

    double M1[Nbins-1], M2[Nbins-1], M3[Nbins-1], M4[Nbins-1], M5[Nbins-1], M6[Nbins-1];
    double eM1[Nbins-1], eM2[Nbins-1], eM3[Nbins-1], eM4[Nbins-1], eM5[Nbins-1], eM6[Nbins-1];
    double S1[Nbins-1], S2[Nbins-1], S3[Nbins-1], S4[Nbins-1], S5[Nbins-1], S6[Nbins-1];
    double eS1[Nbins-1], eS2[Nbins-1], eS3[Nbins-1], eS4[Nbins-1], eS5[Nbins-1], eS6[Nbins-1];

    
    TH1D *h1_y, *h2_y, *h3_y, *h4_y, *h5_y, *h6_y;
    TF1 *fit1, *fit2, *fit3, *fit4, *fit5, *fit6;
    double mu1, emu1, sig1, esig1;
    double mu2, emu2, sig2, esig2;
    double mu3, emu3, sig3, esig3;
    double mu4, emu4, sig4, esig4;
    double mu5, emu5, sig5, esig5;
    double mu6, emu6, sig6, esig6;

    
    TF1 *fxn1 = new TF1("fxn1",fitfxn,0,500,3);
    TF1 *fxn2 = new TF1("fxn2",fitfxn,0,500,3);
    TF1 *fxn3 = new TF1("fxn3",fitfxn,0,500,3);
    TF1 *fxn4 = new TF1("fxn4",fitfxn,0,500,3);
    TF1 *fxn5 = new TF1("fxn5",fitfxn,0,500,3);
    TF1 *fxn6 = new TF1("fxn6",fitfxn,0,500,3);

    TCanvas *ctmp = new TCanvas("ctmp","ctmp",600,600);
    TPad *ptmp = new TPad("ptmp","ptmp",0,0,1,1);
    ptmp->SetLeftMargin(0.2);
    ptmp->SetBottomMargin(0.2);
    ctmp->cd();
    ptmp->Draw();
    ptmp->cd();
    
    
    for(int i = 0; i < Nbins-1; i++){

        x[i] = (pt_axis[i+1]+pt_axis[i])/2;
        ept_axis[i] = (pt_axis[i+1] - pt_axis[i])/2;

        cout << "i = " << i << ", x[i] = " << x[i] << endl;

        firstxbin = xaxis->FindBin(pt_axis[i]+0.001);
        lastxbin = xaxis->FindBin(pt_axis[i+1]-0.001);

        h1_y = h1->ProjectionX("h1_y",firstxbin,lastxbin); // allJets
        h2_y = h2->ProjectionX("h2_y",firstxbin,lastxbin); // udJets
        h3_y = h3->ProjectionX("h3_y",firstxbin,lastxbin); // sJets
        h4_y = h4->ProjectionX("h4_y",firstxbin,lastxbin); // gJets
        h5_y = h5->ProjectionX("h5_y",firstxbin,lastxbin); // cJets
        h6_y = h6->ProjectionX("h6_y",firstxbin,lastxbin); // bJets

        

       
        //fxn1->SetParLimits(0,0.0,1e9);
        //fxn1->SetParLimits(1,pt_axis[i]-50,pt_axis[i+1]+50.0);
        //fxn1->SetParLimits(2,0.0,TMath::Sqrt(pt_axis[i+1]-pt_axis[i]));
        
       
        // fxn1->SetParameter(0,1.0);
        // fxn1->SetParameter(1,(pt_axis[i+1]+pt_axis[i])/2);
        // fxn1->SetParameter(2,(pt_axis[i+1]-pt_axis[i]));

        double sigmaParMax = 0.6;
   
     
        double maxBinVal_1 = h1_y->GetMean();
        double maxBinVal_2 = h2_y->GetMean();
        double maxBinVal_3 = h3_y->GetMean();
        double maxBinVal_4 = h4_y->GetMean();
        double maxBinVal_5 = h5_y->GetMean();
        double maxBinVal_6 = h6_y->GetMean();

        double initStdDev_1 = h1_y->GetStdDev();
        double initStdDev_2 = h2_y->GetStdDev();
        double initStdDev_3 = h3_y->GetStdDev();
        double initStdDev_4 = h4_y->GetStdDev();
        double initStdDev_5 = h5_y->GetStdDev();
        double initStdDev_6 = h6_y->GetStdDev();

        
        fxn1->SetParameter(0,1.0);
        fxn1->SetParameter(1,maxBinVal_1);
        fxn1->SetParameter(2,initStdDev_1);
        fxn1->SetParLimits(1,0.8,1.2);
        fxn1->SetParLimits(2,0.0,sigmaParMax);
        h1_y->Fit(fxn1,"R","",0,5);
        mu1 = fxn1->GetParameter(1);
        emu1 = fxn1->GetParError(1);
        sig1 = fxn1->GetParameter(2);
        esig1 = fxn1->GetParError(2);

        fxn2->SetParameter(0,1.0);
        fxn2->SetParameter(1,maxBinVal_2);
        fxn2->SetParameter(2,initStdDev_2);
        fxn2->SetParLimits(1,0.8,1.2);
        fxn2->SetParLimits(2,0.0,sigmaParMax);
        h2_y->Fit(fxn2,"R","",0,5);
        mu2 = fxn2->GetParameter(1);
        emu2 = fxn2->GetParError(1);
        sig2 = fxn2->GetParameter(2);
        esig2 = fxn2->GetParError(2);

        fxn3->SetParameter(0,1.0);
        fxn3->SetParameter(1,maxBinVal_3);
        fxn3->SetParameter(2,initStdDev_3);
        fxn3->SetParLimits(1,0.8,1.2);
        fxn3->SetParLimits(2,0.0,sigmaParMax);
        h3_y->Fit(fxn3,"R","",0,5);
        mu3 = fxn3->GetParameter(1);
        emu3 = fxn3->GetParError(1);
        sig3 = fxn3->GetParameter(2);
        esig3 = fxn3->GetParError(2);

        fxn4->SetParameter(0,1.0);
        fxn4->SetParameter(1,maxBinVal_4);
        fxn4->SetParameter(2,initStdDev_4);
        fxn4->SetParLimits(1,0.8,1.2);
        fxn4->SetParLimits(2,0.0,sigmaParMax);
        h4_y->Fit(fxn4,"R","",0,5);
        mu4 = fxn4->GetParameter(1);
        emu4 = fxn4->GetParError(1);
        sig4 = fxn4->GetParameter(2);
        esig4 = fxn4->GetParError(2);

        fxn5->SetParameter(0,1.0);
        fxn5->SetParameter(1,maxBinVal_5);
        fxn5->SetParameter(2,initStdDev_5);
        fxn5->SetParLimits(1,0.8,1.2);
        fxn5->SetParLimits(2,0.0,sigmaParMax);
        h5_y->Fit(fxn5,"R","",0,5);
        mu5 = fxn5->GetParameter(1);
        emu5 = fxn5->GetParError(1);
        sig5 = fxn5->GetParameter(2);
        esig5 = fxn5->GetParError(2);

        fxn6->SetParameter(0,1.0);
        fxn6->SetParameter(1,maxBinVal_6);
        fxn6->SetParameter(2,initStdDev_6);
        fxn6->SetParLimits(1,0.8,1.2);
        fxn6->SetParLimits(2,0.0,sigmaParMax);
        h6_y->Fit(fxn6,"R","",0,5);
        mu6 = fxn6->GetParameter(1);
        emu6 = fxn6->GetParError(1);
        sig6 = fxn6->GetParameter(2);
        esig6 = fxn6->GetParError(2);

     

        h4_y->GetXaxis()->SetTitle("p_{T}^{recoJet} / p_{T}^{genJet}");
        h4_y->GetYaxis()->SetTitle("Weighted entries");
        h4_y->Draw();
        fxn4->SetLineColor(kRed);
        fxn4->Draw("same");
        h4_y->SetTitle(Form("%3.0f < p_{T}^{genJet} < %3.0f, bJets",pt_axis[i],pt_axis[i+1]));
        ctmp->SaveAs(Form("tmp/bJets/tmp_%3.0f_%3.0f.pdf",pt_axis[i],pt_axis[i+1]));


               
        M1[i] = mu1;
        eM1[i] = emu1; 
        M2[i] = mu2;
        eM2[i] = emu2; 
        M3[i] = mu3;
        eM3[i] = emu3; 
        M4[i] = mu4;
        eM4[i] = emu4; 
        M5[i] = mu5;
        eM5[i] = emu5;
        M6[i] = mu6;
        eM6[i] = emu6;

        if(i==2) cout << mu3;
        

        S1[i] = sig1;
        eS1[i] = esig1;
        S2[i] = sig2;
        eS2[i] = esig2;
        S3[i] = sig3;
        eS3[i] = esig3;
        S4[i] = sig4;
        eS4[i] = esig4;
        S5[i] = sig5;
        eS4[i] = esig5;
        S6[i] = sig6;
        eS6[i] = esig6;
        



    }

    

    double markerSize = 1.1;

    
    
    TGraphErrors *gr2 = new TGraphErrors(Nbins-1,x,M2,ept_axis,eM2);
    
    gr2->SetMarkerStyle(24);
    gr2->SetMarkerColor(kBlack);
    gr2->SetLineColor(kBlack);
    gr2->SetMarkerSize(markerSize);

    TGraphErrors *gr3 = new TGraphErrors(Nbins-1,x,M3,ept_axis,eM3);
    
    gr3->SetMarkerStyle(26);
    gr3->SetMarkerColor(kRed);
    gr3->SetLineColor(kRed);
    gr3->SetMarkerSize(markerSize);

    TGraphErrors *gr4 = new TGraphErrors(Nbins-1,x,M4,ept_axis,eM4);
    
    gr4->SetMarkerStyle(25);
    gr4->SetMarkerColor(kBlue);
    gr4->SetLineColor(kBlue);
    gr4->SetMarkerSize(markerSize);

    TGraphErrors *gr5 = new TGraphErrors(Nbins-1,x,M5,ept_axis,eM5);
    
    gr5->SetMarkerStyle(28);
    gr5->SetMarkerColor(kGreen+2);
    gr5->SetLineColor(kGreen+2);
    gr5->SetMarkerSize(markerSize);

    TGraphErrors *gr6 = new TGraphErrors(Nbins-1,x,M6,ept_axis,eM6);
    
    gr6->SetMarkerStyle(32);
    gr6->SetMarkerColor(kMagenta);
    gr6->SetLineColor(kMagenta);
    gr6->SetMarkerSize(markerSize);

    TGraphErrors *Gr1 = new TGraphErrors(Nbins-1,x,S1,ept_axis,eS1);
    
    Gr1->SetMarkerStyle(24);
    Gr1->SetMarkerColor(kBlue+3);
    Gr1->SetLineColor(kBlue+3);
    Gr1->SetMarkerSize(markerSize);


    TGraphErrors *Gr2 = new TGraphErrors(Nbins-1,x,S2,ept_axis,eS2);
    
    Gr2->SetMarkerStyle(24);
    Gr2->SetMarkerColor(kBlack);
    Gr2->SetLineColor(kBlack);
    Gr2->SetMarkerSize(markerSize);

    TGraphErrors *Gr3 = new TGraphErrors(Nbins-1,x,S3,ept_axis,eS3);
    
    Gr3->SetMarkerStyle(26);
    Gr3->SetMarkerColor(kRed);
    Gr3->SetLineColor(kRed);
    Gr3->SetMarkerSize(markerSize);

    TGraphErrors *Gr4 = new TGraphErrors(Nbins-1,x,S4,ept_axis,eS4);
    
    Gr4->SetMarkerStyle(25);
    Gr4->SetMarkerColor(kBlue);
    Gr4->SetLineColor(kBlue);
    Gr4->SetMarkerSize(markerSize);

    TGraphErrors *Gr5 = new TGraphErrors(Nbins-1,x,S5,ept_axis,eS5);
    
    Gr5->SetMarkerStyle(28);
    Gr5->SetMarkerColor(kGreen+2);
    Gr5->SetLineColor(kGreen+2);
    Gr5->SetMarkerSize(markerSize);

    TGraphErrors *Gr6 = new TGraphErrors(Nbins-1,x,S6,ept_axis,eS6);
    
    Gr6->SetMarkerStyle(32);
    Gr6->SetMarkerColor(kMagenta);
    Gr6->SetLineColor(kMagenta);
    Gr6->SetMarkerSize(markerSize);

    // Gr6->SetMarkerStyle(24);
    // Gr6->SetMarkerColor(kGreen);
    // Gr6->SetLineColor(kGreen);
    // Gr6->SetMarkerSize(markerSize);




   TCanvas *c1 = new TCanvas("c1","c1",500,500);
   c1->cd();
   TPad *p1 = new TPad("p1","p1",0,0,1,1);
   p1->SetLeftMargin(0.2);
   p1->SetBottomMargin(0.15);
   p1->Draw();
   p1->cd();
   TMultiGraph  *mg  = new TMultiGraph();
  
   mg->Add(gr2);
   mg->Add(gr3);
   mg->Add(gr4);
   mg->Add(gr5);
   mg->Add(gr6);

   double labelSize = 0.055;
   double titleSize = 0.065;
   int NdivVal = 508;

   mg->GetYaxis()->SetRangeUser(0.8,1.2);
   mg->GetXaxis()->SetRangeUser(0,500);
   mg->GetYaxis()->SetTitle("#mu (#font[52]{p}_{T}^{recoJet} / #font[52]{p}_{T}^{genJet})");
   mg->GetXaxis()->SetTitle("#font[52]{p}_{T}^{genJet} [GeV]");
   mg->GetXaxis()->SetLabelSize(labelSize);
   mg->GetXaxis()->SetTitleSize(titleSize);
   mg->GetYaxis()->SetLabelSize(labelSize);
   mg->GetYaxis()->SetTitleSize(titleSize);
   mg->GetXaxis()->SetTitleFont(42);
   mg->GetXaxis()->SetLabelFont(42);
   mg->GetYaxis()->SetTitleFont(42);
   mg->GetYaxis()->SetLabelFont(42);
   mg->GetYaxis()->SetNdivisions(NdivVal);
   mg->GetXaxis()->SetNdivisions(NdivVal);

   TLegend *l = new TLegend(0.7,0.7,0.89,0.88);
   l->SetBorderSize(0);
   l->SetTextFont(42);
   l->SetTextSize(0.04);
   l->AddEntry(gr2,"#font[52]{ud} jets","p");
   l->AddEntry(gr3,"#font[52]{s} jets","p");
   l->AddEntry(gr4,"#font[52]{g} jets","p");
   l->AddEntry(gr5,"#font[52]{c} jets","p");
   //l->AddEntry(gr1,"all jets","p");
   l->AddEntry(gr6,"#font[52]{b} jets","p");

    double x_t1 = 0.23;
    double y_t1 = 0.92;
    double y_t2 = 0.85;
    double y_t3 = 0.79;
    double y_t4 = 0.73;

   mg->Draw("AP");
   l->Draw();


    TLine *line = new TLine(pt_axis[0],1,pt_axis[Nbins-1],1);
    line->SetLineStyle(7);
    line->Draw();





     TLatex *t1 = new TLatex();
     t1->SetTextFont(42);
     t1->SetTextSize(0.037);
     //t1->DrawLatexNDC(x_t1-0.02,y_t1,"Jet energy scale");
     t1->SetTextSize(0.044);
     
     if(!ispp){
     if(isC1) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 0-10%");
     else if(isC2) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 10-30%");
     else if(isC3) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 30-50%");
     else if(isC4) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 50-90%");
     }
     else t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA");

     t1->DrawLatexNDC(x_t1,y_t3,"|#eta^{jet}| < 1.6");
     //t1->DrawLatexNDC(x_t1,y_t4,"#font[52]{p}_{T}^{#mu} > 7 GeV, |#eta^{#mu}| < 2.0");
   
    
    
    
    
    
   TCanvas *c2 = new TCanvas("c2","c2",500,500);
   c2->cd();
   TPad *p2 = new TPad("p2","p2",0,0,1,1);
   p2->SetLeftMargin(0.2);
   p2->SetBottomMargin(0.15);
   p2->Draw();
   p2->cd();
   TMultiGraph  *mg2  = new TMultiGraph();
  
   mg2->Add(Gr2);
   mg2->Add(Gr3);
   mg2->Add(Gr4);
   mg2->Add(Gr5);
   mg2->Add(Gr6);
   //mg2->Add(Gr1);

   //mg2->GetYaxis()->SetRangeUser(0,0.54);
   mg2->GetYaxis()->SetRangeUser(0,0.5);
   mg2->GetXaxis()->SetRangeUser(0,500);
   mg2->GetYaxis()->SetTitle("#sigma (#font[52]{p}_{T}^{recoJet} / #font[52]{p}_{T}^{genJet})");
   mg2->GetXaxis()->SetTitle("#font[52]{p}_{T}^{genJet} [GeV]");
   mg2->GetXaxis()->SetLabelSize(labelSize);
   mg2->GetXaxis()->SetTitleSize(titleSize);
   mg2->GetYaxis()->SetLabelSize(labelSize);
   mg2->GetYaxis()->SetTitleSize(titleSize);
   mg2->GetXaxis()->SetTitleFont(42);
   mg2->GetXaxis()->SetLabelFont(42);
   mg2->GetYaxis()->SetTitleFont(42);
   mg2->GetYaxis()->SetLabelFont(42);
   mg2->GetYaxis()->SetNdivisions(NdivVal);
   mg2->GetXaxis()->SetNdivisions(NdivVal);
    
    mg2->Draw("AP");
   l->Draw();
    
   
 

    // t1->DrawLatexNDC(x_t1-0.02,y_t1,"Jet energy resolution");
   
     t1->SetTextSize(0.044);
 
    if(!ispp){
     if(isC1) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 0-10%");
     else if(isC2) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 10-30%");
     else if(isC3) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 30-50%");
     else if(isC4) t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA+HYDJET 50-90%");
     }
     else t1->DrawLatexNDC(x_t1,y_t2,"PYTHIA");
	t1->DrawLatexNDC(x_t1,y_t3,"|#eta^{jet}| < 1.6");
	//t1->DrawLatexNDC(x_t1,y_t4,"#font[52]{p}_{T}^{#mu} > 7 GeV, |#eta^{#mu}| < 2.0");



    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /*
    
    
    TH1D *h2_1_r = (TH1D*) h2_1->Rebin(Nbins-1,"h2_1_r",pt_axis);
    TH1D *h3_1_r = (TH1D*) h3_1->Rebin(Nbins-1,"h3_1_r",pt_axis);
    TH1D *h4_1_r = (TH1D*) h4_1->Rebin(Nbins-1,"h4_1_r",pt_axis);

    TH1D *h2_2_r = (TH1D*) h2_2->Rebin(Nbins-1,"h2_2_r",pt_axis);
    TH1D *h3_2_r = (TH1D*) h3_2->Rebin(Nbins-1,"h3_2_r",pt_axis);
    TH1D *h4_2_r = (TH1D*) h4_2->Rebin(Nbins-1,"h4_2_r",pt_axis);

    // normalize by the new bin widths
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    for(int i = 0; i < h2_1_r->GetSize(); i++){
        x = h2_1_r->GetBinWidth(i);
        y = rebinXval*h2_1_r->GetBinContent(i);
        z = rebinXval*h2_1_r->GetBinError(i);
        if(x != 0.0){
            h2_1_r->SetBinContent(i,y/x);
            h2_1_r->SetBinError(i,z/x);
        }
    }

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for(int i = 0; i < h3_1_r->GetSize(); i++){
        x = h3_1_r->GetBinWidth(i);
        y = rebinXval*h3_1_r->GetBinContent(i);
        z = rebinXval*h3_1_r->GetBinError(i);
        if(x != 0.0){
            h3_1_r->SetBinContent(i,y/x);
            h3_1_r->SetBinError(i,z/x);
        }
    }

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for(int i = 0; i < h4_1_r->GetSize(); i++){
        x = h4_1_r->GetBinWidth(i);
        y = rebinXval*h4_1_r->GetBinContent(i);
        z = rebinXval*h4_1_r->GetBinError(i);
        if(x != 0.0){
            h4_1_r->SetBinContent(i,y/x);
            h4_1_r->SetBinError(i,z/x);
        }
    }

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for(int i = 0; i < h2_2_r->GetSize(); i++){
        x = h2_2_r->GetBinWidth(i);
        y = rebinXval*h2_2_r->GetBinContent(i);
        z = rebinXval*h2_2_r->GetBinError(i);
        if(x != 0.0){
            h2_2_r->SetBinContent(i,y/x);
            h2_2_r->SetBinError(i,z/x);
        }
    }

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for(int i = 0; i < h3_2_r->GetSize(); i++){
        x = h3_2_r->GetBinWidth(i);
        y = rebinXval*h3_2_r->GetBinContent(i);
        z = rebinXval*h3_2_r->GetBinError(i);
        if(x != 0.0){
            h3_2_r->SetBinContent(i,y/x);
            h3_2_r->SetBinError(i,z/x);
        }
    }

    x = 0.0;
    y = 0.0;
    z = 0.0;

    for(int i = 0; i < h4_2_r->GetSize(); i++){
        x = h4_2_r->GetBinWidth(i);
        y = rebinXval*h4_2_r->GetBinContent(i);
        z = rebinXval*h4_2_r->GetBinError(i);
        if(x != 0.0){
            h4_2_r->SetBinContent(i,y/x);
            h4_2_r->SetBinError(i,z/x);
        }
    }
    

    // some style options
    h1->SetTitle("");
    h1->GetXaxis()->SetTitle("Gen jet p_{T} [GeV/c]");
    h1->GetYaxis()->SetTitle("Reco jet p_{T} [GeV/c]");
    h2->SetTitle("");
    h2->GetXaxis()->SetTitle("Gen jet p_{T} [GeV/c]");
    h2->GetYaxis()->SetTitle("(Reco/Gen jet p_{T} [GeV/c]");
    // h3->SetTitle("");
    // h3->GetXaxis()->SetTitle("Ref jet p_{T} [GeV/c]");
    // h3->GetYaxis()->SetTitle("(Reco - Ref jet p_{T}) / Ref jet p_{T}");

    TCanvas *c1 = new TCanvas("c1","Reco jet pt vs Ref jet pt",600,600);
    TCanvas *c2 = new TCanvas("c2","|Reco - Ref| vs Ref jet pt",600,600);
   // TCanvas *c3 = new TCanvas("c3","|Reco - Ref|/Ref vs Ref jet pt",700,500);
    TCanvas *c4 = new TCanvas("c4","Mean from fitSlicesY",600,600);
    TCanvas *c5 = new TCanvas("c5","Std Dev from fitSlicesY",600,600);

    c1->SetLogz(1);
    c2->SetLogz(1);
   // c3->SetLogz(1);

    c1->cd();
    h1->GetYaxis()->SetRangeUser(60,500);
    h1->GetXaxis()->SetRangeUser(0,500);
    //h1->GetZaxis()->SetRangeUser(1e-7,0.7e1);
    h1->SetStats(0);
    h1->Draw("colz");

    c2->cd();
    h2->GetYaxis()->SetRangeUser(0,5);
    h2->GetXaxis()->SetRangeUser(0,500);
    //h2->GetZaxis()->SetRangeUser(1e-7,0.7e1);
    h2->SetStats(0);
    h2->Draw("colz");

    // c3->cd();
    // h3->GetXaxis()->SetRangeUser(50,500);
    // h3->SetStats(0);
    // h3->Draw("colz");

    
    
    

    h2_1_r->SetLineColor(kBlack);
    h2_1_r->SetMarkerColor(kBlack);
    h2_1_r->SetMarkerStyle(24);
    h2_1_r->SetMarkerSize(0.8);
    h2_1_r->SetStats(0);

    h2_2_r->SetLineColor(kBlack);
    h2_2_r->SetMarkerColor(kBlack);
    h2_2_r->SetMarkerStyle(24);
    h2_2_r->SetMarkerSize(0.8);
    h2_2_r->SetStats(0);

    h3_1_r->SetLineColor(kBlue);
    h3_1_r->SetMarkerColor(kBlue);
    h3_1_r->SetMarkerStyle(25);
    h3_1_r->SetMarkerSize(0.8);

    h3_2_r->SetLineColor(kBlue);
    h3_2_r->SetMarkerColor(kBlue);
    h3_2_r->SetMarkerStyle(25);
    h3_2_r->SetMarkerSize(0.8);

    h4_1_r->SetLineColor(kRed);
    h4_1_r->SetMarkerColor(kRed);
    h4_1_r->SetMarkerStyle(26);
    h4_1_r->SetMarkerSize(0.8);

    h4_2_r->SetLineColor(kRed);
    h4_2_r->SetMarkerColor(kRed);
    h4_2_r->SetMarkerStyle(26);
    h4_2_r->SetMarkerSize(0.8);

    c4->cd();
    TPad *pad4 = new TPad("pad4","pad4",0,0,1,1);
    double lm = 0.2;
    double bm = 0.18;
    double tm = 0.0;
    pad4->SetLeftMargin(lm);
    pad4->SetBottomMargin(bm);
    pad4->Draw();            
    pad4->cd();
    double xt = 0.055;
    double xl = 0.04;
    double yt = 0.055;
    double yl = 0.04;
    h2_1_r->GetXaxis()->SetTitleSize(xt);
    h2_1_r->GetXaxis()->SetLabelSize(xl);
    h2_1_r->GetYaxis()->SetTitleSize(yt);
    h2_1_r->GetYaxis()->SetLabelSize(yl);
    h2_1_r->SetTitle("");
    h2_1_r->GetXaxis()->SetTitle("Gen jet p_{T} [GeV]");
    h2_1_r->GetYaxis()->SetTitle("#mu (reco p_{T} / gen p_{T})");
    h2_1_r->GetXaxis()->SetRangeUser(60,500);
    h2_1_r->GetYaxis()->SetRangeUser(0.8,1.2);
    h2_1_r->Draw();
    h3_1_r->Draw("same");
    h4_1_r->Draw("same");
    TLegend *leg = new TLegend(0.6,0.73,0.85,0.88);
    leg->SetBorderSize(0);
    leg->AddEntry(h2_1_r,"Inclusive");
    leg->AddEntry(h3_1_r,"Muon-tagged");
    leg->AddEntry(h4_1_r,"Untagged");
    leg->Draw();

    TLatex *l = new TLatex();
    double textSize = 0.03;
    double textX = 0.24;
    double textY1 = 0.85;
    double textY2 = 0.80;
    double textY3 = 0.75;
    l->SetTextFont(42);
    l->SetTextSize(textSize);
    if(isPeriph && isC0) l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 30-50%");
    else if(isPeriph && isC1)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 50-90%");
    else if(isNom && isC1)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 10-30%");
    else if(isNom && isC0)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 0-10%");
    l->DrawLatexNDC(textX,textY2,"|#eta^{jet}| < 1.6");
   // l->DrawLatexNDC(textX,textY3,"L3Mu5 HLT");
    //h2_1->Fit("pol9");
    //TF1 *func = new TF1("func",JES_fit,30.,500.,6);
    //h2_1->Fit(func,"M R WL","",100.,500.);


    c5->cd();
    TPad *pad5 = new TPad("pad5","pad5",0,0,1,1);
    pad5->SetLeftMargin(lm);
    pad5->SetBottomMargin(bm);
    pad5->Draw();           
    pad5->cd();
    h2_2_r->GetXaxis()->SetTitleSize(xt);
    h2_2_r->GetXaxis()->SetLabelSize(xl);
    h2_2_r->GetYaxis()->SetTitleSize(yt);
    h2_2_r->GetYaxis()->SetLabelSize(yl);
    h2_2_r->SetTitle("");
    h2_2_r->GetXaxis()->SetTitle("Gen jet p_{T} [GeV]");
    h2_2_r->GetYaxis()->SetTitle("#sigma (reco p_{T} / gen p_{T})");
    h2_2_r->GetXaxis()->SetRangeUser(60,500);
    h2_2_r->GetYaxis()->SetRangeUser(0.0,0.35);
    h2_2_r->Draw();
    h3_2_r->Draw("same");
    h4_2_r->Draw("same");
    //h2_2->Fit("pol9");
    TF1 *func = new TF1("func",JES_fit,30.,500.,7);
   // h2_2->Fit(func,"M R WL","",70.,500.);
   leg->Draw();
   if(isPeriph && isC0) l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 30-50%");
    else if(isPeriph && isC1)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 50-90%");
    else if(isNom && isC1)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 10-30%");
    else if(isNom && isC0)  l->DrawLatexNDC(textX,textY1,"PYTHIA+HYDJET 0-10%");
    l->DrawLatexNDC(textX,textY2,"|#eta^{jet}| < 1.6");
   // l->DrawLatexNDC(textX,textY3,"L3Mu5 HLT");




    // TString output;
    // if(ispp) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA/JER.root";
    // else{
    //     if(isC0) {
    //         if(isCent) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/0to5/JER.root";
    //         else if(isNom) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/0to10/JER.root";
    //         else if(isPeriph) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/30to50/JER.root";
    //         else;
    //     }

    //     else if(isC1){
    //         if(isCent) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/5to10/JER.root";
    //         else if(isNom) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/10to30/JER.root";
    //         else if(isPeriph) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/50to70/JER.root";
    //         else ;
    //     } 

    //     else if(isC2) {
    //         if(isCent) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/10to30/JER.root";
    //         else if(isNom) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/30to90/JER.root";
    //         else if(isPeriph) output = "~/Analysis/code/RecoRefCompare/JER/PYTHIA+HYDJET/70to90/JER.root";
    //     }

    //     else ;
    // }
    



*/



}
