void LogEnergyMuon(){

  cout<<"Function test...  "<<endl;
  
  Double_t betaf=0.9;
  
  TF1 *fPHI    = new TF1("fPHI",  PHI,   15, 22, 0);
  
  TF1 *fIT     = new TF1("fIT",   EFIT,  15, 22, 2);
  TF1 *fE      = new TF1("fE",    FE,    15, 22, 1);
  TF1 *fEFI    = new TF1("fEFI",  FEFI,  15, 22, 1);
  TF1 *fE2     = new TF1("fE2",   FE2,   15, 22, 1);
  
  TF1 *fITMU   = new TF1("fITMU", MUFIT,  0, 11, 2);
  //TF1 *fITMU   = new TF1("fITMU", FMU,  0, 11, 1);
  TF1 *fMU     = new TF1("fMU",   FMU,    0, 11, 1);
  TF1 *fMUFI   = new TF1("fMUFI", FMUFI,  0, 11, 1);
  TF1 *fMU2    = new TF1("fMU2",  FMU2,   0, 11, 1);
  
  fIT->SetNpx(22000);
  fITMU->SetNpx(22000);
  
  fE->SetNpx(22000);
  fE2->SetNpx(22000);  
  fEFI->SetNpx(22000);
  
  fMU->SetNpx(22000);
  fMU2->SetNpx(22000);  
  fMUFI->SetNpx(22000);
  
  fMU->SetParameter(0,7);
  fMUFI->SetParameter(0,7);
  
  TCanvas *ctest =new TCanvas ("ctest","ctest",1200,600);
  ctest->Divide(2,1);
  ctest->cd(1);
  fMU->Draw();
  ctest->cd(2);
  fITMU->SetParameter(1,3e14);
  fMUFI->Draw();
  
  ctest->Update();
  
  cout<<fMU->Eval(7)/fMUFI->Eval(7)<<endl;
  
  fE->SetParameter(0,18.5);
  fEFI->SetParameter(0,18.5);
 
  TCanvas *ctest1 =new TCanvas ("ctest1","ctest1",1200,600);
  ctest1->Divide(2,1);
  ctest1->cd(1);
  fE->Draw();
  ctest1->cd(2);
  fEFI->Draw();
  
  ctest1->Update();
  
  cout<<fE->Eval(18.5)/fEFI->Eval(18.5)<<endl;
  
  
  TF1 *fSIGMAMUONS = new TF1("fSIGMAMUONS", SIGMAMUONS, 15, 22, 0);
  TF1 *fMUONS = new TF1("fMUONS", MUONS, 15, 22, 1);
  fMUONS->FixParameter(0,betaf);
  
  Double_t BL=17.85, BU=20.25;
  Double_t BS=0.1;
  Int_t binx = TMath::Ceil((BU-BL)/BS);
  
  
  TH1D *hPHI = new TH1D("hPHI", "hPHI", binx, BL, BU);
  TH2D *hE   = new TH2D("hE", "hE",     binx, BL, BU, 1000000, TMath::Power(10,15.),TMath::Power(10,21.));
  TH2D *hE2  = new TH2D("hE2", "hE2",   binx, BL, BU, 1000, 15, 25);
  
  TH2D *hMU  = new TH2D("hMU","hMU",    binx, BL, BU, 1000000, TMath::Power(10,4.),TMath::Power(10,11.));  
  TH2D *hMU2  = new TH2D("hMU2","hMU2", binx, BL, BU, 1000, 1,11);
  
  Double_t logEneT=0, logEneR=0;
  Double_t EneT=0, EneR=0;
  Double_t flux=0;
  Double_t MuonTrue=0, MuonRec=0, logMuonRec=0, SigmaMuonT=0;
  
  
  Int_t k=0;
  while (k<10000000){
    
    //lognormal fluctuation
    logEneT = fPHI->GetRandom(17.,21.);
    logEneR = logEneT + gRandom->Gaus(0,log10(1.2));
    EneT = TMath::Power(10,logEneT);
    EneR = TMath::Power(10,logEneR);
    //EneR = EneT + gRandom->Gaus(0,EneT*0.2);
    //logEneR = TMath::Log10(EneR);
    
    MuonTrue   = fMUONS->Eval(logEneT);
    SigmaMuonT = fSIGMAMUONS->Eval(logEneT);
    //MuonTrue   = MuonTrue + gRandom->Gaus(0,SigmaMuonT);
    //MuonTrue   = MuonTrue + gRandom->Gaus(0,MuonTrue*0.2);
    MuonRec    = MuonTrue;
    
    logMuonRec = TMath::Log10(MuonRec);
    
    
    //cout<<MuonRec<<"   "<<logMuonRec<<endl;
    //getchar();
    flux = fPHI->Eval(logEneT);
    hPHI->Fill(logEneT);
    
    Double_t e = Round(logEneR);
    
    hE->Fill(e,EneT);
    hE2->Fill(e,logEneT);
    
    hMU->Fill(e,MuonRec);
    hMU2->Fill(e,logMuonRec);
    
    
    k++;
  }
  
  TCanvas *c21 = new TCanvas ("c21", "c21", 800,500);
  TCanvas *c22 = new TCanvas ("c22", "c22", 800,500);
  
  
  const int dim = binx;
  Int_t yent = 0, myent = 0, sument = 0;
  Double_t MeanH[dim], RMSH[dim], EnH[dim], MuH[dim];
  Double_t MeanMuH[dim], RMSMuH[dim];
  Double_t B = TMath::Power(10,-10);
  Double_t beta = 0.90;

  
  for(int i=1;i<binx+1;i++){
    
    EnH[i-1] = hE->GetBinCenter(i);
    MuH[i-1] = TMath::Log10(B * TMath::Power(TMath::Power(10,EnH[i-1]),beta));
    
    cout<<"MuH for fit:  "<<MuH[i-1]<<endl;
     
    fIT->FixParameter(0,   EnH[i-1]);
    fITMU->FixParameter(0, MuH[i-1]);
    fIT->SetParameter(1, 1e32);
    fITMU->SetParameter(1, 1e40);
    
    TH1D *py  = hE->ProjectionY("py",   i, i);
    TH1D *py1 = hE2->ProjectionY("py1", i, i);
    py1->Fit(fIT,"R");
    
    TH1D *my  = hMU->ProjectionY("my",   i, i);
    TH1D *my1 = hMU2->ProjectionY("my1", i, i);
    my1->Fit(fITMU,"R");

    myent = my->GetEntries();
    yent  = py->GetEntries();
    
    c21->cd();
    py1->Draw();
    c21->Update();
    c22->cd();
    my1->Draw();
    //fITMU->Draw("SAME");
    c22->Update();
    //getchar();
     
    MeanH[i-1] = py->GetMean();
    RMSH[i-1]  = py->GetRMS();
    
    MeanMuH[i-1] = my->GetMean();
    RMSMuH[i-1]  = my->GetRMS();
    cout<<endl<<"-----------------------------------------------------------------"<<endl;
    //cout<<yent<<"  "<<EnH[i-1]<<"    "<<log10(MeanH[i-1])<<"    "<<RMSH[i-1]<<"    "<<log10(MeanMuH[i-1])<<"    "<<RMSMuH[i-1]<<endl;
    cout<<"--------------------------------------------------------"<<endl<<endl;
    
    sument = sument + yent;
    delete py;
    delete py1;
    delete my;
    delete my1;

  }
  
  cout<<sument<<endl;
  
  
  
  Double_t Ei=18.0, Ef=22.0, Es=0.025;
  Int_t dim1 = (Ef-Ei)/Es;
  const int dim = dim1;
  Double_t MeanE[dim],  MeanE2[dim], RMSE[dim], VRMSE[dim];
  Double_t MeanMU[dim], MeanMU2[dim], RMSMU[dim], VRMSMU[dim];
  Double_t Ene[dim];
  Int_t ke = 0;
  Double_t mu=0;
  
  for(Double_t e=Ei; e<Ef; e+=Es){   
    
    Ene[ke] = e;
    mu = TMath::Log10(B * TMath::Power(TMath::Power(10,e),beta));
    
    //cout<<mu<<endl;
    
    fE->SetParameter(0,e);
    fE2->SetParameter(0,e);
    fEFI->SetParameter(0,e);
    
    fMU->SetParameter(0,mu);
    fMU2->SetParameter(0,mu);
    fMUFI->SetParameter(0,mu);
    
    
    //TH1D *hTEST  = new TH1D("hTEST", "hTEST", 120, 16, 22);
    //TCanvas *cl = new TCanvas("cl","cl",700,600);
    //cl->cd();
    
    Double_t nE1=0, nE2=0, dE1=0;  
    Double_t nMU1=0, nMU2=0, dMU1=0;  
    Double_t imu=0;
    for(double ie=16; ie<22.;ie+=0.001){
      
      imu = TMath::Log10(B * TMath::Power(TMath::Power(10,ie),beta));
      //cout<<imu<<endl;
      
      nE1 += fE->Eval(ie);
      nE2 += fE2->Eval(ie);
      dE1 += fEFI->Eval(ie);
      
      nMU1 += fMU->Eval(imu);
      nMU2 += fMU2->Eval(imu);
      dMU1 += fMUFI->Eval(imu);
      
      
      //hTEST->Fill(log10(nE1/dE1));
      
    }
    // cout<<nMU1/dMU1<<endl;
    
    //cl->cd();
    //hTEST->Draw();
    //cl->Update();
    //getchar();
    
    //delete hTEST;
    //delete cl;
    
    MeanE[ke]  = nE1/dE1;
    MeanE2[ke] = nE2/dE1;
    RMSE[ke]   = MeanE2[ke] - TMath::Power(MeanE[ke],2); 
    VRMSE[ke]  = TMath::Sqrt(RMSE[ke]);
    
    cout<<"Values from Sum En Mean:  "<<log10(MeanE[ke])<<"   "<<MeanE2[ke]<<"   "<<TMath::Sqrt(RMSE[ke])<<endl;
    
    
    MeanMU[ke]  = nMU1/dMU1;
    MeanMU2[ke] = nMU2/dMU1;
    RMSMU[ke]   = MeanMU2[ke] - TMath::Power(MeanMU[ke],2); 
    VRMSMU[ke]  = TMath::Sqrt(RMSMU[ke]);
    
    cout<<"------ Values from Sum MU Mean:  "<<log10(MeanMU[ke])<<"   "<<MeanMU[ke]<<"   "<<TMath::Sqrt(RMSMU[ke])<<endl;
    
    
    ke++;
    
  }
  
  TGraph *gE = new TGraph(dim1,  Ene, MeanE);
  formatTGraph(gE, 20, 0.5, 2, "E_{Rec}", "<E_{True}>");
  
  TGraph *gRE = new TGraph(dim1, Ene, VRMSE);
  formatTGraph(gRE, 20, 0.5, 2, "E_{Rec}", "RMS(E_{True})");
  
  TGraph *GHEr1 = new TGraph(binx, EnH, MeanH);
  formatTGraph(GHEr1, 21, 0.7, 4, "E_{Rec}", "<E_{True}>");
  
  TGraph *GHEr2 = new TGraph(binx, EnH, RMSH);
  formatTGraph(GHEr2, 21, 0.7, 4, "E_{Rec}", "RMS(E_{True})");
   
  TGraph *gMU = new TGraph(dim1,  Ene, MeanMU);
  formatTGraph(gMU, 20, 0.5, 2, "E_{Rec}", "<N_{#mu}>");
  
  TGraph *gRMU = new TGraph(dim1, Ene, VRMSMU);
  formatTGraph(gRMU, 20, 0.5, 2, "E_{Rec}", "RMS(N_{#mu})");
  
  TGraph *GHMUr1 = new TGraph(binx, EnH, MeanMuH);
  formatTGraph(GHMUr1, 21, 0.7, 4,  "E_{Rec}", "<N_{#mu}>");
  
  TGraph *GHMUr2 = new TGraph(binx, EnH, RMSMuH);
  formatTGraph(GHMUr2, 21, 0.7, 4,  "E_{Rec}", "RMS(N_{#mu})");
  
  
  TCanvas *c0 = new TCanvas ("c0", "c0", 600,600);
  c0->cd()->SetLogy();
  hE->Draw();
  
  TCanvas *c1 = new TCanvas ("c1", "c1", 600,600);
  c1->cd()->SetLogy();
  hPHI->Draw();
  
  TCanvas *c3 = new TCanvas ("c3", "c3", 1400,600);
  c3->Divide(2,1);
  
  c3->cd(1)->SetLogy();
  GHEr1->Draw("AP");
  gE->Draw("SAMEP");
  
  c3->cd(2)->SetLogy();
  GHEr2->Draw("AP");
  gRE->Draw("SAMEP");
  
  
  TCanvas *c4 = new TCanvas ("c4", "c4", 1400,600);
  c4->Divide(2,1);
  
  c4->cd(1)->SetLogy();
  GHMUr1->Draw("AP");
  gMU->Draw("SAMEP");

  c4->cd(2)->SetLogy();
  GHMUr2->Draw("AP");
  gRMU->Draw("SAMEP");  
  

  
  
  cout<<"Everything Done!"<<endl;
}






//---------------//
//   FUNCTIONS   //
//---------------//

void formatTH1D(TH1D *h, int col, int lwidth, const char *xtitle, const char *ytitle)
{
  h->SetLineColor(col);
  h->SetLineWidth(lwidth);
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  
}

void formatTGraph(TGraph *g, int m, double ms, int mc, const char *xtitle, const char *ytitle)
{
  g->SetMarkerStyle(m);
  g->SetMarkerSize(ms);
  g->SetMarkerColor(mc);
  g->GetXaxis()->SetTitle(xtitle);
  g->GetYaxis()->SetTitle(ytitle);
  
}

Double_t Round(double X){
  
  Double_t E = floor(X*10 +0.5);
  return E/10;
  
}


//spectrum power law function
Double_t PHI(Double_t *x, Double_t *par){
  Double_t E = TMath::Power(10,x[0]);
  Double_t A = TMath::Power(10,18);
  Double_t gamma = 2.7;
  
  return TMath::Log(10) * A * TMath::Power(E,-gamma+1);
}

//Number of Muons function
Double_t MUONS(Double_t *x, Double_t *par){
  Double_t E = TMath::Power(10,x[0]);
  Double_t B = TMath::Power(10,-10);
  Double_t beta1 = 0.90;
  Double_t beta2 = par[0];
  
  Double_t logENLIM = 18.8;
  Double_t ENLIM = TMath::Power(10,logENLIM);
   
  Double_t ret1=B*TMath::Power(ENLIM,beta1);
  Double_t ret2=B*TMath::Power(ENLIM,beta2);
  
  if(x[0]<logENLIM)
    return B*TMath::Power(E,beta1);
  else
    return B*TMath::Power(E,beta2)/ret2*ret1;
}


//Sigma Muons function
Double_t SIGMAMUONS(Double_t *x, Double_t *par){
  Double_t E = x[0];
  Double_t C = -16.;
  Double_t delta = 1.6;
  
  return TMath::Exp(C + delta*E);
}

//Energy Functions
Double_t EnergyFunction(Double_t Et, Double_t Er){

  Double_t res  = 0; 
  //Double_t ET = TMath::Power(10,Et);
  //Double_t ER = TMath::Power(10,Er);
  
  Double_t ET = TMath::Log10(TMath::Power(10,Et));
  Double_t ER = TMath::Log10(TMath::Power(10,Er));
  //Double_t expET = TMath::Exp(TMath::Power(10,Et));
  Double_t powET = TMath::Power(10,Et);
  
  Double_t A = TMath::Power(10,18)*TMath::Power(powET,-2.7+1);
  
  res = A/powET * TMath::Exp( -TMath::Power(((ET-ER)/(TMath::Log10(1.2))),2)/2  );
  
  return res;
}

Double_t EFIT(Double_t *x, Double_t *par){

  Double_t Energy  = 0;
  
  Double_t Et = x[0];
  Double_t Er = par[0];
  Double_t Norm = par[1];
  Double_t ET = TMath::Power(10,Et);
  
  Energy = EnergyFunction(Et, Er);
  
  return  Norm * TMath::Log(10) * ET * Energy;
}

Double_t FEFI(Double_t *x, Double_t *par){
  
  Double_t Energy  = 0;
  
  Double_t Et = x[0];
  Double_t Er = par[0];
  Double_t ET = TMath::Power(10,Et);
  
  Energy = EnergyFunction(Et, Er);
  
  return  TMath::Log(10) * ET * Energy;
}

Double_t FE(Double_t *x, Double_t *par){
    
  Double_t Energy  = 0;
  
  Double_t Et = x[0];
  Double_t Er = par[0];
  Double_t ET = TMath::Power(10,Et);
  
  Energy = EnergyFunction(Et, Er);
  
  return  TMath::Log(10) * ET * ET * Energy;
}

Double_t FE2(Double_t *x, Double_t *par){
    
  Double_t Energy  = 0;
  
  Double_t Et = x[0];
  Double_t Er = par[0];
  Double_t ET = TMath::Power(10,Et);
  
  Energy = EnergyFunction(Et, Er);
  
  return  TMath::Log(10) * ET * ET * ET * Energy;
}
//Sigma Muons function
Double_t SIGMAMUONS(Double_t *x, Double_t *par){
  Double_t E = x[0];
  Double_t C = -16.;
  Double_t delta = 1.6;
  
  return TMath::Exp(C + delta*E);
}
//Muons Functions
Double_t MuonFunction (Double_t beta, Double_t  Et, Double_t Er){
  
  Double_t Res  = 0;
  Double_t B = TMath::Power(10,-10);
  
  Double_t MU = TMath::Power(10,Et);
  Double_t powNT = TMath::Power(TMath::Power(10,Et)/B,1/beta);
  //Double_t NR = TMath::Power(TMath::Power(10,Er)/B,1/beta);
  
  Double_t NT = TMath::Log10(TMath::Power(TMath::Power(10,Et)/B,1/beta));
  Double_t NR = TMath::Log10(TMath::Power(TMath::Power(10,Er)/B,1/beta));
  
  Double_t A  = TMath::Power(10,18) * TMath::Power(powNT,-2.7+1);
  
  //Double_t C = -16.;
  //Double_t delta = 1.6;
  //Double_t Sigma = TMath::Power(10,TMath::Exp(C + delta*Et));
  
  Res = A/powNT * TMath::Exp( -TMath::Power(((NT-NR)/(TMath::Log10(1.2))),2)/2  );
  
  
  return Res;
}

Double_t MUFIT(Double_t *x, Double_t *par){
  
  Double_t Result  = 0;
  Double_t beta = 0.90;
  
  Double_t ET = x[0];
  Double_t ER = par[0];
  Double_t Norm = par[1];
  Double_t MU = TMath::Power(10,ET);
  
  Result = MuonFunction(beta, ET, ER);
  
  return  Norm * Result * TMath::Log(10) * MU;
}


Double_t FMU(Double_t *x, Double_t *par){
  
  Double_t Result  = 0;
  Double_t beta = 0.90;
  
  Double_t ET = x[0];
  Double_t ER = par[0];
  Double_t MU = TMath::Power(10,ET);
  
  Result = MuonFunction(beta, ET, ER);
  
  return  Result * TMath::Log(10) * MU * MU;
}

Double_t FMU2(Double_t *x, Double_t *par){
  
  Double_t Result  = 0;
  Double_t beta = 0.90;
  
  Double_t ET = x[0];
  Double_t ER = par[0];
  Double_t MU = TMath::Power(10,ET);
  
  Result = MuonFunction(beta, ET, ER);
  
  return  Result * TMath::Log(10) * MU * MU * MU;
}

Double_t FMUFI(Double_t *x, Double_t *par){
  
  Double_t Result  = 0;
  Double_t beta = 0.90;
  
  Double_t ET = x[0];
  Double_t ER = par[0];
  Double_t MU = TMath::Power(10,ET);
  
  Result = MuonFunction(beta, ET, ER);
  
  return  Result * TMath::Log(10) * MU;
}

 
