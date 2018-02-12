#include "TLatex.h"

void anaBkgd(){
    THStack *hLmass = new THStack();
    THStack *hXmass = new THStack();

    TLegend *lL = new TLegend(.15, .6, .5, .85);
    TLegend *lX = new TLegend(.15, .6, .5, .85);
    lL -> SetFillColor(kWhite);
    lL -> SetBorderSize(0);
    lX -> SetFillColor(kWhite);
    lX -> SetBorderSize(0);
    
    TPaveText *pt1 = new TPaveText(.6, .7, .9, .85, "NDC");
    pt1 -> SetFillColor(0);
    pt1 -> SetBorderSize(0);
    pt1 -> SetTextSize(.07);
    TPaveText *pt2 = new TPaveText(.6, .7, .9, .85, "NDC");
    pt2 -> SetFillColor(0);
    pt2 -> SetBorderSize(0);
    pt2 -> SetTextSize(.07);
    TPaveText *pt3 = new TPaveText(.6, .7, .9, .85, "NDC");
    pt3 -> SetFillColor(0);
    pt3 -> SetBorderSize(0);
    pt3 -> SetTextSize(.07);
    
    char chanNo[64];
    int chan[] = {90, 10, 21, 22, 23, 62, 100};
    Double_t cr_sec[] = {4.8, 600., 100., 30., 20., 30., 20.};
    string reac[] = {"signal", "pp2#pi^{+}2#pi^{-}", "p#LambdaK^{0}_{s}#pi^{+}", "n#LambdaK^{0}_{s}2#pi^{+}", "p#Sigma^{0}K^{0}_{s}#pi^{+}", "p#LambdaK^{+}#pi^{+}#pi^{-}", "ppK^{0}_{s}K^{0}_{s}"}
    
    TCanvas *c1 = new TCanvas("c1", "All channels - Lambda(1115) & Ksi mass");
    c1 -> Divide(7,3);
    gStyle -> SetOptStat(0);
    
    for(int i = 0; i < 7; i++){
	sprintf(chanNo, "output_%03d.root", chan[i]);
//	cout << chanNo << endl;
      	TFile *f1 = TFile::Open(chanNo, "READ");
	
	TH1F *hLambda = (TH1F*)f1->Get("hMLAll");
        TH1F *hXi = (TH1F*)f1->Get("hKmassall");
	TH1F *hXiDist = (TH1F*)f1->Get("hKmassdist");
	TH1F *clLambda, *clXi, *clXiDist;

	hXi -> GetXaxis() -> SetRangeUser(1200,1400);
	hLambda -> SetLineColor(i+1);
	hXi -> SetLineColor(i+1);
	hXiDist -> SetLineColor(i+1);
	if(i == 0){
	    hLambda -> SetLineWidth(2);	
	    hXi -> SetLineWidth(2);
	    hXiDist -> SetLineWidth(2);

	    pt1 -> AddText("#Lambda\(1115\)");
	    pt2 -> AddText("#Xi^{-} no cuts");
	    pt3 -> AddText("#Xi^{-} dist cut");
	}

	c1 -> cd(i+1);
	hLambda -> Draw();
	if(i == 0) pt1 -> Draw();
	c1 -> cd(i+8);
	hXi -> Draw();
	if(i == 0) pt2 -> Draw();
	c1 -> cd(i+15);
	hXiDist -> Draw();
	if(i == 0) pt3 -> Draw();

	
	clLambda = (TH1F*)hLambda -> Clone();
	clXi = (TH1F*)hXi -> Clone();
	clXiDist = (TH1F*)hXiDist -> Clone();
	
	clLambda -> Scale(cr_sec[i]);
	clXi -> Scale(cr_sec[i]);
	clXiDist -> Scale(cr_sec[i]);

	int cnt = clXiDist -> Integral();
	printf("chann: %03d, counts: %d\n", chan[i], cnt);

	hLmass -> Add(clLambda);
	hXmass -> Add(clXiDist);
	lL -> AddEntry(clLambda, reac[i].c_str(), "L");
	lX -> AddEntry(clXi, reac[i].c_str(), "L");	
    }

    TCanvas *cL = new TCanvas("cLambda", "Reconstruction of #Lambda(1115)");
    cL -> cd();
    hLmass -> Draw();
    hLmass -> SetTitle("#Lambda1115 mass reconstruction");
    hLmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass -> GetYaxis() -> SetTitle("counts");
    cL -> Modified();
    lL -> Draw("same");

    TCanvas *cX = new TCanvas("cXi", "Reconstruction of #Xi^{-}");
    cX -> cd();
    hXmass -> Draw();
    hXmass -> SetTitle("#Xi^{-} mass reconstruction");
    hXmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass -> GetYaxis() -> SetTitle("counts");
    cX -> Modified();
    lX -> Draw("same");

}
