#include "TLatex.h"

void anaBkgd(){
    THStack *hLmass = new THStack(); //Lambda before cuts and scaling
    THStack *hXmass = new THStack(); //Xi before cuts and scaling
    THStack *hLmass_sc = new THStack(); //Lambda scaled with cr_sec
    THStack *hXmass_dist_sc = new THStack(); //Xi, cut dist, scaled with cr_sec
    THStack *hXmass_vert_sc = new THStack(); //Xi, cut vertex_z, scaled with cr_sec

    TLegend *lL = new TLegend(.15, .6, .25, .85);
    TLegend *lX = new TLegend(.15, .6, .25, .85);
    TLegend *lLsc = new TLegend(.15, .6, .25, .85);
    TLegend *lXdistsc = new TLegend(.15, .6, .25, .85);
    TLegend *lXvertsc = new TLegend(.15, .6, .25, .85);
    lL -> SetFillStyle(0);
    lL -> SetBorderSize(0);
    lL -> SetTextSize(.04);
    lX -> SetFillStyle(0);
    lX -> SetBorderSize(0);
    lX -> SetTextSize(.04);
    lLsc -> SetFillStyle(0);
    lLsc -> SetBorderSize(0);
    lLsc -> SetTextSize(.04);
    lXdistsc -> SetFillStyle(0);
    lXdistsc -> SetBorderSize(0);
    lXdistsc -> SetTextSize(.04);
    lXvertsc -> SetFillStyle(0);
    lXvertsc -> SetBorderSize(0);
    lXvertsc -> SetTextSize(.04);

    TPaveText *pt1 = new TPaveText(.2, .7, .3, .85, "NDC");
    pt1 -> SetFillColor(0);
    pt1 -> SetBorderSize(0);
    pt1 -> SetTextSize(.07);
    TPaveText *pt2 = new TPaveText(.2, .7, .3, .85, "NDC");
    pt2 -> SetFillColor(0);
    pt2 -> SetBorderSize(0);
    pt2 -> SetTextSize(.07);
    TPaveText *pt3 = new TPaveText(.2, .7, .3, .85, "NDC");
    pt3 -> SetFillColor(0);
    pt3 -> SetBorderSize(0);
    pt3 -> SetTextSize(.07);
    TPaveText *pt4 = new TPaveText(.2, .7, .3, .85, "NDC");
    pt4 -> SetFillColor(0);
    pt4 -> SetBorderSize(0);
    pt4 -> SetTextSize(.07);
    
    char chanNo[64];
    int chan[] = {90, 10, 21, 22, 23, 62, 100};
    Double_t cr_sec[] = {4.8, 600., 100., 30., 20., 30., 20.};
    string reac[] = {"#Xi^{-}K^{+}K^{+}p",
		     "pp2#pi^{+}2#pi^{-}",
		     "p#LambdaK^{0}_{s}#pi^{+}",
		     "p#LambdaK^{+}#pi^{+}#pi^{-}",
		     "n#LambdaK^{0}_{s}2#pi^{+}",
		     "p#Sigma^{0}K^{0}_{s}#pi^{+}",
		     "ppK^{0}_{s}K^{0}_{s}"}
    
    TCanvas *c1 = new TCanvas("c1", "All channels - Lambda(1115) & Ksi mass", 2200, 1000);
    c1 -> Divide(7,3);
    gStyle -> SetOptStat(0);
    TCanvas *c2 = new TCanvas("c2", "Xi & Lambda mass --- before and after cuts and scaling", 2200, 1000);
    c2 -> Divide(2,2);
    
    for(int i = 0; i < 7; i++){
	sprintf(chanNo, "output_%03d.root", chan[i]);
//	cout << chanNo << endl;
      	TFile *f1 = TFile::Open(chanNo, "READ");
	
	TH1F *hLambda = (TH1F*)f1->Get("hMLAll"); //Lambda beefore cuts
        TH1F *hXi = (TH1F*)f1->Get("hKmassall"); //Xi before cuts
	TH1F *hLambdaDist = (TH1F*)f1->Get("hMLDist"); //Lambda after dist cut
	TH1F *hXiDist = (TH1F*)f1->Get("hKmassdist"); //Xi after dist cut
	TH1F *hXiVert = (TH1F*)f1->Get("hKmassvert"); //Xi after vertex_z > 50mm cut
	
	TH1F *clLambda, *clXi, *clLambdaDist, *clXiDist, *clXiVert;

	int col = i+1;
	if(i == 2) col = 8;
	if(i == 4) col = 42;
	if(i == 6) col = 9;
	hLambda -> SetLineColor(col);
	hXi -> SetLineColor(col);
	hLambdaDist -> SetLineColor(col);
	hXiDist -> SetLineColor(col);
	hXiVert -> SetLineColor(col);
	if(i == 0){
	    hLambda -> SetLineWidth(2);	
	    hXi -> SetLineWidth(2);
	    hLambdaDist -> SetLineWidth(2);	
	    hXiDist -> SetLineWidth(2);
	    hXiVert -> SetLineWidth(2);
	  	    
	    pt1 -> AddText("#Lambda\(1115\)");
	    pt2 -> AddText("#Xi^{-} no cuts");
	    pt3 -> AddText("#Xi^{-} dist cut");
	    pt4 -> AddText("#Xi^{-} vertex_z cut");
	}

	hLambda -> GetXaxis() -> SetLabelSize(.08);	
	hXi -> GetXaxis() -> SetLabelSize(.08);
	hXiDist -> GetXaxis() -> SetLabelSize(.08);
	hXiVert -> GetXaxis() -> SetLabelSize(.08);
	hLambda -> GetYaxis() -> SetLabelSize(.08);	
	hXi -> GetYaxis() -> SetLabelSize(.08);
	hXiDist -> GetYaxis() -> SetLabelSize(.08);
	hXiVert -> GetYaxis() -> SetLabelSize(.08);
	hLambda -> GetXaxis() -> SetNdivisions(5);
	hXi -> GetXaxis() -> SetNdivisions(5);
	hXiDist -> GetXaxis() -> SetNdivisions(5);
	hXiVert -> GetXaxis() -> SetNdivisions(5);
	    
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
	clLambdaDist = (TH1F*)hLambdaDist -> Clone();
	clXiDist = (TH1F*)hXiDist -> Clone();
	clXiVert = (TH1F*)hXiDist -> Clone();
	
	clLambdaDist -> Scale(cr_sec[i]);
	clXiDist -> Scale(cr_sec[i]);
	clXiVert -> Scale(cr_sec[i]);

	int cnt = hXi -> Integral();
	int cnt_dist = hXiDist -> Integral();
	int cnt_dist_scale = clXiDist -> Integral();
	int cnt_vert = hXiVert -> Integral();
	int cnt_vert_scale = clXiVert -> Integral();
	
	printf("chann: %03d, cr_sec: %.02f, counts: %d, dist: %d, dist_scale: %d, vert: %d, vert_scale: %d\n", chan[i], cr_sec[i], cnt, cnt_dist, cnt_scale, cnt_vert, cnt_vert);

	hLmass -> Add(hLambda);
	hXmass -> Add(hXi);
	hLmass_sc -> Add(clLambdaDist);
	hXmass_dist_sc -> Add(clXiDist);
	hXmass_vert_sc -> Add(clXiVert);
	lL -> AddEntry(clLambda, reac[i].c_str(), "L");
	lX -> AddEntry(clXi, reac[i].c_str(), "L");	
	lLsc -> AddEntry(clLambdaDist, reac[i].c_str(), "L");
	lXdistsc -> AddEntry(clXiDist, reac[i].c_str(), "L");
	lXvertsc -> AddEntry(clXiVert, reac[i].c_str(), "L");	
    }

/*    TCanvas *cL = new TCanvas("cLambda", "Reconstruction of #Lambda(1115)");
    cL -> cd();
    hLmass -> Draw("nostack");
    hLmass -> SetTitle("#Lambda1115 mass reconstruction");
    hLmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass -> GetYaxis() -> SetTitle("counts");
    cL -> Modified();
    lL -> Draw("same");
*/

    TCanvas *cX_dist = new TCanvas("cXi_dist", "Reconstruction of #Xi^{-} (dist cut)", 1000, 600);
    cX_dist -> cd();
    hXmass_dist_sc -> Draw("nostack");
    hXmass_dist_sc -> SetTitle("#Xi^{-} mass reconstruction");
    hXmass_dist_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_dist_sc -> GetYaxis() -> SetTitle("counts");
    cX_dist -> Modified();
    lXdistsc -> Draw("same");

    TCanvas *cX_vert = new TCanvas("cXi_vert", "Reconstruction of #Xi^{-} (vertex_z cut)", 1000, 600);
    cX_vert -> cd();
    hXmass_vert_sc -> Draw("nostack");
    hXmass_vert_sc -> SetTitle("#Xi^{-} mass reconstruction");
    hXmass_vert_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_vert_sc -> GetYaxis() -> SetTitle("counts");
    cX_vert -> Modified();
    lXvertsc -> Draw("same");

    c2 -> cd(1);
    hLmass -> Draw("nostack");
    hLmass -> SetTitle("#Lambda\(1115\), no cuts, no scaling");
    hLmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    lL -> Draw("same");

    c2 -> cd(2);
    hXmass -> Draw("nostack");
    hXmass -> SetTitle("#Xi^{-}, no cuts, no scaling");
    hXmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    lX -> Draw("same");

    c2 -> cd(3);
    hLmass_sc -> Draw("nostack");
    hLmass_sc -> SetTitle("#Lambda\(1115\), dist cut, crsec scaling");
    hLmass_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    lLsc -> Draw("same");

    c2 -> cd(4);
    hXmass_sc -> Draw("nostack");
    hXmass_sc -> SetTitle("#Xi^{-},dist cut, crsec scaling");
    hXmass_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    lXsc -> Draw("same");

    TFile *f = new TFile("out_anaBkgd.root", "RECREATE");
    cX_dist -> Write();
    cX_vert -> Write();
    c1 -> Write();
    c2 -> Write();
    f -> Close();
/*
    cX -> SaveAs("xiMass.eps");
    c1 -> SaveAs("massSpec.eps");
    c2 -> SaveAs("xiLRec.eps");*/
}
