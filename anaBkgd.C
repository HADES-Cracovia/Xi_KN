#include "TLatex.h"

void anaBkgd(){
    THStack *hLmass = new THStack(); //Lambda before cuts and scaling
    THStack *hXmass = new THStack(); //Xi before cuts and scaling
    THStack *hLmass_sc = new THStack(); //Lambda scaled with cr_sec, no cuts
    THStack *hLmass_dist_sc = new THStack(); //Lambda scaled with cr_sec, with dist cut
    THStack *hXmass_sc = new THStack(); //Xi, scaled with cr_sec, no cuts
    THStack *hXmass_dist_sc = new THStack(); //Xi, cut dist, scaled with cr_sec
    THStack *hXmass_vert_sc = new THStack(); //Xi, cut vertex_z, scaled with cr_sec
    THStack *hXmass_sum_bkgd = new THStack(); //Xi and sum of bkgd channels

    TLegend *lL = new TLegend(.15, .3, .5, .95);
    TLegend *lX = new TLegend(.15, .6, .25, .85);
    TLegend *lLsc = new TLegend(.15, .6, .25, .85);
    TLegend *lXdistsc = new TLegend(.15, .5, .25, .85);
    TLegend *lXvertsc = new TLegend(.15, .5, .25, .85);
    TLegend *lXsumbkgd = new TLegend(.15, .75, .25, .85);
    lL -> SetFillStyle(0);
    lL -> SetBorderSize(0);
    lL -> SetTextSize(.1);
    lX -> SetFillStyle(0);
    lX -> SetBorderSize(0);
    lX -> SetTextSize(.05);
    lLsc -> SetFillStyle(0);
    lLsc -> SetBorderSize(0);
    lLsc -> SetTextSize(.05);
    lXdistsc -> SetFillStyle(0);
    lXdistsc -> SetBorderSize(0);
    lXdistsc -> SetTextSize(.05);
    lXvertsc -> SetFillStyle(0);
    lXvertsc -> SetBorderSize(0);
    lXvertsc -> SetTextSize(.05);
    lXsumbkgd -> SetFillStyle(0);
    lXsumbkgd -> SetBorderSize(0);
    lXsumbkgd -> SetTextSize(.04);

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
    
    gStyle -> SetTitleOffset(1.3, "xy");
    gStyle -> SetTitleSize(0.04,"xy");
    gStyle->SetLabelSize(0.04,"xy");
      
    char chanNo[64];
    int chan[] = {90, 10, 21, 22, 23, 62, 100};
    Double_t cr_sec[] = {4.8, 600., 100., 30., 30., 20., 20.};
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
    c2 -> Divide(3,2);
    
    TH1F *hXi_sum_bkgd = new TH1F("hXi_sum_bkgd", "Reconstructed Xi- and sum of all bkgd channels", 200, 1200, 1400); //Xi with bkgd
    // TH1F *clXi_sum_bkgd = new TH1F("clXi_sum_bkgd", "Reconstructed Xi- and sum of all bkgd channels (scaled)", 200, 1200, 1400); //Xi with bkgd after cesec scaling

    TF1 *fit = new TF1("fit", "gaus");
    double center_fit, sigma_fit, a, b, s_b;
    double cnt_peak_XiDistSc, cnt_peak_sumBkgd;

    for(int i = 0; i < 7; i++){
	sprintf(chanNo, "output_%03d_all.root", chan[i]);
//	cout << chanNo << endl;
      	TFile *f1 = TFile::Open(chanNo, "READ");
	
	TH1F *hLambda = (TH1F*)f1->Get("hMLAll"); //Lambda before cuts
	TH1F *hLambdaDist = (TH1F*)f1->Get("hMLDist"); //Lambda after dist cut
        TH1F *hXi = (TH1F*)f1->Get("hKmassall"); //Xi before cuts
	TH1F *hXiDist = (TH1F*)f1->Get("hKmassdist"); //Xi after dist cut
	TH1F *hXiVert = (TH1F*)f1->Get("hKmassvert"); //Xi after vertex_z > 50mm cut
	
	TH1F *clLambda, *clLambdaSc, *clLambdaDistSc, *clXi, *clXiSc, *clXiDistSc,  *clXiVertSc;

	int col = i+1;
	if(i == 2) col = 8;
	if(i == 4) col = 42;
	if(i == 6) col = 9;
	hLambda -> SetLineColor(col);
	hLambdaDist -> SetLineColor(col);
	hXi -> SetLineColor(col);
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

	TPaveText *ptch = new TPaveText(.2, .8, .3, .9, "NDC");
	ptch -> SetFillColor(0);
	ptch -> SetBorderSize(0);
	ptch -> SetTextSize(.07);
	ptch -> AddText(reac[i].c_str());
	
	c1 -> cd(i+1);
	hLambda -> Draw();
	ptch -> Draw("same");
	if(i == 0) pt1 -> Draw();
	c1 -> cd(i+8);
	hXi -> Draw();
	if(i == 0) pt2 -> Draw();
	c1 -> cd(i+15);
	hXiDist -> Draw();
	if(i == 0) pt3 -> Draw();

	clLambda = (TH1F*)hLambda -> Clone();
	clLambdaSc = (TH1F*)hLambdaDist -> Clone();
	clLambdaDistSc = (TH1F*)hLambdaDist -> Clone();
	clXi = (TH1F*)hXi -> Clone();
	clXiSc = (TH1F*)hXiDist -> Clone();
	clXiDistSc = (TH1F*)hXiDist -> Clone();
	clXiVertSc = (TH1F*)hXiVert -> Clone();

	clLambdaSc -> Scale(cr_sec[i]);
	clLambdaDistSc -> Scale(cr_sec[i]);
	clXiSc -> Scale(cr_sec[i]);
	clXiDistSc -> Scale(cr_sec[i]);
	clXiVertSc -> Scale(cr_sec[i]);

	clXiDistSc -> Fit(fit);
	center_fit = fit -> GetParameter(1);
	sigma_fit = fit -> GetParameter(2);
	a = center_fit - 3*sigma_fit;
	b = center_fit + 3*sigma_fit;

	if(i == 0){
	    cnt_peak_XiDistSc = clXiDistSc -> Integral(a,b);
	    s_b = 1;
	}
	else{
	    cnt_peak_sumBkgd = clXiDistSc -> Integral(a,b);
	    s_b = cnt_peak_XiDistSc/cnt_peak_sumBkgd;
	}
		
	int cnt = hXi -> Integral();
	int cnt_scale = clXiSc -> Integral();
	int cnt_dist_scale = clXiDistSc -> Integral();
	int cnt_vert_scale = clXiVertSc -> Integral();
	
	printf("chann: %03d, cr_sec: %.02f, counts: %d, scale: %d, dist_scale: %d, vert_scale: %d, S/B: %f\n", chan[i], cr_sec[i], cnt, cnt_scale, cnt_dist_scale, cnt_vert_scale, s_b);

	hLmass -> Add(hLambda);
	hXmass -> Add(hXi);
	hLmass_sc -> Add(clLambdaSc);
	hLmass_dist_sc -> Add(clLambdaDistSc);
	hXmass_sc -> Add(clXiSc);
	hXmass_dist_sc -> Add(clXiDistSc);
	hXmass_vert_sc -> Add(clXiVertSc);
	lL -> AddEntry(clLambda, reac[i].c_str(), "L");
	lX -> AddEntry(clXi, reac[i].c_str(), "L");	
	lLsc -> AddEntry(clLambdaSc, reac[i].c_str(), "L");
	lXdistsc -> AddEntry(clXiDistSc, reac[i].c_str(), "L");
	lXvertsc -> AddEntry(clXiVertSc, reac[i].c_str(), "L");

	if(i == 0){
	    hXmass_sum_bkgd -> Add(clXiDistSc);
	    lXsumbkgd -> AddEntry(clXiDistSc);
	}
	else
	    hXi_sum_bkgd -> Add(clXiDistSc);
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
    
    TCanvas *cX_dist_sc = new TCanvas("cXi_dist_sc", "Reconstruction of #Xi^{-} (dist cut, crsec scaling)", 1000, 600);
    cX_dist_sc -> cd();
    hXmass_dist_sc -> Draw("nostack");
    hXmass_dist_sc -> SetTitle("#Xi^{-} mass reconstruction");
    hXmass_dist_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_dist_sc -> GetYaxis() -> SetTitle("counts");
    cX_dist_sc -> Modified();
    lXdistsc -> Draw("same");

    TCanvas *cX_vert_sc = new TCanvas("cXi_vert_sc", "Reconstruction of #Xi^{-} (vertex_z cut, crsec scaling)", 1000, 600);
    cX_vert_sc -> cd();
    hXmass_vert_sc -> Draw("nostack");
    hXmass_vert_sc -> SetTitle("#Xi^{-} mass reconstruction");
    hXmass_vert_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_vert_sc -> GetYaxis() -> SetTitle("counts");
    cX_vert_sc -> Modified();
    lXvertsc -> Draw("same");

    c2 -> cd(1);
    hLmass -> Draw("nostack");
    hLmass -> SetTitle("#Lambda\(1115\), no cuts, no scaling");
    hLmass -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    //lL -> Draw("same");

    c2 -> cd(2);
    hXmass_sc -> Draw("nostack");
    hXmass_sc -> SetTitle("#Xi^{-}, #Lambda mass cut, with scaling");
    hXmass_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    //lX -> Draw("same");

    c2 -> cd(4);
    hLmass_dist_sc -> Draw("nostack");
    hLmass_dist_sc -> SetTitle("#Lambda\(1115\), dist cut, crsec scaling");
    hLmass_dist_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hLmass_dist_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    //lLsc -> Draw("same");

    c2 -> cd(5);
    hXmass_dist_sc -> Draw("nostack");
    hXmass_dist_sc -> SetTitle("#Xi^{-}, #Lambda mass & dist cut, crsec scaling");
    hXmass_dist_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_dist_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();
    //lXdistsc -> Draw("same");

    /*c2 -> cd(6);
    hXmass_vert_sc -> Draw("nostack");
    hXmass_vert_sc -> SetTitle("#Xi^{-}, vertex_z cut, crsec scaling");
    hXmass_vert_sc -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_vert_sc -> GetYaxis() -> SetTitle("counts");
    c2 -> Modified();*/
    //lXdistsc -> Draw("same");

    c2 -> cd(3);
    lL -> Draw();


    cnt_peak_sumBkgd = hXi_sum_bkgd -> Integral(a,b);
    
    
    TCanvas *cXi_sum_bkgd = new TCanvas("cXi_sum_bkgd", "Reconstruction of #Xi^{-} & sum of all bkgd channels", 1000, 600);
    hXmass_sum_bkgd -> Add(hXi_sum_bkgd);
    hXi_sum_bkgd -> SetLineColor(4);;
    lXsumbkgd -> AddEntry(hXi_sum_bkgd);
    cXi_sum_bkgd -> cd();
    hXmass_sum_bkgd -> Draw("nostack");
    hXmass_sum_bkgd -> SetTitle("#Xi^{-} mass reconstruction & bkgd");
    hXmass_sum_bkgd -> GetXaxis() -> SetTitle("Mass [MeV]");
    hXmass_sum_bkgd -> GetYaxis() -> SetTitle("counts");
    lXsumbkgd -> Draw("same");
    cXi_sum_bkgd -> Modified();
    lXsumbkgd -> Draw("same");
    
    TFile *f = new TFile("out_anaBkgd_test.root", "RECREATE");
    cX_dist_sc -> Write();
    cX_vert_sc -> Write();
    cXi_sum_bkgd -> Write();
    c1 -> Write();
    c2 -> Write();
    f -> Close();

    cX_dist_sc -> SaveAs("xiMass_dist.eps");
    cX_vert_sc -> SaveAs("xiMass_vert.eps");
    cXi_sum_bkgd ->  SaveAs("xiMass_sumbkgd.eps");
    c1 -> SaveAs("massSpec.eps");
    c2 -> SaveAs("xiLRec.eps");
}
