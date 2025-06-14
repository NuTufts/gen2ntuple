import os,sys
from math import fabs
import ROOT as rt

rt.gStyle.SetOptStat(0)

filepath="flashprediction_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"

inputfile = rt.TFile(filepath)
tree = inputfile.Get("FlashPredictionTree")

out = rt.TFile("temp.root","recreate")

h_sinkdiv_v_dist2true = rt.TH2D("h_sinkdiv_v_dist2true","",100,0,100,100,0,100)

h_sinkdiv_v_fracerr_good = rt.TH2D("h_sinkdiv_v_fracerr_good",";fractional error;sinkhorn divergence",100,-2.0,2.0,100,0,100)
h_sinkdiv_v_fracerr_bad  = rt.TH2D("h_sinkdiv_v_fracerr_bad", ";fractional error;sinkhorn divergence",100,-2.0,2.0,100,0,100)
h_sinkdiv_v_fracerr_good.SetLineColor(rt.kRed)
h_sinkdiv_v_fracerr_bad.SetLineColor(rt.kBlack)

hsinkdiv_good = rt.TH1D("hsinkdiv_good",";sinkhorn divergence",100,0,500)
hsinkdiv_bad  = rt.TH1D("hsinkdiv_bad", ";sinkhorn divergence",100,0,500)
hsinkdiv_good.SetLineColor(rt.kRed)
hsinkdiv_bad.SetLineColor(rt.kBlack)

hfracerr_good = rt.TH1D("hfracerr_good",";fractional error",100,-2.0,2.0)
hfracerr_bad  = rt.TH1D("hfracerr_bad", ";fractional error",100,-2.0,2.0)
hfracerr_good.SetLineColor(rt.kRed)
hfracerr_bad.SetLineColor(rt.kBlack)

hobspe_good = rt.TH1D("hobspe_good",";observed pe total",1000,0.0,10000.0)
hobspe_bad  = rt.TH1D("hobspe_bad", ";observed pe total",1000,0.0,10000.0)
hobspe_good.SetLineColor(rt.kRed)
hobspe_bad.SetLineColor(rt.kBlack)

h_predpe_v_obspe_good = rt.TH2D("h_predpe_v_obspe_good",";observed flash pe;predicted flashe pe",100,0,1.0e4,100,0,1.0e4)
h_predpe_v_obspe_bad  = rt.TH2D("h_predpe_v_obspe_bad", ";observed flash pe;predicted flashe pe",100,0,1.0e4,100,0,1.0e4)
h_predpe_v_obspe_good.SetLineColor(rt.kRed)
h_predpe_v_obspe_bad.SetLineColor(rt.kBlack)

nentries = tree.GetEntries()

for ientry in range(nentries):

    if ientry>0 and ientry%10000==0:
        print("entry ",ientry)
        
    tree.GetEntry(ientry)

    obs_total_pe = tree.obs_total_pe
    nvertices = tree.pred_total_pe_all.size()

    if obs_total_pe<5.0:
        continue

    for ivtx in range(nvertices):
        dist2true = tree.vtx_dist_to_true.at(ivtx)
        sinkdiv   = tree.sinkhorn_div_all.at(ivtx)

        sinkdiv_reg1 = sinkdiv[1]
        pred_total_pe_all = tree.pred_total_pe_all.at(ivtx)*2.0
        
        fracerr = (pred_total_pe_all-obs_total_pe)/(0.1+obs_total_pe)

        h_sinkdiv_v_dist2true.Fill( dist2true, sinkdiv_reg1 )
        if dist2true<3.0:
            h_sinkdiv_v_fracerr_good.Fill( fracerr, sinkdiv_reg1 )
            hfracerr_good.Fill( fracerr )
            hsinkdiv_good.Fill( sinkdiv_reg1 )
            h_predpe_v_obspe_good.Fill( obs_total_pe, pred_total_pe_all )
            hobspe_good.Fill( obs_total_pe )
        else:
            h_sinkdiv_v_fracerr_bad.Fill( fracerr, sinkdiv_reg1 )
            hfracerr_bad.Fill( fracerr )
            hsinkdiv_bad.Fill( sinkdiv_reg1 )
            h_predpe_v_obspe_bad.Fill( obs_total_pe, pred_total_pe_all )
            hobspe_bad.Fill( obs_total_pe )

       

csinkdiv_v_dist2true = rt.TCanvas("csinkdiv_v_dist2true","Sinkhorn Divergence vs. Distance to True Vertex",1200,1000)
h_sinkdiv_v_dist2true.Draw("colz")
csinkdiv_v_dist2true.SetGridx(1)
csinkdiv_v_dist2true.SetGridy(1)
csinkdiv_v_dist2true.Update()

cpredpe_v_obspe = rt.TCanvas("cpredpe_vs_obspe","Predicted PE vs. Observed PE",1800,800)
cpredpe_v_obspe.Divide(2,1)
cpredpe_v_obspe.cd(1)
cpredpe_v_obspe.cd(1).SetGridx(1)
cpredpe_v_obspe.cd(1).SetGridy(1)
h_predpe_v_obspe_good.Draw("colz")
cpredpe_v_obspe.cd(2)
cpredpe_v_obspe.cd(2).SetGridx(1)
cpredpe_v_obspe.cd(2).SetGridy(1)
h_predpe_v_obspe_bad.Draw("colz")
cpredpe_v_obspe.Update()
cpredpe_v_obspe.SaveAs("cpredpe_vs_obspe.png")

csinkdiv_v_fracerr = rt.TCanvas("csinkdiv_v_fracerr","Sinkhorn Divergence vs. Frac Err",1200,1000)
h_sinkdiv_v_fracerr_bad.Draw("box")
h_sinkdiv_v_fracerr_good.Draw("boxsame")
csinkdiv_v_fracerr.SetGridx(1)
csinkdiv_v_fracerr.SetGridy(1)
csinkdiv_v_fracerr.SetLogz(1)
csinkdiv_v_fracerr.Update()
csinkdiv_v_fracerr.SaveAs("csinkdiv_v_fracerr.png")

csinkdiv = rt.TCanvas("csinkdiv","Sinkhorn Divergence",1000,800)
hsinkdiv_bad.Draw("hist")
hsinkdiv_good.Draw("histsame")
csinkdiv.SetGridx(1)
csinkdiv.SetGridy(1)
csinkdiv.Update()
csinkdiv.SaveAs("csinkdiv.png")

cfracerr = rt.TCanvas("cfracerr","Fractional Error",1000,800)
hfracerr_bad.Draw("hist")
hfracerr_good.Draw("histsame")
cfracerr.SetGridx(1)
cfracerr.SetGridy(1)
cfracerr.Update()
cfracerr.SaveAs("cfracerr.png")

cobspe = rt.TCanvas("cobspe","Observed Total PE",1000,800)
hobspe_bad.Draw("hist")
hobspe_good.Draw("histsame")
cobspe.SetGridx(1)
cobspe.SetGridy(1)
cobspe.Update()
cobspe.SaveAs("cobspe.png")


input()

out.Write()


    
