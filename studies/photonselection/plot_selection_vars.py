import os,sys
import ROOT as rt

inputfile = "photonsel_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"
ftfile = rt.TFile(inputfile)
ttree  = ftfile.Get("PhotonSelectionTree")

# Make an output file: a place for histograms to be retrieved from the ROOT system
tout = rt.TFile("temp.root","recreate")

norm = True


sig_def = "dwall_true_edep>5.0 && dist2photonedep<5.0"
bg_def  = "dwall_true_edep<0.0 || dist2photonedep>10.0"

# total PE fractional error
cfracerr = rt.TCanvas("cfracerr","Total PE Fractional Error",1000,800)
ttree.Draw("(1.0+totpefracerr)>>hfracerr_bg(150,-1.5,4.0)",bg_def,"hist")
ttree.Draw("(1.0+totpefracerr)>>hfracerr_sig(150,-1.5,4.0)",sig_def,"histsame")


hfracerr_bg  = tout.Get("hfracerr_bg")
hfracerr_sig = tout.Get("hfracerr_sig")
hfracerr_sig.SetLineColor(rt.kRed)

binwidth_fracerr = (4.0+1.5)/150.0
n_hfracerr_bg  = hfracerr_bg.Integral()
n_hfracerr_sig = hfracerr_sig.Integral()

# approximating likelihood
llfracerr = rt.TF1("llfracerr","landau",0,10.0)
llfracerr.SetParameter(0,8.00)
llfracerr.SetParameter(1,0.30)
llfracerr.SetParameter(2,0.18)
llfracerr.Draw("same")

if norm:
    hfracerr_bg.Scale(  1.0/(n_hfracerr_bg*binwidth_fracerr) )
    hfracerr_sig.Scale( 1.0/(n_hfracerr_sig*binwidth_fracerr) )


cfracerr.Update()


# sinkhorn divergence
csink = rt.TCanvas("csink","Sinkhorn Divergence",1000,800)
ttree.Draw("sinkhorndiv>>hsink_sig(100,0,40)",sig_def,"hist")
ttree.Draw("sinkhorndiv>>hsink_bg(100,0,40)",bg_def,"histsame")
hsink_bg  = tout.Get("hsink_bg")
hsink_sig = tout.Get("hsink_sig")
hsink_sig.SetLineColor(rt.kRed)

binwidth_sink = 40.0/100.0
n_hsink_bg  = hsink_bg.Integral()
n_hsink_sig = hsink_sig.Integral()

# approximating likelihood
llsink = rt.TF1("llsink","landau",0,40.0)
llsink.SetParameter(0,1.0)
llsink.SetParameter(1,1.45)
llsink.SetParameter(2,1.45)
llsink.Draw("same")

if norm:
    hsink_bg.Scale( 1.0/(n_hsink_bg*binwidth_sink) )
    hsink_sig.Scale( 1.0/(n_hsink_sig*binwidth_sink) )


csink.Update()

print("[enter] to exit")
input()

tout.Write()


