import os,sys
import ROOT as rt
from math import sqrt

inputfile = "photonsel_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"
ftfile = rt.TFile(inputfile)
ttree  = ftfile.Get("PhotonSelectionTree")

# Make an output file: a place for histograms to be retrieved from the ROOT system
tout = rt.TFile("temp.root","recreate")

norm = True


sig_def = "out_is_target==1 && dwall_reco_nuvtx>20.0"
bg_def  = "out_num_event_vertices>1 && out_has_visible_primary_photon==1 && out_is_target==0"

# total PE fractional error
cfracerr = rt.TCanvas("cfracerr","Total PE Fractional Error",1000,800)
ttree.Draw("(1.0+totpefracerr)>>hfracerr_bg(150, 0.0,5.0)",bg_def,"hist")
ttree.Draw("(1.0+totpefracerr)>>hfracerr_sig(150,0.0,5.0)",sig_def,"histsame")


hfracerr_bg  = tout.Get("hfracerr_bg")
hfracerr_sig = tout.Get("hfracerr_sig")
hfracerr_sig.SetLineColor(rt.kRed)

binwidth_fracerr = (5.0)/150.0
n_hfracerr_bg  = hfracerr_bg.Integral()
n_hfracerr_sig = hfracerr_sig.Integral()

# approximating likelihood
fL_farwall_fracerr  = rt.TF1( "farwall_fracerr", "gaus(0) + [3]*exp(-[4]*x)", 0.0, 5.0 )
#farwall_frac_const = 190.0
farwall_frac_const = 190.0
farwall_frac_mean  = 1.0
farwall_frac_sigma = 0.7
#farwall_frac_expconst  = 148.0
farwall_frac_expconst  = 0.0
farwall_frac_lambda    = 1.0
farwall_frac_gaus_norm = farwall_frac_sigma*sqrt(2.0*3.14159)*farwall_frac_const
farwall_frac_exp_norm  = farwall_frac_expconst/farwall_frac_lambda

farwall_frac_totalnorm = farwall_frac_gaus_norm+farwall_frac_exp_norm
fL_farwall_fracerr.SetParameter(0, farwall_frac_const/farwall_frac_totalnorm)
fL_farwall_fracerr.SetParameter(1, farwall_frac_mean)
fL_farwall_fracerr.SetParameter(2, farwall_frac_sigma)
fL_farwall_fracerr.SetParameter(3, farwall_frac_expconst/farwall_frac_totalnorm)
fL_farwall_fracerr.SetParameter(4, farwall_frac_lambda)
fL_farwall_fracerr.Draw("same")

# llfracerr = rt.TF1("llfracerr","landau",0,10.0)
# llfracerr.SetParameter(0,8.00)
# llfracerr.SetParameter(1,0.30)
# llfracerr.SetParameter(2,0.18)
# llfracerr.Draw("same")

if norm:
    hfracerr_bg.Scale(  1.0/(n_hfracerr_bg*binwidth_fracerr) )
    hfracerr_sig.Scale( 1.0/(n_hfracerr_sig*binwidth_fracerr) )
    #hfracerr_bg.Scale(  1.0/(n_hfracerr_bg) )
    #hfracerr_sig.Scale( 1.0/(n_hfracerr_sig) )


cfracerr.Update()

# =====================================================================
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
fL_farwall_sinkdiv  = rt.TF1( "farwall_sinkdiv", "gaus(0) + [3]*exp(-[4]*x)", 0.0, 40.0 )
#farwall_sink_const = 96.0
farwall_sink_const = 0.0
farwall_sink_mean  = 1.4
farwall_sink_sigma = 0.7
farwall_sink_expconst  = 245.0
#farwall_sink_lambda    = 0.17
farwall_sink_lambda    = 0.10
farwall_sink_gaus_norm = farwall_sink_sigma*sqrt(2.0*3.14159)*farwall_sink_const
farwall_sink_exp_norm  = farwall_sink_expconst/farwall_sink_lambda
farwall_sink_totalnorm = farwall_sink_gaus_norm+farwall_sink_exp_norm
fL_farwall_sinkdiv.SetParameter(0, farwall_sink_const/farwall_sink_totalnorm)
fL_farwall_sinkdiv.SetParameter(1, farwall_sink_mean)
fL_farwall_sinkdiv.SetParameter(2, farwall_sink_sigma)
fL_farwall_sinkdiv.SetParameter(3, farwall_sink_expconst/farwall_sink_totalnorm)
fL_farwall_sinkdiv.SetParameter(4, farwall_sink_lambda)


# llsink = rt.TF1("llsink","landau",0,40.0)
# llsink.SetParameter(0,1.0)
# llsink.SetParameter(1,1.45)
# llsink.SetParameter(2,1.45)
# llsink.Draw("same")
fL_farwall_sinkdiv.Draw("same")

if norm:
    hsink_bg.Scale( 1.0/(n_hsink_bg*binwidth_sink) )
    hsink_sig.Scale( 1.0/(n_hsink_sig*binwidth_sink) )


csink.Update()

print("[enter] to exit")
input()

tout.Write()


