import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

rfile = rt.TFile(sys.argv[1])
ttree = rfile.Get("FlashPredictionTree")

scale_factor=1.0

nentries = ttree.GetEntries()

c = rt.TCanvas("c","c",800,600)
c.Draw()

hobs   = rt.TH1D("hobs","",32,0,32)
hsiren = rt.TH1D("hsiren","",32,0,32)
hublm  = rt.TH1D("hublm","",32,0,32)
hobs.SetLineColor(rt.kBlack)
hsiren.SetLineColor(rt.kRed)
hublm.SetLineColor(rt.kBlue)

for ientry in range(nentries):
    ttree.GetEntry(ientry)

    hobs.Reset()    
    for ipmt in range(32):
        hobs.SetBinContent( ipmt+1, ttree.obs_pe_per_pmt.at(ipmt) )
        
    for ivtx in range(ttree.n_vertices):        

        print('Event: ')
        print('  run: ',ttree.run)
        print('  subrun: ',ttree.subrun)
        print('  event: ',ttree.event)
        print('  vertex: ',ivtx)

        hsiren.Reset()
        hublm.Reset()
        for ipmt in range(32):
            hublm.SetBinContent(  ipmt+1, ttree.ubpred_pe_per_pmt_all.at(ivtx).at(ipmt) )
            hsiren.SetBinContent( ipmt+1, ttree.siren_pe_per_pmt_all.at(ivtx).at(ipmt)*scale_factor )

        sinkdiv_siren     = ttree.siren_sinkhorn_div_all.at(ivtx).at(0)
        sinkdiv_ublm      = ttree.ub_sinkhorn_div_all.at(ivtx).at(0)
        unb_sinkdiv_siren = ttree.siren_unbalanced_sinkhorn_div_all.at(ivtx).at(0)
        unb_sinkdiv_ublm  = ttree.ub_unbalanced_sinkhorn_div_all.at(ivtx).at(0)
        ttext_siren     = rt.TText(0.2,0.80,f"Siren: balanced sinkdiv {sinkdiv_siren:0.2e}")
        ttext_ublm      = rt.TText(0.2,0.75,f"UB LM: balanced sinkdiv {sinkdiv_ublm:0.2e}")    
        ttext_unb_siren = rt.TText(0.2,0.70,f"Siren: balanced sinkdiv {unb_sinkdiv_siren:0.2e}")
        ttext_unb_ublm  = rt.TText(0.2,0.65,f"UB LM: balanced sinkdiv {unb_sinkdiv_ublm:0.2e}")    
        ttext_siren.SetNDC(True)
        ttext_ublm.SetNDC(True)   
        ttext_unb_siren.SetNDC(True)
        ttext_unb_ublm.SetNDC(True)  


        hmax = hobs
        maxmax  = 0.0
        for h in [hobs, hsiren, hublm]:
            h.SetLineWidth(2)
            if h.GetMaximum()>maxmax:
                maxmax = h.GetMaximum()
                hmax = h
        hmax.Draw()
        hublm.Draw("same")
        hsiren.Draw("same")
        hobs.Draw("same")
        ttext_siren.Draw()
        ttext_ublm.Draw()
        ttext_unb_siren.Draw()
        ttext_unb_ublm.Draw()
        c.Update()
        print("enter for next vertex")
        input()

    



