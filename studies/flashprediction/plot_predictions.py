import os,sys
import ROOT as rt
rt.gStyle.SetOptStat(0)

rfile = rt.TFile(sys.argv[1])
ttree = rfile.Get("FlashPredictionTree")

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
            hsiren.SetBinContent( ipmt+1, ttree.siren_pe_per_pmt_all.at(ivtx).at(ipmt)*3.0 )


        hmax = hobs
        maxmax  = 0.0
        for h in [hobs, hsiren, hublm]:
            if h.GetMaximum()>maxmax:
                maxmax = h.GetMaximum()
                hmax = h
        hmax.Draw()
        hobs.Draw("same")
        hublm.Draw("same")
        hsiren.Draw("same")
        c.Update()
        print("enter for next vertex")
        input()

    



