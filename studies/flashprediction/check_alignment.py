import os, sys
import ROOT as rt

#samplename="mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune"
samplename="mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune"
#samplename="mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune"

ntuple_gen2=f"../../ntuple_{samplename}.root"
gen2file = rt.TFile(ntuple_gen2)
gen2tree = gen2file.Get("EventTree")
nentries = gen2tree.GetEntries()

flashpred_filename=f"flashprediction_{samplename}.root"
flashpred_file = rt.TFile(flashpred_filename)
flashpred_tree = flashpred_file.Get("FlashPredictionTree")
nentries_flashpred = flashpred_tree.GetEntries()

print("ntuple file: ",ntuple_gen2)
print("flash prediction: ",flashpred_filename)
print("Gen2 ntuple nentries: ",nentries)
print("FlashPred ntuple nentries: ",nentries_flashpred)

# check where breakdown is, if there is a breakdown
pass_check = True
for ientry in range( nentries ):
    if ientry>0 and ientry%10000==0:
        print("checking entry ",ientry)
        
    gen2tree.GetEntry(ientry)
    flashpred_tree.GetEntry(ientry)
    gen2_rse = (gen2tree.run,gen2tree.subrun,gen2tree.event)
    flashpred_rse = (flashpred_tree.run,flashpred_tree.subrun,flashpred_tree.event)
    if gen2_rse!=flashpred_rse:
        print("Mis-alignment occurs at gen2-entry ",ientry)
        print("gen2_rse: ",gen2_rse)
        print("flash prediction rse: ",flashpred_rse)
        pass_check = False
        break

if pass_check:
    print("gen2 ntuple and flash prediction tree are aligned")
else:
    print("gen2 ntuple and flash prediction tree are mis-aligned")
