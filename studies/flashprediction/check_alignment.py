import os, sys
import ROOT as rt

ntuple_goodlist="../../output_goodlist_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.txt"

ntuplefile = open(ntuple_goodlist,'r')
lines = ntuplefile.readlines()

output_dir="./output/"
output_list = os.listdir(output_dir)

ntuple_gen2="../../ntuple_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"
gen2file = rt.TFile(ntuple_gen2)
gen2tree = gen2file.Get("EventTree")
nentries = gen2tree.GetEntries()

flashpred_filename="flashprediction_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"
flashpred_file = rt.TFile(flashpred_filename)
flashpred_tree = flashpred_file.Get("FlashPredictionTree")
nentries_flashpred = flashpred_tree.GetEntries()

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
