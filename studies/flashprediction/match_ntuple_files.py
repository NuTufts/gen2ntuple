import os, sys
import ROOT as rt

#samplename="mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune"
#samplename="mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune"
#samplename="mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune"
samplename="mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune"

#ntuple_goodlist="../../output_goodlist_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.txt"
#ntuple_goodlist=f"../../output_goodlist_{samplename}.txt"
ntuple_goodlist=f"../../{samplename}_test_complete.txt"

ntuplefile = open(ntuple_goodlist,'r')
lines = ntuplefile.readlines()

output_dir=f"./output/{samplename}/"
output_list = os.listdir(output_dir)

#ntuple_gen2="../../ntuple_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.root"
#gen2file = rt.TFile(ntuple_gen2)
#gen2tree = gen2file.Get("EventTree")
#nentries = gen2tree.GetEntries()

matching_files = []
for ll in lines:
    ll = ll.strip()

    # get the fileid number from the list
    fileid = int(ll.split("_")[-1].split(".")[0])
    print(ll,": ",fileid)

    gen2_tfile = rt.TFile("./../../"+ll)
    gen2_ttree = gen2_tfile.Get("EventTree")
    nentries_gen2 = gen2_ttree.GetEntries()
    
    ntuple_subfile="./output/%s/ntuplefile%05d.root"%(samplename,fileid)
    if not os.path.exists(ntuple_subfile):
        continue
    fpred_tfile = rt.TFile( ntuple_subfile )
    fpred_ttree = fpred_tfile.Get("FlashPredictionTree")
    nentries_fpred = fpred_ttree.GetEntries()

    if nentries_gen2!=nentries_fpred:
        print("number of event mismatch for ",ll)
        print("  gen2: ",nentries_gen2)
        print("  flash-pred: ",nentries_fpred)
        continue
    else:

        pass_check = True
        for ientry in range( nentries_gen2 ):
            if ientry>0 and ientry%10000==0:
                print("checking entry ",ientry)
        
            gen2_ttree.GetEntry(ientry)
            fpred_ttree.GetEntry(ientry)
            gen2_rse  = (gen2_ttree.run,gen2_ttree.subrun,gen2_ttree.event)
            fpred_rse = (fpred_ttree.run,fpred_ttree.subrun,fpred_ttree.event)
            if gen2_rse!=fpred_rse:
                print("Mis-alignment occurs at gen2-entry ",ientry)
                print("gen2_rse: ",gen2_rse)
                print("flash prediction rse: ",fpred_rse)
                pass_check = False
                break
        if not pass_check:
            continue
        
    
    matching_files.append( ntuple_subfile )

flist_matching = open(f"matching_flashpred_{samplename}.txt",'w')
for f in matching_files:
    print(f,file=flist_matching)
flist_matching.close()

os.system(f"hadd -f flashprediction_{samplename}.root @matching_flashpred_{samplename}.txt")
