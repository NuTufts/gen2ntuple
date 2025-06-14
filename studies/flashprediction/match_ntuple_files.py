import os, sys

ntuple_goodlist="../../output_goodlist_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.txt"

ntuplefile = open(ntuple_goodlist,'r')
lines = ntuplefile.readlines()

output_dir="./output/"
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

    ntuple_subfile="./output/ntuplefile%05d.root"%(fileid)
    matching_files.append( ntuple_subfile )

flist_matching = open("matching_flashpred.txt",'w')
for f in matching_files:
    print(f,file=flist_matching)
flist_matching.close()
