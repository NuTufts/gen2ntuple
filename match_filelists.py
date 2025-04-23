import os,sys

recolist="/cluster/tufts/wongjiradlabnu/mrosen25/filelists/larflowreco_v2_me_05/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_filelist.txt"
mergedlist="/cluster/tufts/wongjiradlabnu/mrosen25/filelists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_filelist.txt"
samplename="mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
# parse merged list, find

matched = {}
hashlist = []

# parse merged filenames
fmerged = open(mergedlist,'r')
lmerged = fmerged.readlines()

i = 0
for l in lmerged:
    l = l.strip()
    fname = os.path.basename(l)
    hashname = fname.split("_")[-1].split(".")[0]
    #print(hashname)
    matched[hashname] = {"mergedfile":l}
    hashlist.append(hashname)
    i += 1
    if i%1000==0:
        print("processing file: ",i)
    #if i>=10:
    #    break
print("number of merged list: ",i)

freco = open(recolist,'r')
lreco = freco.readlines()
for l in lreco:
    l = l.strip()
    fname = os.path.basename(l)
    hashname = fname.split("_")[2]
    print(hashname)
    if hashname in matched:
        matched[hashname]["recofile"] = l

fout_merged = open("%s_matched_mergedlist.txt"%(samplename),'w')
fout_reco   = open("%s_matched_recolist.txt"%(samplename),'w')

hashlist.sort()
nmatched = 0
for k in hashlist:
    filedict = matched[k]
    if "recofile" in filedict:
        nmatched += 1
        print(filedict["mergedfile"],file=fout_merged)
        print(filedict["recofile"],file=fout_reco)
print("number matched: ",nmatched," out of ",i)



