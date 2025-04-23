import os,sys

sample="mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
fcomplete = open("%s_complete.txt"%(sample),'r')
lcomplete = fcomplete.readlines()

haddout = sample+"_reco_v2me06_gen2ntuple_photon_edep_fix_preview_v2.root"
command = "hadd -f %s "%(haddout)

fdict = {}
arr_v = []
for l in lcomplete:
    l = l.strip()
    arrayid = int(l.split("_")[-1].split(".")[0])
    #command += " "+l
    fdict[arrayid] = l
    arr_v.append(arrayid)
arr_v.sort()
for k in arr_v:
    command += " "+fdict[k]
    
print(command)
os.system(command)
