import os,sys

sample="mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
fmissing = open("%s_incomplete.txt"%(sample),'r')
lmissing = fmissing.readlines()

for l in lmissing:
    l = l.strip()
    command = "rm "+l
    print(command)
    os.system(command)
