import os,sys
import ROOT as rt

narrayids = 155
completed_arrayids = []
incomplete_arrayids = []
missing_arrayids = []

#sample="mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge"
sample="mcc9_v40_bnb_nu_overlay_500k_CV_run3b"

fmissing_out = open("%s_incomplete.txt"%(sample),'w')
fcomplete_out = open("%s_complete.txt"%(sample),'w')

flist = os.listdir("./out/")
for f in flist:
    #print(f.strip())
    fpath = "./out/"+f.strip()
    rfile = rt.TFile( fpath, "open" )
    pottree = rfile.Get("potTree")
    nentries = pottree.GetEntries()
    arrayid = int(f.strip().split("_")[-1].split(".")[0])
    if nentries!=100 and arrayid!=narrayids:
        print("Incomplete output file. expected=100, number of entries in the pottree is ",nentries)
        print(fpath,file=fmissing_out)
        incomplete_arrayids.append(arrayid)
    else:
        completed_arrayids.append(arrayid)
        print(fpath,file=fcomplete_out)
    rfile.Close()
fmissing_out.close()

completed_arrayids.sort()
incomplete_arrayids.sort()
print("Completed. N=",len(completed_arrayids),": ",completed_arrayids)
print("Incomplete. N=",len(incomplete_arrayids),": ",incomplete_arrayids)


for i in range(narrayids+1):
    if i in completed_arrayids or i in incomplete_arrayids:
        continue
    missing_arrayids.append(i)
arraystr = ""
lastid = -10
seqstart = -10
seqend = -10
for i in missing_arrayids:
    if lastid<0:
        # first in list
        lastid = i
        continue
    if seqstart>=0:
        # existing sequence defined
        if (i-1)==lastid:
            # continue sequence
            seqend = i
            lastid = i
        else:
            # existing seq broken, write sequence, start new possible seq
            if seqstart==seqend:
                arraystr += "%d,"%(seqstart)
            else:
                arraystr += "%d-%d,"%(seqstart,seqend)
            seqstart = -10
            lastid = i
    else:
        # existing sequence not defined yet
        if (i-1)==lastid:
            # start of sequence
            seqstart = lastid
            seqend = i
            lastid = i
        else:
            # individual entry, start seq
            #arraystr += "%d,"%(i)
            seqstart = i
            seqend = i
            lastid = i
    print("[",i,"]: seqstart=",seqstart," seqend=",seqend," lastid=",lastid)            
            
#arraystr = arraystr
print("missing array string: ",arraystr)


