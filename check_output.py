import os,sys
import ROOT as rt

completed_arrayids = []
incomplete_arrayids = []
missing_arrayids = []

#sample="mcc9_v40_bnb_nu_overlay_500k_CV_run3b"
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/kpsana_list_mcc9_v40_bnb_nu_overlay_500k_CV_run3b.txt"
#narrayids = 249
#nfiles=20

#sample="mcc9_v28_wctagger_bnb5e19"
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v28_wctagger_bnb5e19.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnb5e19.txt"
#narrayids = 583
#nfiles=20

#sample="mcc9_v29e_dl_run3_G1_extbnb_dlana"
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt"
#narrayids=584
#nfiles=20

#sample="v3dev_reco_retune"
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v40a_dl_run1_bnb_intrinsic_nue_overlay_CV.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v40a_dl_run1_bnb_intrinsic_nue_overlay_CV_v3dev_reco_retune.txt"
#nfiles=5
#narrayids=1

sample="mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune"
nfiles=5
narrayids=98
bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v28_wctagger_nueintrinsics.txt"
recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune.txt"

#sample="mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune"
#nfiles=30
#narrayids=315
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v28_wctagger_bnboverlay.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.txt"

#sample="mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune"
#nfiles=100
#narrayids=235
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v29e_dl_run1_C1_extbnb.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune.txt"

#sample="mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune"
#nfiles=40
#narrayids=293
#bookkeeping="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/bookkeeping/fileinfo_mcc9_v28_wctagger_bnb5e19.txt"
#recolist="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune.txt"


# open and parse the book keeping file
fileidmap = {}
with open(bookkeeping,'r') as f:
    ll = f.readlines()
    for l in ll:
        l = l.strip()
        info = l.split()
        fileid = int(info[0])
        nentries = int(info[1])
        fileidmap[fileid] = nentries

# open the recolist in order to determine the expected num entries for each output ntuple file
ntuple_num_expected = {}
ntuple_num_inputfiles = {}
with open(recolist,'r') as f:
    ll = f.readlines()
    i = 0
    iout = 0
    nexpected = 0
    ninfiles = 0
    for l in ll:
        l = l.strip()
        fbase = os.path.basename(l)        
        info = fbase.split("_")
        fid = int(info[1][len("fileid"):])
        print(fbase,": ",fid)
        if fid in fileidmap:
            nexpected += fileidmap[fid]
            ninfiles += 1
        else:
            pass
        i+=1
        if i==nfiles:
            ntuple_num_expected[iout] = nexpected
            ntuple_num_inputfiles[iout] = ninfiles
            # reset
            i=0
            nexpected = 0
            ninfiles = 0
            iout += 1
        else:
            pass
        
fmissing_out = open("%s_test_incomplete.txt"%(sample),'w')
fcomplete_out = open("%s_test_complete.txt"%(sample),'w')

output_dir="./out_test/%s/"%(sample)
print("Scan output ntuple directory: ",output_dir)
ntotentries = 0
outputmap = {}
flist = os.listdir(output_dir)
for f in flist:
    f = f.strip()
    fpath = "./out_test/%s/"%(sample)+f

    if os.path.isdir(fpath):
        # if the file is a directory

        # the directories have the array id in its name
        arrayid = int(f.split("_")[-1])
        print("fpath-directory [arrayid=%d]: "%(arrayid),fpath)

        # create the sub-ntuple file name that combines the files in the directory
        subntuple = "./out_test/%s/subntuple_%s_%03d.root"%(sample,sample,arrayid)

        # if the subntuple does not exist, we run hadd to make it
        if not os.path.exists(subntuple):
            os.system(f"hadd -f {subntuple} {fpath}/*.root")
            #raise(" did not find subntuple file for this directory: ",subntuple)
    else:

        # check if the file is one of our subntuple files
        if len(f)<9 or f[:9]!="subntuple" or f[-5:]!=".root":
            print("not a subntuple file: ",f)
            continue

        #print(f.strip())
        # get the arrayid from the name of the file
        arrayid = int(f.strip().split("_")[-1].split(".")[0])

        # open the file to get the num of entries
        rfile = rt.TFile( fpath, "open" )
        eventtree = rfile.Get("EventTree")
        nentries  = eventtree.GetEntries()
        rfile.Close()
        print("file[arrayid=%d]: "%(arrayid),f," nentries=",nentries)

        # get the number of entries that are expected
        # from summing the block of input ntuple files
        if arrayid not in ntuple_num_expected:
            # if the arrayid is not in the compiled arrayid
            incomplete_arrayids.append(arrayid)        
            continue

        # get number of expected events from summing the ntuple files
        nexpected   = ntuple_num_expected[arrayid]
        # and also the number of files
        ninputfiles = ntuple_num_inputfiles[arrayid]
        
        if sample in ["mcc9_v40_bnb_nu_overlay_500k_CV_run3b"]:
            pottree = rfile.Get("potTree")    
            npotentries = 0
            for ii in range(pottree.GetEntries()):
                pottree.GetEntry(ii)
                if pottree.totGoodPOT>0:
                    npotentries += 1
        else:
            # force a pass here?
            npotentries=ninputfiles

        # check if the number of entries matches that in the summed ntuple files
        if nexpected!=nentries or npotentries!=ninputfiles:
            print("output[",arrayid,"] Incomplete output file. EventTree entries=",nentries,
                  " vs expected=",nexpected," pottree entries=",npotentries," vs. expected=",ninputfiles)
            print("  ",f)
            print(fpath,file=fmissing_out)
            incomplete_arrayids.append(arrayid)
        else:
            print("output[",arrayid,"] EventTree match with nentries=",nentries," and pottree match with nentries=",npotentries)
            ntotentries += nentries
            completed_arrayids.append(arrayid)
            outputmap[arrayid] = fpath
            print(fpath,file=fcomplete_out)

fmissing_out.close()

completed_arrayids.sort()
incomplete_arrayids.sort()
print("Completed. N=",len(completed_arrayids),": ",completed_arrayids)
print("Incomplete. N=",len(incomplete_arrayids),": ",incomplete_arrayids)
print("Num total entries completed=",ntotentries)

with open('output_goodlist_%s.txt'%(sample),'w') as f:
    for arrayid in completed_arrayids:
        print(outputmap[arrayid],file=f)
    

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


