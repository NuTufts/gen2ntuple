import os,sys

samplename="mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune"
#samplename="mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune"
#samplename="mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune"
#samplename="mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune"

output_dir=f"./output/{samplename}/"
output_list = os.listdir(output_dir)
for outdir in output_list:
    foutdir=f"{output_dir}/{outdir}"    
    if not os.path.isdir(foutdir):
        print("not a dir: ",foutdir)
        continue
    
    arrayid = int(outdir[len("ntuplefile"):])
    #if arrayid not in [99]:
    #    continue

    #print(outdir)
    ntupleout = f"{outdir}.root"
    filedict = {}
    fileid_v = []
    print(f"make {ntupleout} from files in {outdir} directory")
    subntuple_files = os.listdir( foutdir )
    for f in subntuple_files:
        fileid = int(f.split("_")[1][len("fileid"):])
        print("[",fileid,"]: ",f)
        fileid_v.append(fileid)
        filedict[fileid] = f

    fileid_v.sort()
    hadd_flist=""
    for fileid in fileid_v:
        hadd_flist += " %s/%s"%( foutdir, filedict[fileid] )    
    os.system(f"hadd -f {output_dir}/{ntupleout} {hadd_flist}")

