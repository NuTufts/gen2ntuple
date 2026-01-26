import os,sys

INPUT_FILELISTS = {
    "run3_opendata_1e19":"/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/filelist_mcc9_v28_wctagger_run3_bnb1e19.txt",
}

KPSRECO_DIR = {
    "run3_opendata_1e19":"/cluster/tufts/wongjiradlabnu/nutufts/data/v2_me_06_03_prod/mcc9_v28_wctagger_run3_bnb1e19/larflowreco/ana/"
}


def get_filelist( sample ):
    """
    
    """

    fileid_dict = {}
    
    # Parse dlmerged list
    inputlist = INPUT_FILELISTS[sample]
    finputlist = open(inputlist,'r')
    
    fileid=0
    for linput in finputlist.readlines():
        dlmerged = linput.strip()
        fileid_dict[fileid] = {'dlmerged':dlmerged}
        fileid += 1

    # Parse reco files
    kpsreco_dir = KPSRECO_DIR[sample]
    precolist = os.popen(f"find {kpsreco_dir} -type f | grep root")
    lrecolist = precolist.readlines()
    ilreco = 0
    for lreco in lrecolist:
        lreco = lreco.strip()
        lrecobase = os.path.basename( lreco )

        recosplit = lrecobase.split("_")
        fileidsplit = recosplit[1]
        fileid = int(fileidsplit[len("fileid"):])
        #print(f"{fileid}: {lreco}")
        fileid_dict[fileid]['reco'] = lreco        
        ilreco += 1

    print(f"Number of dlreco/reco files: {len(fileid_dict)}")

    return fileid_dict
