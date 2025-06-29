#!/bin/bash

#SBATCH --job-name=phselect
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0
##SBATCH --partition=wongjiradlab
#SBATCH --partition=batch
##SBATCH --exclude=i2cmp006,s1cmp001,s1cmp002,s1cmp003,p1cmp005,p1cmp041,c1cmp003,c1cmp004i
##SBATCH --exclude=p1cmp075
#SBATCH --error=err/griderr_photon_selection.sub00.%A.%a.node%N.err
#SBATCH --output=log/stdout_photon_selection.sub00.%A.%a.node%N.log

# Container to run in -- needs to be same as the one used to build UBDL
CONTAINER=/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/u20.04_cu111_cudnn8_torch1.9.0_minkowski_npm.sif

# Location where this script is located
REPO_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/photonselection/

# Location of repo copy that was used to run the lantern reco
# -- provides the filelists to run on
LMRECO_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/

# BELOW ARE CONFIGURATIONS FOR DIFFERENT SAMPLES
# We follow parameters used by the gen2ntuple maker in order to help align the flash prediction
# tree with the ntuple tree.

# mcc9_v28_wctagger_bnboverlay: run 1 BNB nu overlay
# n good reco files: 9462
# njobs=316
NFILES=30
SAMPLENAME=mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune
RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnboverlay_v3dev_reco_retune.txt
TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v28_wctagger_bnboverlay.txt
MCFLAG="-mc"

# mcc9_v28_wctagger_nueintrinsics: run 1 BNB nue intrinsics overlay
# 490 good reco files
#NFILES=5
# njobs=98
#SAMPLENAME=mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune
#RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_nueintrinsics_v3dev_reco_retune.txt
#TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v28_wctagger_nueintrinsics.txt
#MCFLAG="-mc"

# mcc9_v29e_dl_run1_C1_extbnb : run 1 EXTBNB
#NFILES=100
## total files in good list: 23431
## number of jobs is 234
#SAMPLENAME=mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune
#RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v29e_dl_run1_C1_extbnb_v3dev_reco_retune.txt
#TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v29e_dl_run1_C1_extbnb.txt
#MCFLAG=""

# mcc9_v28_wctagger_bnb5e19 : run 1 open data
# number of files: 11681
#NFILES=40
## number of jobs=293
#SAMPLENAME=mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune
#RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnb5e19_v3dev_reco_retune.txt
#TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v28_wctagger_bnb5e19.txt
#MCFLAG=""

#RECOFILELIST=/cluster/tufts/wongjiradlabnu/mrosen25/filelists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_filelist.txt
#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/mrosen25/filelists/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_filelist.txt

#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_matched_mergedlist.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/mcc9_v29e_dl_run3b_bnb_nu_overlay_nocrtremerge_matched_recolist.txt

#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/dlmerged_list_mcc9_v40_NC_Pi0_run3b.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/kpsana_list_mcc9_v40_NC_Pi0_run3b.txt

#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/dlmerged_list_mcc9_v40_NC_Pi0_run3b.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/kpsana_list_mcc9_v40_NC_Pi0_run3b.txt

# 223 files
#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/input_dlrecolist_epem_darknu_benchD.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/kpsana_list_epem_DarkNu_BenchmarkD.txt

# a lot of files
#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/filelist_mcc9_v40a_dl_run3b_bnb_nu_overlay_500k_CV.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/kpsana_list_mcc9_v40_bnb_nu_overlay_500k_CV_run3b.txt
#SAMPLENAME=mcc9_v40_bnb_nu_overlay_500k_CV_run3b

# MCC9 Run 3 G1 EXTBNB: for cosmic estimates. ED1CNT_wcut is 39195178 for full sample. Total number of entries is 234610.
# good list number of files: 10970 --> number of jobs 10970/20=548
# requires -ana flag
#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/filelist_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt
#RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt
#SAMPLENAME=mcc9_v29e_dl_run3_G1_extbnb_dlana

# MCC9 v40
#TRUTHFILELIST=${LMRECO_DIR}/filelists/filelist_mcc9_v40a_dl_run1_bnb_intrinsic_nue_overlay_CV.txt
#RECOFILELIST=${LMRECO_DIR}/goodoutput_lists/goodoutput_list_mcc9_v40a_dl_run1_bnb_intrinsic_nue_overlay_CV_v3dev_reco_retune.txt
#SAMPLENAME=v3dev_reco_retune

# MCC9 v28 BNB 5e19
# Num reco files: 11671
# Njobs with 20 files/job: 584
#TRUTHFILELIST="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/filelist_mcc9_v28_wctagger_bnb5e19.txt"
#RECOFILELIST="/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v28_wctagger_bnb5e19.txt"
#SAMPLENAME="mcc9_v28_wctagger_bnb5e19"


OUTTAG=v3dev_reco_retune
BINDING=/cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab

#module load singularity/3.5.3
module load apptainer/1.2.4-suid

apptainer exec --bind ${BINDING} ${CONTAINER} bash -c "cd ${REPO_DIR} && ./tufts_run_photon_selection.sh ${RECOFILELIST} ${TRUTHFILELIST} ${SAMPLENAME} ${NFILES} ${MCFLAG}"
