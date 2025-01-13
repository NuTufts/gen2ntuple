#!/bin/bash

#SBATCH --job-name=nuNTuple
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0
##SBATCH --partition=wongjiradlab
#SBATCH --partition=batch
##SBATCH --exclude=i2cmp006,s1cmp001,s1cmp002,s1cmp003,p1cmp005,p1cmp041,c1cmp003,c1cmp004i
#SBATCH --exclude=p1cmp075
#SBATCH --error=err/griderr_ntuple_mcc9_v29e_dl_run3_G1_extbnb_dlana_500k_CV_sub00.%A.%a.node%N.err
#SBATCH --output=log/stdout_ntuple_mcc9_v29e_dl_run3_G1_extbnb_dlana_500k_CV_sub00.%A.%a.node%N.log

NFILES=20

CONTAINER=/cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowski_u20.04.cu111.torch1.9.0_jupyter_xgboost.sif

VALSCRIPT=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/tufts_run_ntuple_maker.sh

#RECOFILELIST=/cluster/tufts/wongjiradlabnu/mrosen25/filelists/larflowreco_v2_me_06/mcc9_v28_wctagger_bnboverlay_filelist.txt
#TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/mrosen25/filelists/mcc9_v28_wctagger_bnboverlay_filelist.txt
#WEIGHTFILE=weights_forCV_v48_Sep24_bnb_nu_run1.pkl

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
# good list number of files: 10970
TRUTHFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/filelists/filelist_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt
RECOFILELIST=/cluster/tufts/wongjiradlabnu/twongj01/gen2/dlgen2prod/larmatch_and_reco_scripts/goodoutput_lists/goodoutput_list_mcc9_v29e_dl_run3_G1_extbnb_dlana.txt
SAMPLENAME=mcc9_v29e_dl_run3_G1_extbnb_dlana


OUTTAG=reco_v3lmshowerkp_gen2ntuple_evisvertex

CNNMODEL=run3bOverlays_quadTask_plAll_2inChan_5ClassHard_minHit10_b64_oneCycleLR_v2me05_noPCTrainOrVal/ResNet34_recoProng_5class_epoch20.pt

#WEIGHTFILE=weights_forCV_v48_Sep24_bnb_nu_run1.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_bnb_nu_run2.pkl
WEIGHTFILE=weights_forCV_v48_Sep24_bnb_nu_run3.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_dirt_nu_run1.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_dirt_nu_run3.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_intrinsic_nue_run1.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_intrinsic_nue_run2.pkl
#WEIGHTFILE=weights_forCV_v48_Sep24_intrinsic_nue_run3.pkl

module load singularity/3.5.3

#singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${CONTAINER} bash -c "source $VALSCRIPT make_dlgen2_flat_ntuples.py $RECOFILELIST $TRUTHFILELIST $WEIGHTFILE $CNNMODEL $OUTTAG $NFILES ${SAMPLENAME} -mc"
singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${CONTAINER} bash -c "source $VALSCRIPT make_dlgen2_flat_ntuples.py $RECOFILELIST $TRUTHFILELIST $WEIGHTFILE $CNNMODEL $OUTTAG $NFILES ${SAMPLENAME} -ana"
