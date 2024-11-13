#!/bin/bash

#SBATCH --job-name=nuNTuple
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000
#SBATCH --array=0-99
##SBATCH --partition=wongjiradlab
#SBATCH --partition=batch
#SBATCH --exclude=i2cmp006,s1cmp001,s1cmp002,s1cmp003,p1cmp005,p1cmp041,c1cmp003,c1cmp004

NFILES=100

CONTAINER=/cluster/tufts/wongjiradlabnu/larbys/larbys-container/singularity_minkowskiengine_u20.04.cu111.torch1.9.0_comput8.sif

VALSCRIPT=/path/to/tufts_run_ntuple_maker.sh
RECOFILELIST=/path/to/larflowreco_filelist.txt
TRUTHFILELIST=/path/to/merged_dlreco_filelist.txt
WEIGHTFILE=weights_file_name.pkl
OUTTAG=reco_v2me06_gen2ntuple_v5
CNNMODEL=run3bOverlays_quadTask_plAll_2inChan_5ClassHard_minHit10_b64_oneCycleLR_v2me05_noPCTrainOrVal/ResNet34_recoProng_5class_epoch20.pt

module load singularity/3.5.3

singularity exec --bind /cluster/tufts/wongjiradlabnu:/cluster/tufts/wongjiradlabnu,/cluster/tufts/wongjiradlab:/cluster/tufts/wongjiradlab ${CONTAINER} bash -c "source $VALSCRIPT make_dlgen2_flat_ntuples.py $RECOFILELIST $TRUTHFILELIST $WEIGHTFILE $CNNMODEL $OUTTAG $NFILES -mc"

