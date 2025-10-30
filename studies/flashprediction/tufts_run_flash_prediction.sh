#!/bin/bash

JOBSTARTDATE=$(date)

RECO_FILELIST=$1
DLMERGED_FILELIST=$2
SAMPLENAME=$3
NFILES=$4
MCFLAG=$5

ubdlDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/ubdl/
outDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/output/${SAMPLENAME}/
FLASH_PREDICTION_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/
FLASHMATCH_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/flashmatchdata_petastorm/
#SIREN_MODEL_FILE=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/siren_model_extbnb_good_sun_115k.pt
#SIREN_MODEL_FILE=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/flashmlp_model_extbnb_deft_universe_iteration_075k.pt
SIREN_MODEL_FILE=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/siren_model_extbnb_floral_shape_059k.pt

localDir=`printf /tmp/calc_flash_prediction_jobarrayid%05d ${SLURM_ARRAY_TASK_ID}`
mkdir -p ${localDir}
mkdir -p ${outDir}

source ${ubdlDir}/setenv_py3_container.sh
source ${ubdlDir}/configure_container.sh
export PYTHONPATH=${PYTHONPATH}:${scriptDir}
export PATH=${FLASH_PREDICTION_DIR}/build/installed/bin:${PATH}
cd ${FLASHMATCH_DIR}
source setenv_flashmatchdata.sh

echo "LIBTORCH DIR: ${LIBTORCH_DIR}"
echo "LD_LIBRARY_PATH: ${LD_LIBRARY_PATH}"

cd ${localDir}

maxFileCount=`wc -l < ${RECO_FILELIST}`
let firstfile="${SLURM_ARRAY_TASK_ID}*${NFILES}+1"
let lastfile="${firstfile}+$NFILES-1"
echo "filest to run between first=${firstfile} to last=${lastfile}"
for n in $(seq $firstfile $lastfile); do
  if (($n > $maxFileCount)); then
    break
  fi

  fileidx=`sed -n ${n}p ${RECO_FILELIST} | awk '{ print $1 }'`
  recofile=`sed -n ${n}p ${RECO_FILELIST} | awk '{ print $2 }'`
  let mergedfile_lineno="${fileidx}+1"
  mergedfile=`sed -n ${mergedfile_lineno}p ${DLMERGED_FILELIST}`
  recobase=`basename ${recofile}`
  outfile=`echo ${recobase} | sed 's|larflowreco|flashprediction|g'`

  fileOutDir=`printf ${outDir}/ntuplefile%05d ${SLURM_ARRAY_TASK_ID}`
  mkdir -p ${fileOutDir}
  
  echo "mergedfile: ${mergedfile}"
  echo "recobase: ${recobase}"
  echo "outfile: ${outfile}"
  CMD="calculate_flash_predictions --dlmerged ${mergedfile} --reco ${recofile} --output ${outfile} -v -tb --siren-model-file ${SIREN_MODEL_FILE} --dataset-type ${MCFLAG}"
  echo ${CMD}
  ${CMD}
  
  echo "copy outfile to ${fileOutDir}"
  cp ${outfile} ${fileOutDir}/${outfile}
done

# clean-up
cd /tmp
rm -r ${localDir}

JOBENDDATE=$(date)

echo "Job began at $JOBSTARTDATE"
echo "Job ended at $JOBENDDATE"

