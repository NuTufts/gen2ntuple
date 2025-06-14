#!/bin/bash

JOBSTARTDATE=$(date)

RECO_FILELIST=$1
DLMERGED_FILELIST=$2
NFILES=$3
MCFLAG=$4
ubdlDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/ubdl/
outDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/output/
FLASH_PREDICTION_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/flashprediction/

localDir=`printf /tmp/calc_flash_prediction_jobarrayid%05d ${SLURM_ARRAY_TASK_ID}`
mkdir -p ${localDir}
mkdir -p ${outDir}

source ${ubdlDir}/setenv_py3_container.sh
source ${ubdlDir}/configure_container.sh
export PYTHONPATH=${PYTHONPATH}:${scriptDir}
export PATH=${FLASH_PREDICTION_DIR}/build/installed/bin:${PATH}

cd ${localDir}

maxFileCount=`wc -l < $kpsRecoFiles`
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
  CMD="calculate_flash_predictions --dlmerged ${mergedfile} --reco ${recofile} --output ${outfile} -tb ${MCFLAG}"
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

