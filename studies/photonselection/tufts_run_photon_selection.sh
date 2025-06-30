#!/bin/bash

JOBSTARTDATE=$(date)

RECO_FILELIST=$1
DLMERGED_FILELIST=$2
SAMPLENAME=$3
NFILES=$4
MCFLAG=$5
ubdlDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/ubdl/
outDir=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/photonselection/output/${SAMPLENAME}/
PHOTON_SELECTION_DIR=/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/studies/photonselection/

localDir=`printf /tmp/photonselection_jobarrayid%05d ${SLURM_ARRAY_TASK_ID}`
mkdir -p ${localDir}
mkdir -p ${outDir}

source ${ubdlDir}/setenv_py3_container.sh
source ${ubdlDir}/configure_container.sh
export PYTHONPATH=${PYTHONPATH}:${scriptDir}
export PATH=${PHOTON_SELECTION_DIR}/build/installed/bin:${PATH}

cd ${localDir}

maxFileCount=`wc -l < ${RECO_FILELIST}`
let firstfile="${SLURM_ARRAY_TASK_ID}*${NFILES}+1"
let lastfile="${firstfile}+$NFILES-1"
echo "files to run between first=${firstfile} to last=${lastfile} maxfilecount=${maxFileCount}"
for n in $(seq $firstfile $lastfile); do
  if (($n > $maxFileCount)); then
    break
  fi

  fileidx=`sed -n ${n}p ${RECO_FILELIST} | awk '{ print $1 }'`
  recofile=`sed -n ${n}p ${RECO_FILELIST} | awk '{ print $2 }'`
  let mergedfile_lineno="${fileidx}+1"
  mergedfile=`sed -n ${mergedfile_lineno}p ${DLMERGED_FILELIST}`
  recobase=`basename ${recofile}`
  outfile=`echo ${recobase} | sed 's|larflowreco|photonselection|g'`
  cp ${recofile} ${recobase}

  fileOutDir=`printf ${outDir}/ntuplefile%05d ${SLURM_ARRAY_TASK_ID}`
  mkdir -p ${fileOutDir}
  
  echo "mergedfile: ${mergedfile}"
  echo "recobase: ${recobase}"
  echo "outfile: ${outfile}"
  CMD="make_photon_selection_study_tree ${mergedfile} ${recobase} ${outfile} --tickbackward"
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

