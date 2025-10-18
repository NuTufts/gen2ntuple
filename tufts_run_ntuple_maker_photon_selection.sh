#!/bin/bash

JOBSTARTDATE=$(date)

scriptDir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple"
weightDir="/cluster/tufts/wongjiradlabnu/mrosen25/gen2ntuple/event_weighting/"

kpsRecoFiles=$1 # list of reco files
mdlRecoFiles=$2 # list of marged dlreco files

weightFile="${weightDir}/$3"

modelPath="/cluster/tufts/wongjiradlabnu/nutufts/larpid_weights/Scripted_LArPID_default_network_weights.pt"
outTag=$4
nfiles=$5
sampleName=$6

echo "Reco File List: ${kpsRecoFiles}"
echo "DLMerged File List: ${mdlRecoFiles}"
echo "weightFile: ${weightFile}"
echo "out tag: ${outTag}"
echo "nfiles: ${nfiles}"
echo "sampleName: ${sampleName}"

ubdlDir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/photon_analysis/ubdl"
gen2ntuple_dir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/"

outDir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/out_test/${sampleName}_${outTag}/"
logDir="/cluster/tufts/wongjiradlabnu/twongj01/gen2/gen2ntuple/log/${sampleName}_${outTag}/"

mkdir -p ${outDir}
mkdir -p ${logDir}

echo "running array ID $SLURM_ARRAY_TASK_ID (sample output tag: $outTag) on node $SLURMD_NODENAME"

cd ${ubdlDir}
source ${ubdlDir}/setenv_py3_container.sh
source ${ubdlDir}/configure_container.sh
cd ${gen2ntuple_dir}/gen2ntuple
source set_gen2ntuple_env.sh

local_jobdir=`printf /tmp/gen2ntuple_${sampleName}_jobid%d_%04d ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
#rm -rf $local_jobdir
mkdir -p $local_jobdir
cd $local_jobdir  


logFile="${logDir}/ntuple_maker_${outTag}_${SLURM_ARRAY_TASK_ID}.log"

maxFileCount=`wc -l < $kpsRecoFiles`
let firstfile="${SLURM_ARRAY_TASK_ID}*${nfiles}+1"
let lastfile="$firstfile+$nfiles-1"

outputs=""

echo "filest to run between first=${firstfile} to last=${lastfile}"
files=""

((iF = 0))
for n in $(seq $firstfile $lastfile); do
  if (($n > $maxFileCount)); then
    break
  fi

  suboutDir=`printf ${outDir}/arrayid_%04d ${SLURM_ARRAY_TASK_ID}`
  mkdir -p $suboutDir

  fileid=`sed -n ${n}p ${kpsRecoFiles} | awk '{ print $1 }'`
  reco_input=`sed -n ${n}p ${kpsRecoFiles} | awk '{ print $2 }'`

  let dlreco_lineno="${fileid}+1"
  dlmerged_input=`sed -n ${dlreco_lineno}p ${mdlRecoFiles}`
  dlmerged_base=`basename ${dlmerged_input}`
  reco_base=`basename ${reco_input}`

  echo "dlmerged: ${dlmerged_input}" >> ${logFile}
  echo "reco: ${reco_input}" >> ${logFile}
  echo "dlmerged (base): ${dlmerged_base}" >> ${logFile}  
  echo "reco (base): ${reco_base}" >> ${logFile}

  cp ${dlmerged_input} ${dlmerged_base}
  cp ${reco_input} ${reco_base}
  
  outFile=`printf ntuple_${outTag}_fileid%05d_output%05d_arrayid%04d_%03d.root ${fileid} ${n} ${SLURM_ARRAY_TASK_ID} ${iF}`
  echo "outFile: ${outFile}" >> ${logFile}
  cmd="make_gen2_ntuples -f $reco_base -t $dlmerged_base -m ${modelPath} -w ${weightFile} -o ${outFile} -x photon_vertex_selection --mc"
  echo $cmd >> $logFile
  $cmd >> $logFile

  cp $outFile $suboutDir/  

  outputs="$outputs $outFile"
  ((iF = iF + 1))  
done

mergedOutput="ntuple_${sampleName}_${outTag}_output_${SLURM_ARRAY_TASK_ID}.root"
echo "output files: $outputs" >> ${logFile}
echo "merging into: $mergedOutput" >> ${logFile}
hadd -f $mergedOutput $outputs
cp $mergedOutput $outDir/

# clean-up
#rm $outputs
#rm $mergedOutput
cd /tmp
#rm -r $local_jobdir

JOBENDDATE=$(date)

echo "Job began at $JOBSTARTDATE" >> $logFile
echo "Job ended at $JOBENDDATE" >> $logFile

