#!/bin/bash

JOBSTARTDATE=$(date)

scriptDir="/path/to/repo"
pythonScript="${scriptDir}/$1"

kpsRecoFiles=$2
mdlRecoFiles=$3

weightFile="${scriptDir}/event_weighting/$4"

checkpointDir="/cluster/tufts/wongjiradlabnu/mrosen25/prongCNN/models/checkpoints"
modelPath="${checkpointDir}/$5"

outTag=$6

nfiles=$7

echo "running array ID $SLURM_ARRAY_TASK_ID (sample output tag: $outTag) on node $SLURMD_NODENAME"

ubdlDir="/path/to/ubdl"
source ${ubdlDir}/setenv_py3.sh
source ${ubdlDir}/configure.sh
export PYTHONPATH=${PYTHONPATH}:${scriptDir}

maxFileCount=`wc -l < $kpsRecoFiles`
let firstfile="${SLURM_ARRAY_TASK_ID}*${nfiles}+1"
let lastfile="$firstfile+$nfiles-1"
files=""

for n in $(seq $firstfile $lastfile); do
  if (($n > $maxFileCount)); then
    break
  fi
  newfile=`sed -n ${n}p ${kpsRecoFiles}`
  files="$files $newfile"
done

scriptName=`echo $1 | sed s/.py//g`
outDir="/path/to/output/directory"
logDir="/path/to/logfile/directory"
logFile="${logDir}/${scriptName}_${outTag}_${SLURM_ARRAY_TASK_ID}.log"

local_jobdir=`printf /tmp/run_selection_jobid%d_%04d ${SLURM_JOB_ID} ${SLURM_ARRAY_TASK_ID}`
rm -rf $local_jobdir
mkdir -p $local_jobdir
cd $local_jobdir  

((iF = 0))
outputs=""

for file in $files; do
  
  outFile="${scriptName}_${outTag}_output_${SLURM_ARRAY_TASK_ID}_${iF}.root"
  echo "inputfile path: $file" >> ${logFile}
  echo "outFile: $outFile" >> ${logFile}

  if (($# > 8)); then
    echo "python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath $8 $9 >> $logFile"
    python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath $8 $9 >> $logFile
  elif  (($# > 7)); then
    echo "python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath $8 >> $logFile"
    python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath $8 >> $logFile
  else
    echo "python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath >> $logFile"
    python3 $pythonScript -f $file -t $mdlRecoFiles -o $outFile -w $weightFile -m $modelPath >> $logFile
  fi

  ((iF = iF + 1))
  outputs="$outputs $outFile"

done

mergedOutput="${scriptName}_${outTag}_output_${SLURM_ARRAY_TASK_ID}.root"
echo "output files: $outputs" >> ${logFile}
echo "merging into: $mergedOutput" >> ${logFile}
hadd $mergedOutput $outputs
cp $mergedOutput $outDir

# clean-up
rm $outputs
rm $mergedOutput
cd /tmp
rm -r $local_jobdir

JOBENDDATE=$(date)

echo "Job began at $JOBSTARTDATE" >> $logFile
echo "Job ended at $JOBENDDATE" >> $logFile

