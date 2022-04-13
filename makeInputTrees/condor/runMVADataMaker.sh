#!/bin/bash
set -x
source /cvmfs/cms.cern.ch/cmsset_default.sh

dir=/afs/cern.ch/work/r/rchudasa/private/hiforest_1034/photonID_PbPb2018/makeInputTrees
pwd 
uname -a

#root setup
#source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-centos7-gcc48-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.24.06/x86_64-centos7-gcc48-opt/bin/thisroot.sh

inputPath=""
outputPath=""
sampleName=""

basePath="/eos/cms/store/group/phys_diffraction/lbyl_2018"
suffix="_13April22"

if [ $2 -eq 0 ]
then
  sampleName="flatPtPi0" # 10400 files
  inputPath="${basePath}/mc_flat_pt_pi0/flatPtPi0_5p02TeV_PbPb/flatPtPi0_HiForest_v2/220411_101918/0000/flatPtMC_HiForestAOD_${1}.root"
  outputPath="${basePath}/mc_flat_pt_pi0/mvaTrees${suffix}"
elif [ $2 -eq 1 ]
then
  sampleName="flatPtPhoton" # 10400 files
  inputPath="${basePath}/mc_flat_pt_photon/flatPtPhoton_5p02TeV_PbPb/flatPtPhoton_HiForest_v2/220411_102132/0000/flatPtMC_HiForestAOD_${1}.root"
  outputPath="${basePath}/mc_flat_pt_photon/mvaTrees${suffix}"
fi

mkdir -p $outputPath
output="${outputPath}/mvaTrees_${1}.root"

cp $dir/mvaDataMaker.exe .
#cp $dir/configs/Pi0Sample.cfg . 
if [ -s ${output} ]
then
  echo "File already exists, skipping"
else
  echo "Running"
  ./mvaDataMaker.exe $dir/configs/Pi0Sample.cfg $inputPath $output $3
fi 
