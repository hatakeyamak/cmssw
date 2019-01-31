
# Quickstart

Run locally on lxplus
~~~
#set up the work area
export SCRAM_ARCH=slc6_amd64_gcc700
cmsrel CMSSW_10_4_0
cd CMSSW_10_4_0
cmsenv

#get the code
git cms-merge-topic jpata:pfvalidation-10_4_X
scram b -j4
cd $CMSSW_BASE/src/Validation/RecoParticleFlow

#make a temporary directory for the output
mkdir tmp

#RECO step, about 30 minutes
make QCD_reco

#DQM step, a few minutes
make QCD_dqm

#Do postprocessing on the DQM histograms
make QCD_post

#Do final HTML plots
make plots
~~~


# Running on condor

The reco sequence takes about 1-2 hours / 100 events on batch. We have prepared condor scripts to facilitate this on lxbatch. 
~~~
cd $CMSSW_BASE/src/Validation/RecoParticleFlow
mkdir -p tmp2/QCD/log
cd tmp2/QCD

condor_submit $CMSSW_BASE/src/Validation/RecoParticleFlow/test/condor_sub.jdl

#wait for jobs to finish, monitor using `condor_q`

#remove dummy output files from jobs that failed
du *.root | grep "^0" | awk '{print $2}' | xargs rm

cd ${CMSSW_BASE}/src/Validation/RecoParticleFlow

make QCD_dqm QCD_post plots
~~~
