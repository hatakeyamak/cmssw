#!/bin/bash


cmsRun EDAnalyzers/TreeMaker/python/ConfFile_cfg.py \
    sourceFile=\
sourceFiles/SingleElectron_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1_FEVT/SingleElectron_PT2to200_Phase2HLTTDRSummer20ReRECOMiniAOD-PU200_111X_mcRun4_realistic_T15_v1_ext2-v1_FEVT.txt \
    genEleFilter=0 \
    isGunSample=1 \
    onRaw=1 \
    debugFile=1 \
    maxEvents=20 \
