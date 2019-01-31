TMPDIR=tmp
RELVALCMD=${CMSSW_BASE}/src/Validation/RecoParticleFlow/test/run_relval.sh

QCD: QCD_reco QCD_dqm

QCD_reco:
	cd ${TMPDIR} && ${RELVALCMD} QCD reco 0

QCD_dqm:
	rm -f ${TMPDIR}/QCD/DQM*.root
	cd ${TMPDIR} && ${RELVALCMD} QCD dqm 0

QCD_post:
	rm -f ${TMPDIR}/QCD/*withResponse.root ${TMPDIR}/QCD/DQMTotal.root
	cd ${TMPDIR} && python ${CMSSW_BASE}/src/Validation/RecoParticleFlow/test/addResponse.py QCD/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO.root
	mv ${TMPDIR}/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_withResponse.root ${TMPDIR}/QCD/DQM_V0001_R000000001__Global__CMSSW_X_Y_Z__RECO_withResponse.root
	hadd -f ${TMPDIR}/QCD/DQMTotal.root ${TMPDIR}/QCD/DQM_*.root
 
plots:
	rm -Rf plots	
	python test/compare.py --sample FlatQCD:${TMPDIR}/QCD/DQMTotal.root:${TMPDIR}/QCD/DQMTotal.root \
		--doResponsePlots

.PHONY: QCD_reco QCD_dqm QCD_post plots
