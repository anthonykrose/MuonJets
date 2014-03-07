MuonJets
========

MuonJet analysis with 20fb-1 


1.- Setup CMSSW release  (so far using CMSSW_5_3_11_patch5 in Brazos)

2.- Fork the repository following instructions in https://help.github.com/articles/fork-a-repo  (you need to have github account)

3.- Clone the repository to your computer

git clone https://github.com/castaned/MuonJets.git

4.- Compile the code in CMSSW 

5.- Folder Structure

    AnalysisAlgos  (Folder including the muon clustering class "MuJetProducer")
    AnalysisDataFormats (Folder including Multimuon class)
    AnalysisTools  (User Analyzer, read PAT sample and produce Ntuples for analysis)
        - CutFlowAnalyzer
               -test  
              cmsRun cutflowanalyzer_cfg.py

For further information https://help.github.com/articles/fork-a-repo





