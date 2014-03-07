import FWCore.ParameterSet.Config as cms

maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
				   
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_10_1_ijK.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_11_1_iVO.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_1_1_hs6.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_12_1_GQE.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_13_1_axo.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_14_1_AbS.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_15_1_WdY.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_16_1_WVb.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_2_1_WHr.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_3_1_5fl.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_4_1_Y71.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_5_1_m4b.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_6_1_ezv.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_7_1_fRX.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_8_1_vGO.root",
"file:/fdata/hepx/store/user/pakhotin/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_GEN_v2/DarkSUSY_mH_125_mGammaD_0400_ctauExp_2_8TeV-madgraph452_bridge224_LHE_pythia6_537p4_PAT_v1/output_9_1_Y0U.root",


   ]);

secFiles.extend( [
    ] )

