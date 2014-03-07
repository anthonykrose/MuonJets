// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// user include files
#include "TTree.h"
#include "TRandom3.h"
#include "TMath.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/TriggerReport.h"
#include "FWCore/Framework/interface/TriggerNamesService.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSegmentMatch.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "AnalysisDataFormats/MuJetAnalysis/interface/MultiMuon.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"


//******************************************************************************
//              Auxiliary function: Order objects by pT                         
//******************************************************************************
bool PtOrder (const reco::GenParticle* p1, const reco::GenParticle* p2) { return (p1->pt() > p2->pt() ); }
bool PtOrderrec (const reco::Muon* p1, const reco::Muon* p2) { return (p1->pt() > p2->pt() ); }

//******************************************************************************
// Auxiliary function: Calculate difference between two angles: -PI < phi < PI  
//******************************************************************************
double My_dPhi (double phi1, double phi2) {
  double dPhi = phi1 - phi2;
  if (dPhi >  M_PI) dPhi -= 2.*M_PI;
  if (dPhi < -M_PI) dPhi += 2.*M_PI;
  return dPhi;
}

//******************************************************************************
//                           Class declaration                                  
//******************************************************************************

class Analysis : public edm::EDAnalyzer {
public:
  explicit Analysis(const edm::ParameterSet&);
  ~Analysis();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:

  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      
  // ---------- Member data ----------
  
  TTree * AnaTree;  // Pointer to Tree
  
  Int_t run;    // run number    
  Int_t lumi;   // lumi number   
  Int_t event;  // event number  
  Int_t trigger;
  Int_t isVtx;

  // ---------- Generator Level ----------
  
  edm::InputTag m_genParticles;   // Label to access generator particles
  

  // Gen branches in ROOT tree

  // ---------- GEN level ---------------

  Int_t genmuons;
  Int_t gengamD;
  Int_t genHiggs;

  Float_t ptgenMuons[10];
  Float_t etagenMuons[10];

  Float_t ptgenMuonsA0[5];
  Float_t ptgenMuonsA1[5];
  Float_t etagenMuonsA0[5];
  Float_t etagenMuonsA1[5];
  Float_t ptgengamD[5];
  Float_t ptgenHiggs[5];

  // ---------- RECO Level ----------
  
  Int_t recmuons;
  Float_t ptrecMuons[30];
  Float_t etarecMuons[30];

  Int_t recmujets;
  Int_t mj1muons;
  Int_t mj2muons;
  Float_t ptmj1muons[4];
  Float_t ptmj2muons[4];
  Float_t etamj1muons[4];
  Float_t etamj2muons[4];

  Float_t dzmj1;
  Float_t dzmj2;

  Float_t massmj1;
  Float_t massmj2;

  Float_t isomj1;
  Float_t isomj2;

  // Labels to access
  edm::InputTag m_muons;  // reconstructed muons
  edm::InputTag m_muJets; // muon jets built from reconstructed muons

  // Auxiliary variables
  TRandom3 m_trandom3;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
Analysis::Analysis(const edm::ParameterSet& iConfig)

{

  
  genmuons=0;
  gengamD=0;
  genHiggs=0;
  recmujets=0;
  mj1muons=0;
  mj2muons=0;
  genmuons=0;
  recmuons=0;
  dzmj1=0.;
  dzmj2=0.;
  massmj1=0.0;
  massmj2=0.0;
  isomj1=0.0;
  isomj2=0.0;
  isVtx=0;

  AnaTree  = NULL;    
  
  // ---------- Generator Level ----------
  m_genParticles = iConfig.getParameter<edm::InputTag>("genParticles");
  
  
  // ---------- RECO Level ----------
  m_muons = iConfig.getParameter<edm::InputTag>("muons");
  m_muJets = iConfig.getParameter<edm::InputTag>("muJets");


  //  m_maxIsoDiMuons  = iConfig.getParameter<double>("maxIsoDiMuons");

}


Analysis::~Analysis()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}

//
// member functions
//

// ------------ method called for each event  ------------
void
Analysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  // get the run, lumi and event number
  //  run   = iEvent.id().run();
  //  lumi  = iEvent.id().luminosityBlock();
  event = iEvent.id().event();


  //  std::cout<<" event  "<<event<<std::endl;

  //****************************************************************************
  //                              GEN LEVEL                                     
  //****************************************************************************
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  iEvent.getByLabel(m_genParticles, genParticles);
  
  // Loop over all genParticles and save prompt muons from particles with codes 36 (a1) or 3000022 (gammaD) in vector genMuons
  std::vector<const reco::GenParticle*> genH;
  std::vector<const reco::GenParticle*> genA;
  std::vector<const reco::GenParticle*> genMuons;
  std::vector<const reco::Candidate*>   genMuonMothers;

  // Loop over all gen particles
  int counterGenParticle = 0;
  for(reco::GenParticleCollection::const_iterator iGenParticle = genParticles->begin();  iGenParticle != genParticles->end();  ++iGenParticle) {
    counterGenParticle++;
    
    if ( fabs( iGenParticle->pdgId() ) == 13 && iGenParticle->status() == 1 ) {
      // Mother of the muon can be muon. Find the last muon in this chain: genMuonCand
      // Example: a1 -> mu+ (status = 3) mu- (status = 3)
      //          mu- (status = 3) -> mu- (status = 2) -> mu- (status = 1)
      const reco::Candidate *genMuonCand = &(*iGenParticle);
      bool isMuonMother = true;
      while(isMuonMother) {
        isMuonMother = false;
        for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
          if ( fabs( genMuonCand->mother(iMother)->pdgId() ) == 13 ) {
            isMuonMother = true;
            genMuonCand = genMuonCand->mother(iMother);
          }
        }
      }
      // Loop over all real (non-muon) mothers of the muon (here we use genMuonCand)
      for ( size_t iMother = 0; iMother < genMuonCand->numberOfMothers(); iMother++ ) {
        // Check if mother is CP-odd Higgs (PdgId = 36) or gamma_Dark (PdgId = 3000022)
        if ( genMuonCand->mother(iMother)->pdgId() == 36 || genMuonCand->mother(iMother)->pdgId() == 3000022 || genMuonCand->mother(iMother)->pdgId() == 443 ) {
          // Store the muon (stable, first in chain) into vector
          genMuons.push_back(&(*iGenParticle));
          // Store mother of the muon into vector. We need this to group muons into dimuons later
          genMuonMothers.push_back(genMuonCand->mother(iMother));
        }
      }
    }

    // Check if gen particle is decaying (status = 3) CP-even Higgs (pdgId = +/-35)
    if ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 35 ) {
      genH.push_back(&(*iGenParticle)); // Store the Higgs into vector
    }
    // Check if gen particle is
    if (    ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 36      )     // decaying (status = 3) CP-odd Higgs (pdgId = 36)
	    || ( iGenParticle->status() == 3 && iGenParticle->pdgId() == 3000022 )     // decaying (status = 3) gamma_Dark (pdgId = 3000022)
	    || ( iGenParticle->status() == 2 && iGenParticle->pdgId() == 443     ) ) { // decaying (status = 2) J/psi (pdgId = 443)
      genA.push_back(&(*iGenParticle));
    }
  }
  

  genmuons = genMuons.size();
  gengamD  = genA.size();
  genHiggs = genH.size();


  std::sort( genMuons.begin(), genMuons.end(), PtOrder );



  //  std::cout<<" gen mouns "<<genmuons<<std::endl;

  for(int k=0;k<genmuons;k++){
    ptgenMuons[k] = genMuons[k]->pt();
    etagenMuons[k] = genMuons[k]->eta();
  }


  if ( genH.size() == 1 ) {
    ptgenHiggs[0] = genH[0]->pt();
  }
  
  if ( genA.size() >= 2 ) {
    // Sort genA by pT (leading pT first)
    std::sort (genA.begin(), genA.end(), PtOrder);
    ptgengamD[0] = genA[0]->pt();
    ptgengamD[1] = genA[1]->pt();
  } else {
    std::cout << "WARNING! genA.size() < 2" << std::endl;
  }
  
  // Group muons from the same mother into dimuons
  std::vector< std::vector<const reco::GenParticle*> > genMuonGroups;
  std::vector<const reco::GenParticle*> genMuonsTMP1       = genMuons;
  std::vector<const reco::Candidate*>   genMuonMothersTMP1 = genMuonMothers;
  unsigned int nMuonGroup = 0;
  while ( genMuonsTMP1.size() > 0 ) {
    std::vector<const reco::GenParticle*> genMuonsTMP2;
    std::vector<const reco::Candidate*>   genMuonMothersTMP2;
    std::vector<const reco::GenParticle*> genMuonsSameMother;
    for ( unsigned int j = 0; j < genMuonsTMP1.size(); j++ ) {
      // Check if mothers are the same particle
      if ( fabs( genMuonMothersTMP1[0]->pt() - genMuonMothersTMP1[j]->pt() ) < 0.00001 ) {
	genMuonsSameMother.push_back( genMuonsTMP1[j] );
      } else {
	genMuonsTMP2.push_back( genMuonsTMP1[j] );
	genMuonMothersTMP2.push_back( genMuonMothersTMP1[j] );
      }
    }
    genMuonGroups.push_back(genMuonsSameMother);
    genMuonsTMP1       = genMuonsTMP2;
    genMuonMothersTMP1 = genMuonMothersTMP2;
    nMuonGroup++;
  }


  //  std::cout<<" gen mouns groups  "<<genMuonGroups.size()<<std::endl;
  
  if ( genMuonGroups.size() == 2 && genMuonGroups[0].size() == 2 && genMuonGroups[1].size() == 2 ) {
    std::sort( genMuonGroups[0].begin(), genMuonGroups[0].end(), PtOrder );
    std::sort( genMuonGroups[1].begin(), genMuonGroups[1].end(), PtOrder );

    ptgenMuonsA0[0] = genMuonGroups[0][0]->pt();
    ptgenMuonsA0[1] = genMuonGroups[0][1]->pt();

    ptgenMuonsA1[0] = genMuonGroups[1][0]->pt();
    ptgenMuonsA1[1] = genMuonGroups[1][1]->pt();

    etagenMuonsA0[0] = genMuonGroups[0][0]->eta();
    etagenMuonsA0[1] = genMuonGroups[0][1]->eta();

    etagenMuonsA1[0] = genMuonGroups[1][0]->eta();
    etagenMuonsA1[1] = genMuonGroups[1][1]->eta();
  }


  
   //================= Reco =======================//
  
  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByLabel(m_muons, muons);

  std::vector<const reco::Muon*> recMuons;

  for (pat::MuonCollection::const_iterator iMuon = muons->begin();  iMuon != muons->end();  ++iMuon) {
    if ( fabs(iMuon->eta()) < 2.4 && ( iMuon->isTrackerMuon() || iMuon->isGlobalMuon() ) ) {
      recMuons.push_back( &(*iMuon));
    }
  }

  recmuons = recMuons.size();


  //  std::cout<<" number of recmuons  "<<recmuons<<std::endl;

  // Sort recMuons by pT (leading pT first)
  if ( recMuons.size() > 1 ) std::sort( recMuons.begin(), recMuons.end(), PtOrderrec );


  for(int k=0;k<recmuons;k++){
    ptrecMuons[k] = recMuons[k]->pt();
    etarecMuons[k] = recMuons[k]->eta();
  }
  

   //================ MuJets  ======================//

   edm::Handle<pat::MultiMuonCollection> muJets;
   iEvent.getByLabel(m_muJets, muJets);

   const pat::MultiMuon *muJet1 = NULL;
   const pat::MultiMuon *muJet2 = NULL;

   unsigned int nMuJets = muJets->size();
   recmujets = muJets->size();

   //   std::cout<<" number of mujets  "<<recmujets<<std::endl;


   if (nMuJets == 2) {

     mj1muons = (*muJets)[0].numberOfDaughters();
     mj2muons = (*muJets)[1].numberOfDaughters();
      
     muJet1 = &((*muJets)[0]);
     muJet2 = &((*muJets)[1]);

     //     std::cout<<" muons in MuJet1  "<<mj1muons<<std::endl;
     //     std::cout<<" muons in MuJet2  "<<mj2muons<<std::endl;
     for(int m=0;m<mj1muons;m++){
       ptmj1muons[m]  = (*muJets)[0].muon(m)->pt();
       etamj1muons[m]  = (*muJets)[0].muon(m)->eta();
     }

     for(int m=0;m<mj2muons;m++){
       ptmj2muons[m]  = (*muJets)[1].muon(m)->pt();
       etamj2muons[m]  = (*muJets)[1].muon(m)->eta();
     }


   }
  

   //   std::cout<<"  beam spot   "<<std::endl;
   edm::Handle<reco::BeamSpot> beamSpot;
   iEvent.getByLabel("offlineBeamSpot",beamSpot);
//    //  std::cout << "Beam spot x: " << beamSpot->position().x() << " y: " << beamSpot->position().y() << " z: " << beamSpot->position().z() << std::endl;
  


   if ( muJet1 != NULL && muJet2 != NULL ) {
     //     std::cout<<" muon jet dz calculation  "<<std::endl;

     dzmj1 = muJet1->dz(beamSpot->position());
     dzmj2 = muJet2->dz(beamSpot->position());
   }
  
    edm::Handle<pat::TriggerEvent> triggerEvent;
    iEvent.getByLabel("patTriggerEvent", triggerEvent);
  
    //    std::cout<<"  trigger   "<<std::endl;
    bool isDiMuonHLTFired = false;
    if (    ( triggerEvent->path("HLT_Mu17_Mu8_v22") && triggerEvent->path("HLT_Mu17_Mu8_v22")->wasAccept() )
         || ( triggerEvent->path("HLT_Mu17_Mu8_v21") && triggerEvent->path("HLT_Mu17_Mu8_v21")->wasAccept() )
         || ( triggerEvent->path("HLT_Mu17_Mu8_v19") && triggerEvent->path("HLT_Mu17_Mu8_v19")->wasAccept() )
         || ( triggerEvent->path("HLT_Mu17_Mu8_v18") && triggerEvent->path("HLT_Mu17_Mu8_v18")->wasAccept() )
         || ( triggerEvent->path("HLT_Mu17_Mu8_v17") && triggerEvent->path("HLT_Mu17_Mu8_v17")->wasAccept() )
         || ( triggerEvent->path("HLT_Mu17_Mu8_v16") && triggerEvent->path("HLT_Mu17_Mu8_v16")->wasAccept() ) ) {
      isDiMuonHLTFired = true;
    }  

    if(isDiMuonHLTFired) trigger=1;
    else trigger=0;

    //    std::cout<<"  mass   "<<std::endl;
    if ( muJet1!= NULL && muJet2!= NULL ) {
      //      std::cout<<" muon jet mass calculation  "<<std::endl;
      massmj1= muJet1->mass();
      massmj2= muJet2->mass();
    }

   edm::Handle<reco::TrackCollection> tracks;
   iEvent.getByLabel("generalTracks", tracks);
  
   if ( muJet1 != NULL && muJet2 != NULL ) {

     const pat::MultiMuon *diMuonTmp = NULL;
     for ( unsigned int i = 1; i <= 2; i++ ) { 
       double isoTmp = 0.0;
       if ( i == 1 ) diMuonTmp = muJet1;
       if ( i == 2 ) diMuonTmp = muJet2;
       for (reco::TrackCollection::const_iterator track = tracks->begin(); track != tracks->end(); ++track) {
 	bool trackIsMuon = false;
 	if ( diMuonTmp->sameTrack( &*track, &*(diMuonTmp->muon(0)->innerTrack()) ) || diMuonTmp->sameTrack( &*track, &*(diMuonTmp->muon(1)->innerTrack()) ) ) trackIsMuon = true;
 	if (!trackIsMuon) {
 	  double dPhi = My_dPhi( diMuonTmp->phi(), track->phi() );
 	  double dEta = diMuonTmp->eta() - track->eta();
 	  double dR = sqrt( dPhi*dPhi + dEta*dEta ); 
 	  if ( dR < 0.4 && track->pt() > 0.5 ) {
 	    double dz = fabs( track->dz(beamSpot->position()) - diMuonTmp->dz(beamSpot->position()) );
 	    if ( dz < 0.1 ) isoTmp += track->pt();
 	  }    
 	}
       }
       if ( i == 1 ) isomj1 = isoTmp;
       if ( i == 2 ) isomj2 = isoTmp;
     }
   }  

  
   //    std::cout<<"vertex   "<<std::endl;
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByLabel("offlinePrimaryVertices", primaryVertices);
  

    bool isVtxOk=false;
    for (reco::VertexCollection::const_iterator vertex = primaryVertices->begin();  vertex != primaryVertices->end();  ++vertex) {
      if (vertex->isValid() && !vertex->isFake() && vertex->tracksSize() > 3 && fabs(vertex->z()) < 24.) isVtxOk=true;
    }

    if(isVtxOk) isVtx=1;
  

    //    std::cout<<" Fill TRee   "<<std::endl;
    AnaTree->Fill();
    //   std::cout<<"after  Fill TRee   "<<std::endl;


}


// ------------ method called once each job just before starting event loop  ------------
void 
Analysis::beginJob()
{
  std::cout << "BEGIN JOB" << std::endl;
  
  edm::Service<TFileService> tFileService;
  AnaTree = tFileService->make<TTree>("Events", "Events");


  // Event variables
  AnaTree->Branch("event",     &event,     "event/I");
  AnaTree->Branch("trigger",   &trigger,   "trigger/I");
  AnaTree->Branch("isVtx",     &isVtx,     "isVtx/I");
  AnaTree->Branch("genmuons",  &genmuons,  "genmuons/I");
  AnaTree->Branch("gengamD",   &gengamD,   "gengamD/I");
  AnaTree->Branch("genHiggs",  &genHiggs,  "genHiggs/I");
  AnaTree->Branch("recmuons",  &recmuons,  "recmuons/I");
  AnaTree->Branch("recmujets", &recmujets, "recmujets/I");

  // Generator Muons
  AnaTree->Branch("ptgengamD0",      &ptgengamD,      "ptgengamD[2]/F");
  AnaTree->Branch("ptgenHiggs",      &ptgenHiggs,     "ptgenHiggs[0]/F");

  AnaTree->Branch("ptgenMuons",     &ptgenMuons,    "ptgenMuons[genmuons]/F");
  AnaTree->Branch("etagenMuons",    &etagenMuons,   "etagenMuons[genmuons]/F");

  AnaTree->Branch("ptgenMuonsA0",    &ptgenMuonsA0,   "ptgenMuonsA0[2]/F");
  AnaTree->Branch("ptgenMuonsA1",    &ptgenMuonsA1,   "ptgenMuonsA1[2]/F");
  AnaTree->Branch("etagenMuonsA0",   &etagenMuonsA0,  "etagenMuonsA0[2]/F");
  AnaTree->Branch("etagenMuonsA1",   &etagenMuonsA1,  "etagenMuonsA1[2]/F");
  

    // Reconstructed Muons

    AnaTree->Branch("ptrecMuons",  &ptrecMuons,    "ptrecMuons[recmuons]/F");
    AnaTree->Branch("etarecMuons", &etarecMuons,   "etarecMuons[recmuons]/F");

    AnaTree->Branch("mj1muons",    &mj1muons,  "mj1muons/I");
    AnaTree->Branch("mj2muons",    &mj2muons,  "mj2muons/I");

    AnaTree->Branch("ptmj1muons",  &ptmj1muons, "ptmj1muons[mj1muons]/F");
    AnaTree->Branch("ptmj2muons",  &ptmj2muons, "ptmj2muons[mj2muons]/F");
    AnaTree->Branch("etamj1muons", &etamj1muons, "etamj1muons[mj1muons]/F");
    AnaTree->Branch("etamj2muons", &etamj2muons, "etamj2muons[mj2muons]/F");

    AnaTree->Branch("dzmj1",       &dzmj1,  "dzmj1/F");
    AnaTree->Branch("dzmj2",       &dzmj2,  "dzmj2/F");

    AnaTree->Branch("massmj1",     &massmj1, "massmj1/F");
    AnaTree->Branch("massmj2",     &massmj2, "massmj2/F");

   AnaTree->Branch("isomj1",     &isomj1, "isomj1/F");
   AnaTree->Branch("isomj2",     &isomj2, "isomj2/F");

}

// ------------ method called once each job just after ending the event loop  ------------
void 
Analysis::endJob() 
{
  
}

// ------------ method called when starting to processes a run  ------------
void 
Analysis::beginRun(edm::Run const & iRun, edm::EventSetup const & iSetup)
{
  
}

// ------------ method called when ending the processing of a run  ------------
void 
Analysis::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
Analysis::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
Analysis::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analysis);
