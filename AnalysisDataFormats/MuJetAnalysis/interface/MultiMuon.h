//
// $Id: MultiMuon.h,v 1.12 2012/05/01 10:00:52 aysen Exp $
//
// Jim Pivarski <pivarski@physics.tamu.edu>
// 

#ifndef AnalysisDataFormats_MuJetAnalysis_MultiMuon_h
#define AnalysisDataFormats_MuJetAnalysis_MultiMuon_h

/**
  \class    pat::MultiMuon MultiMuon.h "AnalysisDataFormats/MuJetAnalysis/interface/MultiMuon.h"
  \brief    Analysis-level particle class

   MultiMuon implements an analysis-level multi-muon class within the 'pat'
   namespace.

  \version  $Id: MultiMuon.h,v 1.12 2012/05/01 10:00:52 aysen Exp $
*/

// uncomment for compiling in strict FWLite
// #define MULTIMUONCANDIDATE_FOR_FWLITE

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "TMath.h"
#include "TTree.h"

#ifdef MULTIMUONCANDIDATE_FOR_FWLITE
typedef int TransientTrackBuilder;
#endif
#ifndef MULTIMUONCANDIDATE_FOR_FWLITE
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#endif

// Define typedefs for convenience
namespace pat {
   class MultiMuon;
   typedef std::vector<MultiMuon>              MultiMuonCollection; 
   typedef edm::Ref<MultiMuonCollection>       MultiMuonRef; 
   typedef edm::RefVector<MultiMuonCollection> MultiMuonRefVector; 
}

// Class definition
namespace pat {
   class MultiMuon : public pat::CompositeCandidate {
      public:
	 /// default constructor
	 MultiMuon() : pat::CompositeCandidate()
		     , m_vertexFitted(false)
		     , m_chi2(0.)
		     , m_ndof(0.)
		     , m_centralTrackIsolationCone(0.)
		     , m_unionTrackIsolationCone(0.)
		     , m_centralTrackThresholdPt(0.)
		     , m_unionTrackThresholdPt(0.)
		     , m_centralCaloIsolationCone(0.)
		     , m_unionCaloIsolationCone(0.)
		     , m_centralNumberAboveThresholdCone(0.)
		     , m_unionNumberAboveThresholdCone(0.)
		     , m_centralNumberAboveThresholdPt(0.)
		     , m_unionNumberAboveThresholdPt(0.)
		     , m_centralTrackIsolation(0.)
		     , m_unionTrackIsolation(0.)
		     , m_centralECALIsolation(0.)
		     , m_unionECALIsolation(0.)
		     , m_centralHCALIsolation(0.)
		     , m_unionHCALIsolation(0.)
		     , m_centralNumberAboveThreshold(0)
		     , m_unionNumberAboveThreshold(0) {};
	 
	 MultiMuon(double phi) : pat::CompositeCandidate()
			       , m_vertexFitted(false)
			       , m_chi2(0.)
			       , m_ndof(0.)
			       , m_centralTrackIsolationCone(0.)
			       , m_unionTrackIsolationCone(0.)
			       , m_centralTrackThresholdPt(0.)
			       , m_unionTrackThresholdPt(0.)
			       , m_centralCaloIsolationCone(0.)
			       , m_unionCaloIsolationCone(0.)
			       , m_centralNumberAboveThresholdCone(0.)
			       , m_unionNumberAboveThresholdCone(0.)
			       , m_centralNumberAboveThresholdPt(0.)
			       , m_unionNumberAboveThresholdPt(0.)
			       , m_centralTrackIsolation(0.)
			       , m_unionTrackIsolation(0.)
			       , m_centralECALIsolation(0.)
			       , m_unionECALIsolation(0.)
			       , m_centralHCALIsolation(0.)
			       , m_unionHCALIsolation(0.)
			       , m_centralNumberAboveThreshold(0)
			       , m_unionNumberAboveThreshold(0) { phi_ = phi; }
	 
	 /// constructor with muons
	 MultiMuon(std::vector<const pat::Muon*> &muons, const TransientTrackBuilder *transientTrackBuilder = NULL, const reco::TrackCollection *tracks = NULL, const pat::MuonCollection *allmuons = NULL, const CaloTowerCollection *caloTowers = NULL, double centralTrackIsolationCone = 0., double unionTrackIsolationCone = 0., double centralTrackThresholdPt = 1e6, double unionTrackThresholdPt = 1e6, double centralCaloIsolationCone = 0., double unionCaloIsolationCone = 0., double centralNumberAboveThresholdCone = 0., double unionNumberAboveThresholdCone = 0., double centralNumberAboveThresholdPt = 1e6, double unionNumberAboveThresholdPt = 1e6);
	 
	 /// constructor from a composite candidate
	 MultiMuon(const pat::MultiMuon & aMultiMuon);
	 
	 /// destructor
	 virtual ~MultiMuon();
	 
	 /// required reimplementation of the Candidate's clone method
	 virtual MultiMuon * clone() const { return new MultiMuon(*this); }
	 
	 /// cast daughters as MultiMuons
	 const pat::Muon *muon(int i) const { return dynamic_cast<const pat::Muon*>(daughter(i)); }
	 
	 /// calculate a vertex from the daughter muons (performed by constructor if transientTrackBuilder != NULL)
	 bool calculateVertex(const TransientTrackBuilder *transientTrackBuilder);
	 
	 /// calculate isolation (performed by constructor if tracks, muons, and caloTowers != NULL)
	 void calculateTrackIsolation(const reco::TrackCollection *tracks, const pat::MuonCollection *allmuons, double centralCone, double unionCone, double centralThreshold, double unionThreshold, TTree *diagnosticTTree = NULL, Float_t *diagnosticdR = NULL, Float_t *diagnosticpT = NULL);
	 void calculateCaloIsolation(const CaloTowerCollection *caloTowers, double centralCone, double unionCone);
	 void calculateNumberAboveThresholdIsolation(const reco::TrackCollection *tracks, const pat::MuonCollection *allmuons, double centralCone, double unionCone, double centralThreshold, double unionThreshold, TTree *diagnosticTTree = NULL, Float_t *diagnosticdR = NULL, Float_t *diagnosticpT = NULL);
	 
	 /// does this MultiMuon overlap another one? or contain a given muon?
	 bool overlaps(const pat::MultiMuon &aMultiMuon) const;
	 bool contains(const pat::Muon &aMuon) const;
	 
	 /// create a new MultiMuon which has muons from this and aMultiMuon
	 pat::MultiMuon merge(const pat::MultiMuon &aMultiMuon, const TransientTrackBuilder *transientTrackBuilder = NULL, const reco::TrackCollection *tracks = NULL, const pat::MuonCollection *allmuons = NULL, const CaloTowerCollection *caloTowers = NULL, double centralTrackIsolationCone = 0., double unionTrackIsolationCone = 0., double centralTrackThresholdPt = 1e6, double unionTrackThresholdPt = 1e6, double centralCaloIsolationCone = 0., double unionCaloIsolationCone = 0., double centralNumberAboveThresholdCone = 0., double unionNumberAboveThresholdCone = 0., double centralNumberAboveThresholdPt = 1e6, double unionNumberAboveThresholdPt = 1e6);
	 
	 /// return vertex results
	 bool vertexValid() const { return m_vertexFitted; }
	 double vertexChi2() const { checkVertex();  return m_chi2; }
	 double vertexNdof() const { checkVertex();  return m_ndof; }
	 double vertexNormalizedChi2() const { checkVertex();  return (m_ndof > 0. ? m_chi2/m_ndof : 0.); }
	 double vertexProb() const { checkVertex();  return (m_ndof > 0. ? TMath::Prob(m_chi2, m_ndof) : 0.); }
	 CovarianceMatrix vertexCovariance() const { checkVertex();  return m_covarianceMatrix; }
	 double vertexCovariance(int i, int j) const { checkVertex();  return m_covarianceMatrix.At(i, j); }
	 
	 /// return position/momentum of each muon closest to vertex
	 GlobalPoint vertexPCA(int i) const { checkVertex();  return m_vertexPCA[i]; }
	 CovarianceMatrix vertexPCACovarianceMatrix(int i) const { checkVertex();  return m_vertexPCACovarianceMatrix[i]; }
	 GlobalVector vertexMomentum(int i) const { LorentzVector v = vertexP4(i);  return GlobalVector(v.x(), v.y(), v.z()); }
	 LorentzVector vertexP4(int i) const { checkVertex();  return m_vertexP4[i]; }
	 
	 /// return position/momentum of multimuon object at vertex
	 /// for position, use virtual const Point& reco::LeafCandidate::vertex();
	 GlobalPoint vertexPoint() const { checkVertex();  Point v = vertex();  return GlobalPoint(v.x(), v.y(), v.z()); }
	 GlobalVector vertexMomentum() const { LorentzVector v = vertexP4();  return GlobalVector(v.x(), v.y(), v.z()); }
	 LorentzVector vertexP4() const {
	    checkVertex();
	    LorentzVector v;
	    for (unsigned int i = 0;  i < numberOfDaughters();  i++) {
	       v += vertexP4(i);
	    }
	    return v;
	 }
	 
	 double vertexMass() const { return vertexP4().mass(); }

	double dz(const Point& myBeamSpot) const {
		if (vertexValid()) {
			GlobalPoint v = vertexPoint();
			GlobalVector p = vertexMomentum();
			double pt = pow(p.x()*p.x()+p.y()*p.y(),0.5);
			return (v.z()-myBeamSpot.z()) - ((v.x()-myBeamSpot.x())*p.x()+(v.y()-myBeamSpot.y())*p.y())/pt * p.z()/pt;
		}
		else return muon(0)->innerTrack()->dz(myBeamSpot);
	}
		 

	 /// return the distance of flight from the primary vertex in the direction of momentum
	 /// in 2D
	 double lxy(GlobalPoint primaryVertex) const {
	    GlobalPoint v = vertexPoint();
	    return ((v.x()-primaryVertex.x())*px() + (v.y()-primaryVertex.y())*py())/pt();
	 };
	 double lxy(double x, double y, double z) const {
	    return lxy(GlobalPoint(x, y, z));
	 };
	 /// in 3D
	 double lxyz(GlobalPoint primaryVertex) const {
	    GlobalPoint v = vertexPoint();
	    return ((v.x()-primaryVertex.x())*px() + (v.y()-primaryVertex.y())*py() + (v.z()-primaryVertex.z())*pz())/p();
	 };
	 double lxyz(double x, double y, double z) const {
	    return lxyz(GlobalPoint(x, y, z));
	 };

	 /// return daughter's momentum vector in the COM frame
	 GlobalVector daughterCOM(int i, bool vertex = false) const;
	 /// return daughter's momentum vector in the COM frame, rotated to the boost axis (+z is parent's direction)
	 GlobalVector daughterCOMrot(int i, bool vertex = false) const;
	 /// return daughter's cos(theta) where theta is the angle between the daughter's momentum and the parent's momentum, in the center-of-mass frame
	 double daughterCOMcosTheta(int i, bool vertex = false) const;

	 /// return opening angles with or without the vertex correction
	 double dphi(int i, int j, bool vertex = false) const {
	    double phi_i, phi_j;
	    if (vertex) {
	       checkVertex();
	       phi_i = vertexP4(i).phi();
	       phi_j = vertexP4(j).phi();
	    }
	    else {
	       phi_i = daughter(i)->phi();
	       phi_j = daughter(j)->phi();
	    }
	    double delta = phi_i - phi_j;
	    while (delta > M_PI) delta -= 2.*M_PI;
	    while (delta < -M_PI) delta += 2.*M_PI;
	    return delta;
	 }
	 double deta(int i, int j, bool vertex = false) const {
	    double eta_i, eta_j;
	    if (vertex) {
	       checkVertex();
	       eta_i = vertexP4(i).eta();
	       eta_j = vertexP4(j).eta();
	    }
	    else {
	       eta_i = daughter(i)->eta();
	       eta_j = daughter(j)->eta();
	    }
	    return eta_i - eta_j;
	 }
	 double dR(int i, int j, bool vertex = false) const {
	    return sqrt(pow(dphi(i, j, vertex), 2) + pow(deta(i, j, vertex), 2));
	 }
	 double dRmax(bool vertex = false) const {
	    double max = 0.;
	    for (unsigned int i = 0;  i < numberOfDaughters();  i++) {
	       for (unsigned int j = i+1;  j < numberOfDaughters();  j++) {
		  double dR_ij = dR(i, j, vertex);
		  if (dR_ij > max) max = dR_ij;
	       }
	    }
	    return max;
	 }

	 /// return the combination of final-state pairs most consistent with being the same mass
	 /// (assumes all cascades decay to muon pairs through a single lowest-mass resonance)
	 /// return value is a set of pairs of indices: the first is the mu+, the second is the mu- of each pair
	 std::vector<std::pair<int,int> > consistentPairs(bool vertex = false) const;
	 /// return just the masses of the most consistent combination
	 std::vector<double> consistentPairMasses(bool vertex = false) const;

	 /// return isolation information
	 double centralTrackIsolationCone() const { return m_centralTrackIsolationCone; }
	 double unionTrackIsolationCone() const { return m_unionTrackIsolationCone; }
	 double centralTrackThresholdPt() const { return m_centralTrackThresholdPt; }
	 double unionTrackThresholdPt() const { return m_unionTrackThresholdPt; }
	 double centralCaloIsolationCone() const { return m_centralCaloIsolationCone; }
	 double unionCaloIsolationCone() const { return m_unionCaloIsolationCone; }
	 double centralNumberAboveThresholdCone() const { return m_centralNumberAboveThresholdCone; }
	 double unionNumberAboveThresholdCone() const { return m_unionNumberAboveThresholdCone; }
	 double centralNumberAboveThresholdPt() const { return m_centralNumberAboveThresholdPt; }
	 double unionNumberAboveThresholdPt() const { return m_unionNumberAboveThresholdPt; }
	 
	 double centralTrackIsolation() const { return m_centralTrackIsolation; }
	 double unionTrackIsolation() const { return m_unionTrackIsolation; }
	 double centralECALIsolation() const { return m_centralECALIsolation; }
	 double unionECALIsolation() const { return m_unionECALIsolation; }
	 double centralHCALIsolation() const { return m_centralHCALIsolation; }
	 double unionHCALIsolation() const { return m_unionHCALIsolation; }
	 double centralCaloIsolation() const { return m_centralECALIsolation + m_centralHCALIsolation; }
	 double unionCaloIsolation() const { return m_unionECALIsolation + m_unionHCALIsolation; }
	 int centralNumberAboveThreshold() const { return m_centralNumberAboveThreshold; }
	 int unionNumberAboveThreshold() const { return m_unionNumberAboveThreshold; }
	 
	 bool sameTrack(const reco::Track *one, const reco::Track *two) const;

      protected:
	 void checkVertex() const;
	 double noiseEcal(const CaloTower &tower) const;
	 double noiseHcal(const CaloTower &tower) const;
	 double noiseHOcal(const CaloTower &tower) const;
	 void buildPermutation(std::vector<std::vector<int> > &results, std::vector<int> working, int where, int value) const;

	 bool m_vertexFitted;
	 double m_chi2;
	 double m_ndof;
	 CovarianceMatrix m_covarianceMatrix;
	 std::vector<GlobalPoint> m_vertexPCA;
	 std::vector<CovarianceMatrix> m_vertexPCACovarianceMatrix;
	 std::vector<LorentzVector> m_vertexP4;
	 
	 double m_centralTrackIsolationCone;
	 double m_unionTrackIsolationCone;
	 double m_centralTrackThresholdPt;
	 double m_unionTrackThresholdPt;
	 double m_centralCaloIsolationCone;
	 double m_unionCaloIsolationCone;
	 double m_centralNumberAboveThresholdCone;
	 double m_unionNumberAboveThresholdCone;
	 double m_centralNumberAboveThresholdPt;
	 double m_unionNumberAboveThresholdPt;
	 double m_centralTrackIsolation;
	 double m_unionTrackIsolation;
	 double m_centralECALIsolation;
	 double m_unionECALIsolation;
	 double m_centralHCALIsolation;
	 double m_unionHCALIsolation;
	 int m_centralNumberAboveThreshold;
	 int m_unionNumberAboveThreshold;
   };
}

#endif // AnalysisDataFormats_MuJetAnalysis_MultiMuon_h
