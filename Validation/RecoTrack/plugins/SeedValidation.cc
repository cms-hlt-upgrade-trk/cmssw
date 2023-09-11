// -*- C++ -*-
//
// Package:    Validation/RecoTrack
// Class:      SeedValidation
//
/**\class SeedValidation SeedValidation.cc Validation/RecoTrack/plugins/SeedValidation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  adrianodif
//         Created:  Thu, 31 Aug 2023 07:23:36 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimTracker/TrackerHitAssociation/interface/ClusterTPAssociation.h"

#include "DataFormats/TrackReco/interface/SeedStopInfo.h"

#include "TTree.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class SeedValidation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SeedValidation(const edm::ParameterSet&);
  ~SeedValidation() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  void clearVariables();
  void fillBeamSpot(const reco::BeamSpot& bs); 

  // ----------member data ---------------------------
  
  // const edm::ESGetToken<GlobalTrackingGeometry, GlobalTrackingGeometryRecord> geoToken_;
  // const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  // const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> mfToken_;
  // const edm::ESGetToken<TransientTrackingRecHitBuilder, TransientRecHitRecord> ttrhToken_;
  const edm::EDGetTokenT<ClusterTPAssociation> tpMap_;

  edm::EDGetTokenT<edm::View<TrajectorySeed> > seedsToken;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken;
  edm::EDGetTokenT<std::vector<SeedStopInfo>> stopToken;

  TTree* t;

  // event
  edm::RunNumber_t ev_run;
  edm::LuminosityBlockNumber_t ev_lumi;
  edm::EventNumber_t ev_event;

  std::vector<int> seed_nHits;
  std::vector<float> seed_pt;
  std::vector<float> seed_tp_frac;
  std::vector<unsigned short> seed_stopReason;
  std::vector<unsigned short> seed_nCands;

  float bsp_x;
  float bsp_y;
  float bsp_z;
  float bsp_sigmax;
  float bsp_sigmay;
  float bsp_sigmaz;

  std::vector<std::vector<float>> seed_hit_x;
  std::vector<std::vector<float>> seed_hit_y;
  std::vector<std::vector<float>> seed_hit_z;
  std::vector<std::vector<int>> seed_hit_rr;
  std::vector<std::vector<int>> seed_hit_tp_id;
  std::vector<std::vector<int>> seed_hit_tp_pdg;
};

void SeedValidation::clearVariables() {
  ev_run = 0;
  ev_lumi = 0;
  ev_event = 0;

  bsp_x = 99999.9;
  bsp_y = 99999.9;
  bsp_z = 99999.9;
  bsp_sigmax = 99999.9;
  bsp_sigmay = 99999.9;
  bsp_sigmaz = 99999.9;

  //seeds
  seed_nHits.clear();
  seed_pt.clear();
  seed_tp_frac.clear();
  seed_hit_x.clear();
  seed_hit_y.clear();
  seed_hit_z.clear();
  seed_hit_rr.clear();
  seed_hit_tp_id.clear();
  seed_hit_tp_pdg.clear();
  seed_stopReason.clear();
  seed_nCands.clear();
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SeedValidation::SeedValidation(const edm::ParameterSet& iConfig)
    : 
      // geoToken_(esConsumes()),
      // tTopoToken_(esConsumes()),
      // mfToken_(esConsumes()),
      tpMap_(consumes(iConfig.getParameter<edm::InputTag>("tpMap")))
{

  // read parametes
  edm::InputTag seedsTag(iConfig.getParameter<edm::InputTag>("src"));
  edm::InputTag beamSpotTag(iConfig.getParameter<edm::InputTag>("beamSpot"));
  edm::InputTag stopTag(iConfig.getParameter<edm::InputTag>("stops"));

  //consumes
  seedsToken = consumes<edm::View<TrajectorySeed> >(seedsTag);
  beamSpotToken = consumes<reco::BeamSpot>(beamSpotTag);
  stopToken = consumes(stopTag);

  usesResource(TFileService::kSharedResource);
  edm::Service<TFileService> fs;
  t = fs->make<TTree>("seedtree", "seedtree");

  t->Branch("event", &ev_event);
  t->Branch("lumi", &ev_lumi);
  t->Branch("run", &ev_run);

  t->Branch("bsp_x", &bsp_x, "bsp_x/F");
  t->Branch("bsp_y", &bsp_y, "bsp_y/F");
  t->Branch("bsp_z", &bsp_z, "bsp_z/F");
  t->Branch("bsp_sigmax", &bsp_sigmax, "bsp_sigmax/F");
  t->Branch("bsp_sigmay", &bsp_sigmay, "bsp_sigmay/F");
  t->Branch("bsp_sigmaz", &bsp_sigmaz, "bsp_sigmaz/F");

  //seeds
  t->Branch("seed_nHits", &seed_nHits);
  t->Branch("seed_pt", &seed_pt);
  t->Branch("seed_tp_frac", &seed_tp_frac);
  t->Branch("seed_stopReason", &seed_stopReason);
  t->Branch("seed_nCands", &seed_nCands);
  t->Branch("seed_hit_x", &seed_hit_x);
  t->Branch("seed_hit_y", &seed_hit_y);
  t->Branch("seed_hit_z", &seed_hit_z);
  t->Branch("seed_hit_rr", &seed_hit_rr);
  t->Branch("seed_hit_tp_pdg", &seed_hit_tp_pdg);
  t->Branch("seed_hit_tp_pdg", &seed_hit_tp_pdg);

  //now do what ever initialization is needed
}

SeedValidation::~SeedValidation() {
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void SeedValidation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;
  using namespace std;

  //beamspot
  Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByToken(beamSpotToken, recoBeamSpotHandle);
  BeamSpot const& bs = *recoBeamSpotHandle;
  fillBeamSpot(bs);

  // input collection
  // Handle<edm::View<TrajectorySeed> > hseeds;
  // iEvent.getByToken(seedsToken, hseeds);
  const auto& seeds = iEvent.get(seedsToken);//*hseeds;
  auto const& tpClust = iEvent.get(tpMap_);
  const auto& stops = iEvent.get(stopToken);
  // create tracks from seeds

  if(stops.size()!=seeds.size())
  throw cms::Exception("LogicError") << "Got " << seeds.size() << " seeds, but " << stops.size()
                                         << " seed stopping infos" << "\n";

  for (size_t iSeed = 0; iSeed < seeds.size(); ++iSeed) {
    auto const& seed = seeds[iSeed];
    auto tsos = seed.startingState();
    auto const& stop = stops[iSeed];
    
    // // try to create a track
    // TrajectoryStateOnSurface state;
    // if (seed.nHits() == 0) {  //this is for deepCore seeds only
    //   const Surface* deepCore_sruface = &geometry_->idToDet(seed.startingState().detId())->specificSurface();
    //   state = trajectoryStateTransform::transientState(seed.startingState(), deepCore_sruface, theMF);
    // } else {
    //   TransientTrackingRecHit::RecHitPointer lastRecHit = tTRHBuilder->build(&*(seed.recHits().end() - 1));
    //   state = trajectoryStateTransform::transientState(seed.startingState(), lastrecHit.surface(), theMF);
    // }
    // TrajectoryStateClosestToBeamLine tsAtClosestApproachSeed =
    //     tscblBuilder(*state.freeState(), *beamSpot);  //as in TrackProducerAlgorithm
    // if (tsAtClosestApproachSeed.isValid()) {
    //   const reco::TrackBase::Point vSeed1(tsAtClosestApproachSeed.trackStateAtPCA().position().x(),
    //                                       tsAtClosestApproachSeed.trackStateAtPCA().position().y(),
    //                                       tsAtClosestApproachSeed.trackStateAtPCA().position().z());
    //   const reco::TrackBase::Vector pSeed(tsAtClosestApproachSeed.trackStateAtPCA().momentum().x(),
    //                                       tsAtClosestApproachSeed.trackStateAtPCA().momentum().y(),
    //                                       tsAtClosestApproachSeed.trackStateAtPCA().momentum().z());
    //   //GlobalPoint vSeed(vSeed1.x()-beamSpot->x0(),vSeed1.y()-beamSpot->y0(),vSeed1.z()-beamSpot->z0());
    //   PerigeeTrajectoryError seedPerigeeErrors =
    //       PerigeeConversions::ftsToPerigeeError(tsAtClosestApproachSeed.trackStateAtPCA());
    //   tracks->emplace_back(0., 0., vSeed1, pSeed, state.charge(), seedPerigeeErrors.covarianceMatrix());
    // } else {
    //   edm::LogVerbatim("SeedValidator") << "TrajectoryStateClosestToBeamLine not valid";
    //   // use magic values chi2<0, ndof<0, charge=0 to denote a case where the fit has failed
    //   // If this definition is changed, change also interface/trackFromSeedFitFailed.h
    //   tracks->emplace_back(
    //       -1, -1, reco::TrackBase::Point(), reco::TrackBase::Vector(), 0, reco::TrackBase::CovarianceMatrix());
    //   nfailed++;
    // }

    // tracks->back().appendHits(seed.recHits().begin(), seed.recHits().end(), ttopo);
    // // store the hits
    // size_t firsthitindex = rechits->size();

    
    auto nHits = 0;
    seed_pt.push_back(tsos.pt());
    std::vector<float> xx,yy,zz;
    std::vector<int> ii,pp,rr;
    

    struct Count {
      int clusters = 0;
      size_t innermostHit = std::numeric_limits<size_t>::max();
    };

    struct TrackTPMatch {
      int key = -1;
      int countClusters = 0;
    };

    std::unordered_map<int, Count> count;

    for (auto const& recHit : seed.recHits()) {

      
      auto range = tpClust.equal_range(dynamic_cast<const SiPixelRecHit*>(&recHit)->firstClusterRef());
        //dynamic_cast<const SiPixelRecHit*>(recHit)->firstClusterRef()).first;
      // (recHit.geographicalId().rawId());

      for (auto ip = range.first; ip != range.second; ++ip) {
        const auto tpKey = ip->second.key();
        auto& elem = count[tpKey];
        ++elem.clusters;
        elem.innermostHit = std::min(int(elem.innermostHit), int(nHits));
      }
      
      xx.push_back(recHit.globalPosition().x());
      yy.push_back(recHit.globalPosition().y());
      zz.push_back(recHit.globalPosition().z());
      rr.push_back((recHit.geographicalId().rawId()));
      

      nHits++;
    }

    TrackTPMatch best;
    int bestCount = 2;  // require >= 3 cluster for the best match
    size_t bestInnermostHit = std::numeric_limits<size_t>::max();
    for (auto& keyCount : count) {
      if (keyCount.second.clusters > bestCount ||
          (keyCount.second.clusters == bestCount && keyCount.second.innermostHit < bestInnermostHit)) {
        best.key = keyCount.first;
        best.countClusters = bestCount = keyCount.second.clusters;
        bestInnermostHit = keyCount.second.innermostHit;
      }
    }
    
    // ii.push_back((particle->second).key());
    // pp.push_back((*particle->second).pdgId());
    
    
    seed_nHits.push_back(nHits);

    // 
    // auto maxCount = 0;
    // for (auto const &i : ii)
    // {
    //   int num_items = std::count(ii.begin(), ii.end(), i);
    //   maxCount = std::max(maxCount,num_items);
    // }

    
    seed_tp_frac.push_back(float(best.countClusters)/float(nHits));
    seed_stopReason.push_back(stop.stopReasonUC());
    seed_nCands.push_back(stop.candidatesPerSeed());  
    seed_hit_x.push_back(xx);
    seed_hit_y.push_back(yy);
    seed_hit_z.push_back(zz);
    seed_hit_tp_id.push_back(ii);
    seed_hit_tp_pdg.push_back(pp);
    seed_hit_tp_pdg.push_back(rr);
    
    // // create a trackextra, just to store the hit range
    // trackextras->push_back(TrackExtra());
    // trackextras->back().setHits(ref_rechits, firsthitindex, rechits->size() - firsthitindex);
    // trackextras->back().setSeedRef(edm::RefToBase<TrajectorySeed>(hseeds, iSeed));
    // // create link between track and trackextra
    // tracks->back().setExtra(TrackExtraRef(ref_trackextras, trackextras->size() - 1));
  }

  
  t->Fill();
}


void SeedValidation::fillBeamSpot(const reco::BeamSpot& bs) {
  bsp_x = bs.x0();
  bsp_y = bs.y0();
  bsp_z = bs.x0();
  bsp_sigmax = bs.BeamWidthX();
  bsp_sigmay = bs.BeamWidthY();
  bsp_sigmaz = bs.sigmaZ();
}

// ------------ method called once each job just before starting event loop  ------------
void SeedValidation::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void SeedValidation::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SeedValidation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("initialStepSeeds"));
  desc.add<edm::InputTag>("stops", edm::InputTag("initialStepTrackCandidates"));
  desc.add<edm::InputTag>("tpMap", edm::InputTag("clusterTpMap"));
  desc.add<edm::InputTag>("beamSpot", edm::InputTag("beamSpot"));
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SeedValidation);
