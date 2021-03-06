#include "MLAnalyzer/RecHitAnalyzer/interface/RecHitAnalyzer.h"
#include <iostream>
using namespace std;

// Fill EB rec hits ////////////////////////////////
// Store event rechits in a vector of length equal
// to number of crystals in EB (ieta:170 x iphi:360)

TProfile2D *hEB_energy;
TProfile2D *hEB_time;
std::vector<float> vEB_energy_;
std::vector<float> vEB_time_;

// Initialize branches _____________________________________________________//
void RecHitAnalyzer::branchesEB ( TTree* tree, edm::Service<TFileService> &fs ) {

  // Branches for images
  tree->Branch("EB_energy", &vEB_energy_);
  tree->Branch("EB_time",   &vEB_time_);

  // Histograms for monitoring
  hEB_energy = fs->make<TProfile2D>("EB_energy", "E(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );
  hEB_time = fs->make<TProfile2D>("EB_time", "t(i#phi,i#eta);i#phi;i#eta",
      EB_IPHI_MAX  , EB_IPHI_MIN-1, EB_IPHI_MAX,
      2*EB_IETA_MAX,-EB_IETA_MAX,   EB_IETA_MAX );

} // branchesEB()

// Fill EB rechits _________________________________________________________________//
void RecHitAnalyzer::fillEB ( const edm::Event& iEvent, const edm::EventSetup& iSetup ) {

  int iphi_, ieta_, idx_; // rows:ieta, cols:iphi
  float energy_;

  vEB_energy_.assign( EBDetId::kSizeForDenseIndexing, 0. );
  vEB_time_.assign( EBDetId::kSizeForDenseIndexing, 0. );

  edm::Handle<EcalRecHitCollection> EBRecHitsH_;
  iEvent.getByToken( EBRecHitCollectionT_, EBRecHitsH_);

  // Fill EB rechits 
  for ( EcalRecHitCollection::const_iterator iRHit = EBRecHitsH_->begin();
        iRHit != EBRecHitsH_->end(); ++iRHit ) {

    energy_ = iRHit->energy();
    if ( energy_ <= zs ) continue;
    // Get detector id and convert to histogram-friendly coordinates
    EBDetId ebId( iRHit->id() );
    iphi_ = ebId.iphi() - 1;
    ieta_ = ebId.ieta() > 0 ? ebId.ieta()-1 : ebId.ieta();
    // Fill histograms for monitoring 
    hEB_energy->Fill( iphi_,ieta_,energy_ );
    hEB_time->Fill( iphi_,ieta_,iRHit->time() );
    // Get Hashed Index: provides convenient 
    // index mapping from [ieta][iphi] -> [idx]
    idx_ = ebId.hashedIndex(); // (ieta_+EB_IETA_MAX)*EB_IPHI_MAX + iphi_
    // Fill vectors for images
    vEB_energy_[idx_] = energy_;
    vEB_time_[idx_] = iRHit->time();

  } // EB rechits
std::cout<<"idx_0: "<<vEB_energy_[0]<<" "<<"idx_1: "<<vEB_energy_[1]<<" "<<"idx_2: "<<vEB_energy_[2]<<" "<<"idx_10: "<<vEB_energy_[10]<<" "<<"idx_11: "<<vEB_energy_[11]<<" "<<"idx_12: "<<vEB_energy_[12]<<" "<<"idx_100: "<<vEB_energy_[100]<<" "<<"idx_101: "<<vEB_energy_[101]<<" "<<"idx_102: "<<vEB_energy_[102]<<" "<<"idx_1000: "<<vEB_energy_[1000]<<" "<<"idx_1001: "<<vEB_energy_[1001]<<" "<<"idx_1002: "<<vEB_energy_[1002]<<" "<<"idx_61197: "<<vEB_energy_[61197]<<" "<<"idx_61198: "<<vEB_energy_[61198]<<" "<<"idx_61199: "<<vEB_energy_[61199]<<" "<<" -> size is "<<vEB_energy_.size()<<endl;
} // fillEB()
