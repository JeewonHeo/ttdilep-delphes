#ifndef DELPHYS_ANALYSIS_BASEANALYSER_H_
#define DELPHYS_ANALYSIS_BASEANALYSER_H_

#include "classes/DelphesClasses.h"

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TString.h"

#include <set>
#include <iostream>

// TODO use enum??
// enum class DelphesBranches {
//   kEvent, kParticle, kGenJet, kGenMissingET,
//   kTrack, kTower, kEFlowTrack, kEFlowPhoton, kEFlowNeutralHadron,
//   kElectron, kPhoton, kMuon, kJet, kFatJet, kMissingET, kScalarHT, kVertex,
// };


class BaseAnalyser {
 public:
  BaseAnalyser();

  // constructor to clone Delphes' tree
  BaseAnalyser(const TString & in_path,
               const TString & out_path);

  BaseAnalyser(const TString & in_path,
               const TString & out_path,
               const TString & out_tree_name);

  // constructor for multiple input files
  BaseAnalyser(const std::vector<TString> & in_paths,
               const TString & out_path,
               const TString & out_tree_name);

  ~BaseAnalyser();

  void initDelphesBranch();
  void setBranchAddress(std::set<TString> branches, Bool_t drop=false);
  void setBranchAddress();

  inline const GenParticle* getParticle(Int_t idx) {
    return dynamic_cast<const GenParticle*>(particles_->At(idx));
  }

  inline const Track* getTrack(Int_t idx) {
    return dynamic_cast<const Track*>(tracks_->At(idx));
  }

  inline const Tower* getTower(Int_t idx) {
    return dynamic_cast<const Tower*>(towers_->At(idx));
  }

  inline const Track* getEFlowTrack(Int_t idx) {
    return dynamic_cast<const Track*>(eflow_tracks_->At(idx));
  }

  inline const Tower* getEFlowNeutralHadron(Int_t idx) {
    return dynamic_cast<const Tower*>(eflow_neutral_hadrons_->At(idx));
  }

  inline const Tower* getEFlowPhoton(Int_t idx) {
    return dynamic_cast<const Tower*>(eflow_photons_->At(idx));
  }

  inline const Jet* getGenJet(Int_t idx) {
    return dynamic_cast<const Jet*>(gen_jets_->At(idx));
  }

  inline const MissingET* getGenMissingET(Int_t idx=0) {
    return dynamic_cast<const MissingET*>(gen_mets_->At(idx));
  }

  inline const Jet* getJet(Int_t idx) {
    return dynamic_cast<const Jet*>(jets_->At(idx));
  }

  inline const Electron* getElectron(Int_t idx) {
    return dynamic_cast<const Electron*>(electrons_->At(idx));
  }

  inline const Muon* getMuon(Int_t idx) {
    return dynamic_cast<const Muon*>(muons_->At(idx));
  }

  inline const Photon* getPhoton(Int_t idx) {
    return dynamic_cast<const Photon*>(photons_->At(idx));
  }

  inline const Jet* getFatJet(Int_t idx) {
    return dynamic_cast<const Jet*>(fat_jets_->At(idx));
  }

  inline const MissingET* getMissingET(Int_t idx=0) {
    return dynamic_cast<const MissingET*>(mets_->At(idx));
  }

  inline const ScalarHT* getScalarHT(Int_t idx=0) {
    return dynamic_cast<const ScalarHT*>(scalar_hts_->At(idx));
  }



  // Member data

  TFile* in_file_;
  TFile* out_file_;
  TTree* in_tree_;
  TTree* out_tree_;

  TClonesArray* events_;
  TClonesArray* particles_; // GenParticle
  TClonesArray* gen_jets_;  // 
  TClonesArray* gen_mets_;
  TClonesArray* tracks_;
  TClonesArray* towers_;
  TClonesArray* eflow_tracks_;
  TClonesArray* eflow_photons_;
  TClonesArray* eflow_neutral_hadrons_;
  TClonesArray* electrons_;
  TClonesArray* muons_;
  TClonesArray* photons_;
  TClonesArray* jets_;
  TClonesArray* fat_jets_;
  TClonesArray* mets_;
  TClonesArray* scalar_hts_;
  TClonesArray* vertices_;

  const std::set<TString> kAllDelphesBranches_ = {
      "Event", "Particle", "GenJet", "GenMissingET",
      "Track", "Tower", "EFlowTrack", "EFlowPhoton", "EFlowNeutralHadron",
      "Electron", "Photon", "Muon",
      "Jet", "FatJet", "MissingET", "ScalarHT", "Vertex"
  };

};


#endif // DELPHYS_ANALYSIS_BASEANALYSER_H_
