#ifndef DELPHYS_ANALYSIS_TTDILEPTONANALYSER_H_
#define DELPHYS_ANALYSIS_TTDILEPTONANALYSER_H_

#include "delphes/analysis/interface/BaseAnalyser.h"

#include "classes/SortableObject.h"

#include "TH1F.h"

#include <numeric>

class TLorentzVector;


namespace pdgid {
  static const Int_t kBottom = 5;
  static const Int_t kTop = 6;
  static const Int_t kElectron = 11;
  static const Int_t kMuon = 13;
  static const Int_t kPhoton = 22;
  static const Int_t kWBoson = 24;

  static const Int_t kAntiElectron = -11;
  static const Int_t kAntiMuon = -13;

  // NOTE
  static const Int_t kWrong = std::numeric_limits<Int_t>::max();
  // In TDatabasePDG, 0 means Rootino, which indicates unidentified particle.
  static const Int_t kNeutralHadron = 0;

} // pdgid

// indices for ParticleFlow objects
// consider Make json file and then read from that file
// In deep learning framework like Keras, embedding layaer takes these indices
// as input.
namespace pfindex {
  static const Int_t kElectron = 1;
  static const Int_t kAntiElectron = 2;
  static const Int_t kMuon = 3;
  static const Int_t kAntiMuon = 4;
  static const Int_t kPositivelyChargedHadron = 5;
  static const Int_t kNegativelyChargedHadron = 6;
  static const Int_t kNeutralHadron = 7;
  static const Int_t kPhoton = 8;
} //pfindex


class TTDileptonAnalyser : private BaseAnalyser {
 public:
//  TTDileptonAnalyser(const TString & in_path, const TString & out_path, bool, float, std::string, float);
  TTDileptonAnalyser(const TString & in_path, const TString & out_path, bool, std::string, float);

  ~TTDileptonAnalyser();
  void loop();

 private:

  //
//  struct Lepton {
//    TLorentzVector P4;
//    int charge;
//    int pdgid;
//  };

  // inherited
  void makeBranch(TTree*, Bool_t);
  void resetBranch();
  void resetMemberData();
  void bookHistograms();

  Bool_t selectEvent();
  void analyse();

  void HadronReconstruction();
  void analyseJet();
  void analyseLepton();

  std::vector<const GenParticle*> getDaughters(const GenParticle*);
  Bool_t isQuark(int);
  const GenParticle* getLast(const GenParticle*);
  std::pair<const GenParticle*, const GenParticle*> getTopDecayProducts(const GenParticle*);
  std::pair<const GenParticle*, const GenParticle*> getWBosonDecayProducts(const GenParticle*);

  Bool_t matchPartonsWithJets();

  bool performTotalMinDist(const std::vector<const GenParticle*>&,
                           const std::vector<const Jet*>&);

  Int_t getPFIndex(Int_t pid, Int_t charge);
  std::pair<Int_t, Int_t> getJetPartonMatchingCode(const GenParticle*);

  //////////////////////////////////////////////////////////////////////////////
  // NOTE member data
  //////////////////////////////////////////////////////////////////////////////

  // NOTE constants
  TTree* ambig_tree_;
  const Bool_t kIsTT_;

  // functions
//  const Lepton toLepton(Electron* p){ struct Lepton l; l.P4() = p->P4(); l.charge = p->Charge; l.pdgid=-11*(p->Charge); return l;}
//  const Lepton toLepton(Muon* p){ struct Lepton l; l.P4() = p->P4(); l.charge = p->Charge; l.pdgid=-13*(p->Charge); return l;}

  //
  std::vector<const Jet*> selected_jets_;
  std::vector<const Muon*> selected_muons_;
  std::vector<const Electron*> selected_electrons_;
//  std::vector<const Lepton*> selected_leptons_;
//  std::vector<const Lepton*> selected_muons_;
//  std::vector<const Lepton*> selected_electrons_;
//
  std::vector<TLorentzVector> lep_p4_;
  std::vector<Int_t> lep_pid_;
  std::vector<Int_t> lep_charge_;

  const GenParticle* t1_;
  const GenParticle* t2_;

  const GenParticle *w1dau0_, *w1dau1_, *w2dau0_, *w2dau1_;
  const GenParticle *q1_, *w1_, *l1_, *nu1_;
  const GenParticle *q2_, *w2_, *l2_, *nu2_;
  const GenParticle *wq1_, *wq2_;

  // partons from top / anti-top quark
  std::vector<const GenParticle*> partons;
  std::vector<const GenParticle*> partons_t_;
  std::vector<const GenParticle*> partons_tbar_;

  std::map<const GenParticle*, const Jet*> parton2jet_;
  std::map<const Jet*, const GenParticle*> jet2parton_;

  // NOTE stats
  TH1F* cutflow_;
  TH1F* gen_weight_;
  Int_t cutflow_match_bin_;

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Branches
  //////////////////////////////////////////////////////////////////////////////
  // per event
  Int_t b_label_;
  Int_t b_weight_;

  Int_t b_qcd_ht_min_;

  Float_t b_met_;
  Float_t b_met_eta_;
  Float_t b_met_phi_;

  // NOTE lepton
  std::vector<TLorentzVector> b_lep_;
  std::vector<Float_t> b_lep_pt_;
  std::vector<Float_t> b_lep_eta_;
  std::vector<Float_t> b_lep_phi_;
  std::vector<Float_t> b_lep_mass_;
  std::vector<Float_t> b_lep_energy_;

  std::vector<Int_t> b_lep_pid_;
  std::vector<Int_t> b_lep_charge_;
  std::vector<Bool_t>  b_lep_isMuon_;

  // NOTE Neutrino
  std::vector<Float_t> b_nu1_pt_;
  std::vector<Float_t> b_nu1_eta_;
  std::vector<Float_t> b_nu1_phi_;
  std::vector<Float_t> b_nu2_pt_;
  std::vector<Float_t> b_nu2_eta_;
  std::vector<Float_t> b_nu2_phi_;


  // NOTE jet
  Int_t b_num_jets_; // jet multiplicity
  Int_t b_num_b_jets_; // b-jet multiplicity
  Int_t b_num_b_jets_2_; // b-jet multiplicity in the highest two jets
  Float_t b_dr_bb_; // NOTE

  //
  std::vector<TLorentzVector> b_jet_;
  std::vector<Float_t> b_jet_pt_;
  std::vector<Float_t> b_jet_eta_;
  std::vector<Float_t> b_jet_phi_;
  std::vector<Float_t> b_jet_mass_;
  std::vector<Float_t> b_jet_energy_;
  std::vector<Float_t> b_jet_lep1_dr_;
  std::vector<Float_t> b_jet_lep2_dr_;

  std::vector<Int_t> b_jet_num_chad_;
  std::vector<Int_t> b_jet_num_nhad_;
  std::vector<Int_t> b_jet_num_electron_;
  std::vector<Int_t> b_jet_num_muon_;
  std::vector<Int_t> b_jet_num_photon_;
  std::vector<Float_t> b_jet_major_axis_;
  std::vector<Float_t> b_jet_minor_axis_;
  std::vector<Float_t> b_jet_ptd_;
  std::vector<Float_t> b_jet_charge_;
  std::vector<Bool_t>  b_jet_b_tag_;
  std::vector<Int_t>   b_jet_flavor_;

  std::vector<std::vector<Float_t>> b_track_pt_;
  std::vector<std::vector<Float_t>> b_track_eta_;
  std::vector<std::vector<Float_t>> b_track_phi_;
  std::vector<std::vector<Float_t>> b_track_deta_;
  std::vector<std::vector<Float_t>> b_track_dphi_;
  std::vector<std::vector<Float_t>> b_track_x_;
  std::vector<std::vector<Float_t>> b_track_ptrel_;
  std::vector<std::vector<Float_t>> b_track_pzrel_;
  std::vector<std::vector<Float_t>> b_track_d0_;
  std::vector<std::vector<Float_t>> b_track_dz_;
  std::vector<std::vector<Float_t>> b_track_delta_x_;
  std::vector<std::vector<Float_t>> b_track_delta_y_;
  std::vector<std::vector<Float_t>> b_track_delta_z_;
  std::vector<std::vector<Float_t>> b_track_pid_;

  std::vector<std::vector<Float_t>> b_tower_pt_;
  std::vector<std::vector<Float_t>> b_tower_eta_;
  std::vector<std::vector<Float_t>> b_tower_deta_;
  std::vector<std::vector<Float_t>> b_tower_dphi_;
  std::vector<std::vector<Float_t>> b_tower_x_;
  std::vector<std::vector<Float_t>> b_tower_ptrel_;
  std::vector<std::vector<Float_t>> b_tower_pzrel_;
  std::vector<std::vector<Float_t>> b_tower_d0_;
  std::vector<std::vector<Float_t>> b_tower_dz_;
  std::vector<std::vector<Float_t>> b_tower_E_;
  std::vector<std::vector<Float_t>> b_tower_Eem_;
  std::vector<std::vector<Float_t>> b_tower_Ehad_;



  // NOTE constituents of jet

  // NOTE kIsTT_
  Bool_t b_is_matched_; // TT only, Jet-Parton matching
  Bool_t b_is_tt_; // TT only, Jet-Parton matching

  // FIXME b_jet_parton_match_;
  // 0: other / 1: t / 2: tbar
  std::vector<Int_t> b_jet_parton_match_;

  // 0: other / 1: b/bbar / 2: s/sbar
  std::vector<Int_t> b_jet_parton_match_detail_;

  std::vector<Int_t> b_jet_parton_match_pid_;

  // 0: other / 1: t / -1: tbar
  std::vector<Int_t> b_jet_parton_match_t_;

  // 0: other / 1: b/bbar / -1: s/sbar
  std::vector<Int_t> b_jet_parton_match_b_;

  // NOTE parton
  Int_t   b_l1_pid_, b_nu1_pid_;
  Int_t   b_l2_pid_, b_nu2_pid_;

  TLorentzVector b_t1_, b_q1_, b_w1_, b_l1_, b_nu1_;
  TLorentzVector b_t2_, b_q2_, b_w2_, b_l2_, b_nu2_;

  // NOTE https://github.com/cms-sw/cmssw/blob/CMSSW_11_2_0_pre2/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L120-L132
  std::vector<Float_t> b_parton_jet_distance_;
  Float_t b_sum_delta_energy_;
  Float_t b_sum_delta_pt_;
  Float_t b_sum_delta_r_;

  // Di-lepton channel
  // 0: ee / 1:emu / 2:mumu
  Int_t b_ch_;
  Int_t b_gen_ch_;



  // RecoHadron
  struct RecoHad {
      TLorentzVector tlv;
      TLorentzVector dau1_tlv;
      TLorentzVector dau2_tlv;
      int idx;
      int dau1_idx;
      int dau2_idx;
      int pdgid;
      int dau1_pdgid;
      int dau2_pdgid;
      int isFrom;
      bool isPreSelected;
      bool isJetMatched;
  };

  std::vector<struct RecoHad> m_recoHad;
  
  //////////////////////////////////////////////////////////////////////////////
  // NOTE selection cut
  // Measurement of the top quark mass in the all-jets final state at
  // \sqrt{s} = 13 TeV and combination with the lepton+jets channel
  // https://link.springer.com/article/10.1140/epjc/s10052-019-6788-2
  //////////////////////////////////////////////////////////////////////////////
  // NOTE Object Selection
  const Float_t kMinJetPT_  = 30.0f;
  const Float_t kMaxJetEta_ = 2.4f;

  const Float_t kMinLeptonPT_  = 20.0f;
  const Float_t kMinFirstLeptonPT_  = 25.0f;
  const Float_t kMaxLeptonEta_ = 2.4f;

  const Float_t kElectronIso_ = 0.12f;
  const Float_t kMuonIso_ = 0.15f;

  const Float_t kMinMET_ = 40.0f;
  const Float_t kMinDileptonMass_ = 20.0f;

  const Float_t kZbosonMass_ = 91.1876f;

  // NOTE HLT
  const UInt_t  kHLTMinNumJets_  = 2;
  const UInt_t  kHLTMinNumBJets_ = 1; // b-tagged jets
  const Float_t  kJetOverlap_ = 0.4f;

  // NOTE Jet Parton Matching
  // https://github.com/cms-sw/cmssw/blob/master/TopQuarkAnalysis/TopTools/python/TtFullHadJetPartonMatch_cfi.py#L28-L35
  const Bool_t kUseMaxDist_ = true;
  const Float_t kMaxDist_ = 0.4f;
  const std::string kAlgorithm_;
  const Float_t kExponent_;

  // ECAL Super Cluster  cut
  const Float_t cut_lowElecEffRegionUp = 1.57f;
  const Float_t cut_lowElecEffRegionDn = 1.44f;
  const Float_t cut_ElectronEtaBarrel  = 1.48f;
  const Float_t cut_ElectronRelIso03All_inBarrel = 0.0588f;
  const Float_t cut_ElectronRelIso03All_inEndCap = 0.0571f;

  // HadronReconstruction cut
  const Float_t tkPTCut_ = 0.95f;
  const Float_t tkEtaCut_ = 2.4f;
  const Float_t tkIPSigXYCut_ = 2.0f;

  // Mass
  const Float_t pionMass_ = 0.13967f;
  const Float_t protonMass_ = 0.938272f;
  const Float_t ksMass_ = 0.49761f;
  const Float_t lambda0Mass_ = 1.11568f;
  const Int_t pionPdgId_ = 211;
  const Int_t protonPdgId_ = 2212;
  const Int_t ksPdgId_ = 310;
  const Int_t lambda0PdgId_ = 3122;



};



#endif //  DELPHYS_ANALYSIS_TTDILEPTONSANALYSER_H_
