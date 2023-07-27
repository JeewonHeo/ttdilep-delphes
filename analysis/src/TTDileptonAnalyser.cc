#include "delphes/analysis/interface/TTDileptonAnalyser.h"
#include "delphes/external/interface/JetDiscrimination.h"

#include "TString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include "TDatabasePDG.h"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <memory>
#include <cassert> // assert


using namespace delphes;


TTDileptonAnalyser::TTDileptonAnalyser(const TString & in_path,
                                           const TString & out_path,
                                           bool is_tt,
//                                           bool use_max_dist,
//                                           float max_dist,
                                           std::string chargeAlgorithm,
                                           float exponent)
    : BaseAnalyser(in_path, out_path, "delphes"),
      kIsTT_(is_tt),
//      kUseMaxDist_(use_max_dist),
//      kMaxDist_(max_dist),
      kAlgorithm_(chargeAlgorithm),
      kExponent_(exponent) {

  std::cout << "ctor begin" << std::endl;

  std::cout << "input path = " << in_path << std::endl;
  std::cout << "output path = " << out_path << std::endl;
  std::cout << "is_tt = " << (is_tt ? "true" : "false") << std::endl;
  std::cout << "use_max_dist = " << (kUseMaxDist_ ? "true"  : "false") << std::endl;
  std::cout << "max_dist = " << kMaxDist_ << std::endl;
  std::cout << "JetChargeAlgorithm = " << kAlgorithm_ << std::endl;
  std::cout << "exponent = " << kExponent_ << std::endl;

  if (kAlgorithm_ == "MomentumWeighted") {
    if (kExponent_ < 0.0) {
      auto err_msg = Form("[ValueError] kAlgorithm_: %s, kExponent_: %.2f", kAlgorithm_.c_str(), kExponent_);
      std::cerr << err_msg << std::endl;
      std::exit(1);
    }

  } else if (kAlgorithm_ == "EnergyWeighted") {
    if (kExponent_ < 0.0) {
      auto err_msg = Form("[ValueError] kAlgorithm_: %s, kExponent_: %.2f", kAlgorithm_.c_str(), kExponent_);
      std::cerr << err_msg << std::endl;
      std::exit(1);
    }

  } else {
    auto err_msg = Form("[ValueError] unknown algotiyhm %s", kAlgorithm_.c_str());
    std::cerr << err_msg << std::endl;
    std::exit(1);
  }

  std::cout << "setBranchAddress begin" << std::endl;
  setBranchAddress({"Vertex"}, /*drop=*/true);
  std::cout << "setBranchAddress end" << std::endl;

  std::cout << "makeBranch begin" << std::endl;
  // NOTE
  makeBranch(out_tree_, kIsTT_);
  std::cout << "makeBranch end" << std::endl;

  if (kIsTT_) {
    std::cout << "tt" << std::endl;

    ambig_tree_ = new TTree("unmatched", "unmatched");
    ambig_tree_->SetDirectory(out_file_);
    makeBranch(ambig_tree_, false);
    //ambig_tree_ = nullptr;

  } else {
    ambig_tree_ = nullptr;
    std::cout << "bkg" << std::endl;
  }

  bookHistograms();
  std::cout << "ctor end" << std::endl;
}


TTDileptonAnalyser::~TTDileptonAnalyser() {
  std::cout << "dtor begin" << std::endl;

  out_tree_->Write();
  if (kIsTT_) {
    ambig_tree_->Write();
  }

  std::cout << "out: " << out_tree_->GetEntries() << std::endl;
  if (kIsTT_) {
    std::cout << "unmatched: " << ambig_tree_->GetEntries() << std::endl;
  }

  out_file_->Write();
  out_file_->Close();

  std::cout << "dtor end" << std::endl;
}


#define BRANCH_(name, suffix) tree->Branch(#name, & b_##name##_ , #name "/" #suffix);
#define BRANCH_I(name) BRANCH_(name, I);
#define BRANCH_F(name) BRANCH_(name, F);
#define BRANCH_O(name) BRANCH_(name, O);

#define BRANCH_A_(name, size, suffix) tree->Branch(#name, & b_##name##_ , #name"["#size"]/"#suffix);
#define BRANCH_AF(name, size)  BRANCH_A_(name, size, F);

#define BRANCH_V_(name, type) tree->Branch(#name, "vector<"#type">", & b_##name##_ );
#define BRANCH_VF(name) BRANCH_V_(name, float);
#define BRANCH_VI(name) BRANCH_V_(name, Int_t);
#define BRANCH_VO(name) BRANCH_V_(name, Bool_t);
#define BRANCH_VVF(name) BRANCH_V_(name, std::vector<float>);

#define BRANCH_TLV(name) tree->Branch(#name, "TLorentzVector", & b_##name##_);
#define BRANCH_ANY(name) tree->Branch(#name, & b_##name##_ );

#define BRANCH_GEN_INFO(name) BRANCH_F(name##_pt); \
                              BRANCH_F(name##_eta); \
                              BRANCH_F(name##_phi); \
                              BRANCH_F(name##_mass); \
                              BRANCH_F(name##_energy); \
                              BRANCH_TLV(name);


void TTDileptonAnalyser::makeBranch(TTree* tree, bool is_matched) {
  resetBranch();

  BRANCH_I(label);
  BRANCH_I(weight);
  BRANCH_I(ch);
  BRANCH_I(gen_ch);
//  BRANCH_O(offline);

  // NOTE MET
  BRANCH_F(met);
  BRANCH_F(met_eta);
  BRANCH_F(met_phi);

  // NOTE Lepton
  BRANCH_ANY(lep);
  BRANCH_VF(lep_pt);
  BRANCH_VF(lep_eta);
  BRANCH_VF(lep_phi);
  BRANCH_VF(lep_mass);
  BRANCH_VF(lep_energy);

  BRANCH_VI(lep_pid);
  BRANCH_VI(lep_charge);
  BRANCH_VO(lep_isMuon);

  // NOTE Jet
  // jet variables
  BRANCH_I(num_jets);
  BRANCH_I(num_b_jets);
  BRANCH_I(num_b_jets_2);

  BRANCH_F(dr_bb);

  BRANCH_ANY(jet);
  BRANCH_VF(jet_pt);
  BRANCH_VF(jet_eta);
  BRANCH_VF(jet_phi);
  BRANCH_VF(jet_mass);
  BRANCH_VF(jet_energy);
  BRANCH_VF(jet_lep1_dr);
  BRANCH_VF(jet_lep2_dr);

  BRANCH_VF(nu1_pt);
  BRANCH_VF(nu1_eta);
  BRANCH_VF(nu1_phi);

  BRANCH_VF(nu2_pt);
  BRANCH_VF(nu2_eta);
  BRANCH_VF(nu2_phi);

  BRANCH_VI(jet_num_chad);
  BRANCH_VI(jet_num_nhad);
  BRANCH_VI(jet_num_electron);
  BRANCH_VI(jet_num_muon);
  BRANCH_VI(jet_num_photon);
  BRANCH_VF(jet_major_axis);
  BRANCH_VF(jet_minor_axis);
  BRANCH_VF(jet_ptd);
  BRANCH_VF(jet_charge);

  BRANCH_VO(jet_b_tag); // my own, dR = 0.3
  BRANCH_VI(jet_flavor); // NOTE from delphes dR=0.5

  BRANCH_VVF(track_pt);
  BRANCH_VVF(track_eta);
  BRANCH_VVF(track_phi);
  BRANCH_VVF(track_deta);
  BRANCH_VVF(track_dphi);
  BRANCH_VVF(track_x);
  BRANCH_VVF(track_ptrel);
  BRANCH_VVF(track_pzrel);
  BRANCH_VVF(track_d0);
  BRANCH_VVF(track_dz);
  BRANCH_VVF(track_delta_x);
  BRANCH_VVF(track_delta_y);
  BRANCH_VVF(track_delta_z);
  BRANCH_VVF(track_pid);

  BRANCH_VVF(tower_pt);
  BRANCH_VVF(tower_eta);
  BRANCH_VVF(tower_deta);
  BRANCH_VVF(tower_dphi);
  BRANCH_VVF(tower_x);
  BRANCH_VVF(tower_ptrel);
  BRANCH_VVF(tower_pzrel);
  BRANCH_VVF(tower_E);
  BRANCH_VVF(tower_Eem);
  BRANCH_VVF(tower_Ehad);

  if (is_matched) {
    BRANCH_O(is_matched);
    BRANCH_O(is_tt);

    BRANCH_VI(jet_parton_match);
    BRANCH_VI(jet_parton_match_detail);
    BRANCH_VI(jet_parton_match_pid);
    BRANCH_VI(jet_parton_match_t);
    BRANCH_VI(jet_parton_match_b);

    BRANCH_TLV(t1);
    BRANCH_TLV(q1);
    BRANCH_TLV(w1);
    BRANCH_TLV(l1);
    BRANCH_TLV(nu1);

    BRANCH_TLV(t2);
    BRANCH_TLV(q2);
    BRANCH_TLV(w2);
    BRANCH_TLV(l2);
    BRANCH_TLV(nu2);

    BRANCH_I(l1_pid);
    BRANCH_I(nu1_pid);
    BRANCH_I(l2_pid);
    BRANCH_I(nu2_pid);

    BRANCH_VF(parton_jet_distance);
    BRANCH_F(sum_delta_energy);
    BRANCH_F(sum_delta_pt);
    BRANCH_F(sum_delta_r);
  } else{
    BRANCH_O(is_tt);
    BRANCH_VI(jet_parton_match);
    BRANCH_VI(jet_parton_match_detail);
  }
  return ;
}


void TTDileptonAnalyser::resetBranch() {
  b_label_ = -1;
  b_weight_ = -99;
  b_ch_ = -1;
//  b_offline_ = false;

  b_met_     = -1.0f;
  b_met_eta_ = 100.0f;
  b_met_phi_ = 100.0f;

  b_num_jets_     = -1;
  b_num_b_jets_   = -1;
  b_num_b_jets_2_ = -1;

  b_dr_bb_ = -1.0f;

  // lepton
  b_lep_.clear();
  b_lep_pt_.clear();
  b_lep_eta_.clear();
  b_lep_phi_.clear();
  b_lep_mass_.clear();
  b_lep_energy_.clear();

  b_lep_pid_.clear();
  b_lep_charge_.clear();
  b_lep_isMuon_.clear();

  // neutrino
  b_nu1_pt_.clear();
  b_nu1_eta_.clear();
  b_nu1_phi_.clear();
  b_nu2_pt_.clear();
  b_nu2_eta_.clear();
  b_nu2_phi_.clear();

  // quark
  b_jet_.clear();
  b_jet_pt_.clear();
  b_jet_eta_.clear();
  b_jet_phi_.clear();
  b_jet_mass_.clear();
  b_jet_energy_.clear();
  b_jet_lep1_dr_.clear();
  b_jet_lep2_dr_.clear();

  b_jet_num_chad_.clear();
  b_jet_num_nhad_.clear();
  b_jet_num_electron_.clear();
  b_jet_num_muon_.clear();
  b_jet_num_photon_.clear();
  b_jet_major_axis_.clear();
  b_jet_minor_axis_.clear();
  b_jet_ptd_.clear();
  b_jet_charge_.clear();

  b_jet_b_tag_.clear();
  b_jet_flavor_.clear();

  b_track_pt_.clear();
  b_track_eta_.clear();
  b_track_phi_.clear();
  b_track_deta_.clear();
  b_track_dphi_.clear();
  b_track_x_.clear();
  b_track_ptrel_.clear();
  b_track_pzrel_.clear();
  b_track_d0_.clear();
  b_track_dz_.clear();
  b_track_delta_x_.clear();
  b_track_delta_y_.clear();
  b_track_delta_z_.clear();
  b_track_pid_.clear();

  b_tower_pt_.clear();
  b_tower_eta_.clear();
  b_tower_deta_.clear();
  b_tower_dphi_.clear();
  b_tower_x_.clear();
  b_tower_ptrel_.clear();
  b_tower_pzrel_.clear();
  b_tower_E_.clear();
  b_tower_Eem_.clear();
  b_tower_Ehad_.clear();

  if (kIsTT_) {
    b_is_matched_ = false;

    b_jet_parton_match_.clear();
    b_jet_parton_match_detail_.clear();
    b_jet_parton_match_pid_.clear();
    b_jet_parton_match_b_.clear();
    b_jet_parton_match_t_.clear();

    b_t1_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_q1_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_w1_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_l1_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_nu1_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);

    b_t2_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_q2_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_w2_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_l2_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);
    b_nu2_.SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0);

    b_l1_pid_ = 0;
    b_nu1_pid_ = 0;
    b_l2_pid_ = 0;
    b_nu2_pid_ = 0;

    b_parton_jet_distance_.clear();
    b_sum_delta_energy_ = -999.0f;
    b_sum_delta_pt_ = -999.0f;
    b_sum_delta_r_ = -999.0f;


  } else {
    b_jet_parton_match_.clear();
    b_jet_parton_match_detail_.clear();
  }
}


void TTDileptonAnalyser::bookHistograms() {
  gen_weight_ = new TH1F("genWeight", "genWeight", 1, 0, 1);
  cutflow_ = new TH1F("cutflow", "CutFlow", 7, 0.5, 7.5);

  TAxis* x_axis = cutflow_->GetXaxis();
  x_axis->SetBinLabel(1, "Generated");
  x_axis->SetBinLabel(2, Form("N_{leptons} == 2"));
  x_axis->SetBinLabel(3, Form("veto Z boson"));
  x_axis->SetBinLabel(4, Form("MET cut #geq %.1f GeV", kMinMET_));
  x_axis->SetBinLabel(5, Form("[HLT] N_{jets} #geq %d", kHLTMinNumJets_));
  x_axis->SetBinLabel(6, Form("[HLT] N_{b-jets} #geq %d", kHLTMinNumBJets_));

  std::string jet_parton_match_label = "[JPM] TotalMinDist";
  if (kUseMaxDist_)
    jet_parton_match_label += Form(" (MaxDist=%.3f)", kMaxDist_);
  else
    jet_parton_match_label += " (no MaxDist)";

  x_axis->SetBinLabel(7, jet_parton_match_label.c_str());
}


void TTDileptonAnalyser::resetMemberData() {
  selected_jets_.clear();
  selected_electrons_.clear();
  selected_muons_.clear();

  lep_p4_.clear();
  lep_pid_.clear();
  lep_charge_.clear();

  if (kIsTT_) {
    t1_ = nullptr;
    t2_ = nullptr;

    q1_ = nullptr;
    w1_ = nullptr;
    l1_ = nullptr;
    nu1_ = nullptr;

    q2_ = nullptr;
    w2_ = nullptr;
    l2_ = nullptr;
    nu2_ = nullptr;

    // w1dau0_ = nullptr;
    // w1dau1_ = nullptr;
    // w2dau0_ = nullptr;
    // w2dau1_ = nullptr;
    // wq1_ = nullptr;
    // wq2_ = nullptr;

    parton2jet_.clear();
    jet2parton_.clear();
  }
}


std::vector<const GenParticle*> TTDileptonAnalyser::getDaughters(const GenParticle* mother) {
  std::vector<const GenParticle*> daughters;
  if (mother->D1 == -1) {
    return daughters;
  }

  for (Int_t dau_idx = mother->D1; dau_idx <= mother->D2; dau_idx++) {
    auto dau = getParticle(dau_idx);
    daughters.push_back(dau);
  }

  return daughters;
}


const GenParticle* TTDileptonAnalyser::getLast(const GenParticle* particle) {
  const GenParticle* last = particle;

  while (true) {
    auto daughters = getDaughters(last);
    const unsigned int num_daughters = daughters.size();

    if (num_daughters == 1) {
      auto dau = daughters.front();
      if (last->PID == dau->PID) {
        last = dau;
        continue;

      } else {
        std::cerr << "[ERROR] 1" << std::endl;
        std::exit(1);

      }

    } else if (num_daughters == 2) {
      auto d0 = daughters[0];
      auto d1 = daughters[1];

      if (last->PID == d0->PID) {
        last = d0;
        continue;

      } else if (last->PID == d1->PID) {
        last = d1;
        continue;

      } else {
        break;

      }

    } else {
      std::cerr << "[ERROR] 2" << std::endl;
      std::exit(1);
    }

  }

  return last;
}


std::pair<const GenParticle*, const GenParticle*>
TTDileptonAnalyser::getTopDecayProducts(const GenParticle* top_quark) {
  if (std::abs(top_quark->PID) != 6) {
    std::cerr << "[ERROR] 3" << std::endl;
    std::exit(1);
  }

  auto daughters = getDaughters(top_quark);
  if (daughters.size() != 2) {
    std::cerr << "[ERROR] 4" << std::endl;
    std::exit(1);
  }

  auto d0 = daughters[0];
  auto d1 = daughters[1];

  const GenParticle* quark = nullptr;
  const GenParticle* w_boson = nullptr;
  if (((std::abs(d0->PID) == 5) or (std::abs(d0->PID) == 3)) and (std::abs(d1->PID) == 24)) { quark = d0;
    w_boson = d1;

  } else if ((std::abs(d0->PID) == 24) and ((std::abs(d1->PID) == 5) or (std::abs(d1->PID) == 3))) {
    quark = d1;
    w_boson = d0;

  } else {
    auto msg = Form("[ERROR] TTDileptonAnalyser::getTopDecayProducts >>> (PID=%d, Status=%d) --> (%d, %d), (%d, %d)",
        top_quark->PID, top_quark->Status,
        d0->PID, d0->Status,
        d1->PID, d1->Status);
    std::cerr << msg << std::endl;
    std::exit(1);

  }

  return std::make_pair(quark, w_boson);
}


bool TTDileptonAnalyser::isQuark(int pid) {
  TString particle_class{TDatabasePDG::Instance()->GetParticle(pid)->ParticleClass()};
  return particle_class.EqualTo("Quark");
}


std::pair<const GenParticle*, const GenParticle*>
TTDileptonAnalyser::getWBosonDecayProducts(const GenParticle* w_boson) {
  if (std::abs(w_boson->PID) != 24) {
    std::cerr << "[ERROR] 5" << std::endl;
    std::exit(1);
  }

  auto daughters = getDaughters(w_boson);
  if (daughters.size() != 2) {
    std::cerr << "[ERROR] 6" << std::endl;
    std::exit(1);
  }

  const GenParticle* dau0 = daughters[0];
  const GenParticle* dau1 = daughters[1];

  if (std::abs(dau0->PID) == 11 or std::abs(dau0->PID) == 13) return std::make_pair(dau0, dau1);
  else if (std::abs(dau1->PID) == 11 or std::abs(dau1->PID) == 13) return std::make_pair(dau1, dau0);
  else return std::make_pair(dau0, dau1);

//  while (dau0->Status != 23) {
//    dau0 = getParticle(dau0->D1);
//  }
//  while (dau1->Status != 23) {
//    dau1 = getParticle(dau1->D1);
//  }


//  } else {
//    auto msg = Form("[ERROR] TTDileptonAnalyser::getWBosonDecayProducts >>> (PID=%d, Status=%d) --> (%d, %d), (%d, %d)",
//        w_boson->PID, w_boson->Status,
//        dau0->PID, dau0->Status,
//        dau1->PID, dau1->Status);
//    std::cerr << msg << std::endl;
//    std::exit(1);
//  }
}


Bool_t TTDileptonAnalyser::matchPartonsWithJets() {
  // NOTE find top quark and anti top quark
  for (Int_t idx = 0; idx < particles_->GetEntries(); idx++) {
    auto particle = getParticle(idx);
    if (particle->Status > 30) continue;

    if ((particle->PID == 6) and (t1_ == nullptr)) {
      t1_ = particle;
    } else if ((particle->PID == -6) and (t2_ == nullptr)) {
      t2_ = particle;
    }

    if ((t1_ != nullptr) and (t2_ != nullptr)) {
      break;
    }
  }

  if ((t1_ == nullptr) or (t2_ == nullptr)) {
    std::cerr << "[ERROR] top quark not found" << std::endl;
    std::exit(1);
  }

  t1_ = getLast(t1_);
  t2_ = getLast(t2_);

  std::tie(q1_, w1_) = getTopDecayProducts(t1_);
  std::tie(q2_, w2_) = getTopDecayProducts(t2_);

  w1_ = getLast(w1_);
  w2_ = getLast(w2_);

  std::tie(l1_, nu1_) = getWBosonDecayProducts(w1_);
  std::tie(l2_, nu2_) = getWBosonDecayProducts(w2_);

  // std::tie(w1dau0_, w1dau1_) = getWBosonDecayProducts(w1_);
  // std::tie(w2dau0_, w2dau1_) = getWBosonDecayProducts(w2_);

  if ((l1_ == nullptr) or (l2_ == nullptr)) return false;

  // if (std::abs(w1dau0_->PID) > 10 && std::abs(w2dau0_->PID) > 10) {
  //   std::cout << "dilepton channel" << std::endl;
  //   l1_ = w1dau0_;
  //   l2_ = w2dau0_;
  //   if (std::abs(l1_->PID)==std::abs(l2_->PID)) {
  //     if (std::abs(l1_->PID) == 11) {
  //       b_gen_ch_ = 0;
  //     } else {
  //       b_gen_ch_ = 2;
  //     }
  //   } else {
  //     b_gen_ch_ = 1;
  //   }
  // } else {
  //   std::cout << "semilepton channel" << std::endl;
  //   if (std::abs(w1dau0_->PID) > 10) {
  //     l1_ = w1dau0_;
  //     b_gen_ch_ = 3;
  //   } else {
  //     wq1_ = w2dau0_;
  //     wq2_ = w2dau1_;
  //   }
  //   if (std::abs(w2dau0_->PID) > 10) {
  //     l2_ = w2dau0_;
  //     b_gen_ch_ = 4;
  //   } else {
  //     wq1_ = w1dau0_;
  //     wq2_ = w1dau1_;
  //   }
  // }

  // std::vector<const GenParticle*> genLeps({l1_, l2_});

  // partons from W

  // matchingJets
  std::vector<const GenParticle*> partons({q1_, q2_});
  // if ((wq1_ == nullptr) or (wq2_ == nullptr)) {
  //   std::vector<const GenParticle*> partons({q1_, q2_});
  //   std::cout << partons.size() << std::endl;
  // } else {
  //   std::vector<const GenParticle*> partons({q1_, q2_, wq1_, wq2_});
  //   std::cout << partons.size() << std::endl;
  // }
  // NOTE initialize
  for (const auto & parton : partons) {
    parton2jet_[parton] = nullptr;
  }
  for (const auto & jet : selected_jets_) {
    jet2parton_[jet] = nullptr;
  }

  bool is_matched = false;
  is_matched = performTotalMinDist(partons, selected_jets_);

  // sanity check
  if (is_matched) {
    std::set<const Jet*> matched_jets;
    for (const auto& parton_and_jet : parton2jet_) {
      if (matched_jets.find(parton_and_jet.second) == matched_jets.end()) {
        matched_jets.insert(parton_and_jet.second);
      } else {
        auto err_msg = Form("[RuntimeError] (dR_max=%.2f): At least two or more partons share a matched jet",
                            kMaxDist_);
        std::cerr << err_msg << std::endl;
        std::exit(1);
      }
    }
  }
  return is_matched;
}



bool TTDileptonAnalyser::performTotalMinDist(
    const std::vector<const GenParticle*>& partons,
    const std::vector<const Jet*>& jets) {

  // return value
  bool is_matched = true;

  // NOTE https://github.com/cms-sw/cmssw/blob/CMSSW_11_2_0_pre2/TopQuarkAnalysis/TopTools/src/JetPartonMatching.cc#L149-L196
  // match parton to jet with shortest distance
  // starting with the shortest distance available
  // apply some outlier rejection if desired

  // prepare vector of pairs with distances between
  // all partons to all jets in the input vectors
  std::vector<std::pair<double, unsigned int> > distances;
  for (unsigned int ip = 0; ip < partons.size(); ++ip) {
    for (unsigned int ij = 0; ij < jets.size(); ++ij) {
      // double dist = distance(jets[ij]->p4(), partons[ip]->p4());
      Float_t dist = partons[ip]->P4().DeltaR(jets[ij]->P4());
      distances.push_back(std::pair<double, unsigned int>(dist, ip * jets.size() + ij));
    }
  }
  std::sort(distances.begin(), distances.end());

  std::vector<std::pair<unsigned int, int> > match;

  while (match.size() < partons.size()) {
    unsigned int partonIndex = distances[0].second / jets.size();
    int jetIndex = distances[0].second - jets.size() * partonIndex;

    // use primitive outlier rejection if desired
    if (kUseMaxDist_ && distances[0].first > kMaxDist_)
      jetIndex = -1;

    // prevent underflow in case of too few jets
    if (distances.empty())
      match.push_back(std::make_pair(partonIndex, -1));
    else
      match.push_back(std::make_pair(partonIndex, jetIndex));

    // remove all values for the matched parton
    // and the matched jet
    for (unsigned int a = 0; a < distances.size(); ++a) {
      unsigned int pIndex = distances[a].second / jets.size();
      int jIndex = distances[a].second - jets.size() * pIndex;
      if ((pIndex == partonIndex) || (jIndex == jetIndex)) {
        distances.erase(distances.begin() + a, distances.begin() + a + 1);
        --a;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////
  // NOTE
  //////////////////////////////////////////////////////////////////////////////
  for (auto&& [parton_index, jet_index] : match) {
    if (jet_index < 0) {
      is_matched = false;
      break;
    }

    const GenParticle* parton = partons[parton_index];
    const Jet* jet = jets[jet_index];

    parton2jet_[parton] = jet;
    jet2parton_[jet] = parton;
  }


  if (is_matched) {
    for (const auto& parton_and_jet : parton2jet_) {
      if (parton_and_jet.second == nullptr) {
        is_matched = false;
        break;
      }
    }
  }
  return is_matched;
}



Bool_t TTDileptonAnalyser::selectEvent() {
  // initial
  cutflow_->Fill(1);

  for (Int_t i = 0; i < events_->GetEntries(); i++) {
    auto event = (HepMCEvent*) events_->At(i);
    if (event->Weight == 0) b_weight_ = 0;
    else b_weight_ = event->Weight > 0 ? 1 : -1;
  }
  gen_weight_->Fill(0.5, b_weight_);

  // lepton selection
  for (Int_t i = 0; i < muons_->GetEntries(); i++) {
    auto muon = dynamic_cast<const Muon*>(muons_->At(i));
    if (muon->PT < kMinLeptonPT_) continue;
    if (std::fabs(muon->Eta) > kMaxLeptonEta_) continue;
    if (muon->IsolationVar > kMuonIso_) continue;
    //if (muon->P4().M() < 0.0f) continue;
    selected_muons_.push_back(muon);
  }
  UInt_t num_selected_muons  = selected_muons_.size();

  for (Int_t i = 0; i < electrons_->GetEntries(); i++) {
    auto electron = dynamic_cast<const Electron*>(electrons_->At(i));
    if (electron->PT < kMinLeptonPT_) continue;
    if (std::fabs(electron->Eta) > kMaxLeptonEta_) continue;
    if (electron->IsolationVar > kElectronIso_) continue;
    if (abs(electron->Eta) < cut_lowElecEffRegionUp && abs(electron->Eta) > cut_lowElecEffRegionDn) continue;
    if (abs(electron->Eta) < cut_ElectronEtaBarrel && electron->IsolationVar > cut_ElectronRelIso03All_inBarrel) continue;
    if (abs(electron->Eta) > cut_lowElecEffRegionUp && electron->IsolationVar > cut_ElectronRelIso03All_inEndCap) continue;
    //if (electron->P4().M() < 0.0f) continue;
    selected_electrons_.push_back(electron);
  }
  UInt_t num_selected_electrons  = selected_electrons_.size();

  // Di-lepton
  if (num_selected_muons + num_selected_electrons != 2) return false;
  // b_ch / 0: ee / 1: emu / 2: mumu
  int charge;

  if (num_selected_electrons == 2) {
    b_ch_ = 0;
    lep_p4_.push_back(selected_electrons_[0]->P4());
    lep_p4_.push_back(selected_electrons_[1]->P4());
    lep_charge_.push_back(selected_electrons_[0]->Charge);
    lep_charge_.push_back(selected_electrons_[1]->Charge);
    lep_pid_.push_back(-11*selected_electrons_[0]->Charge);
    lep_pid_.push_back(-11*selected_electrons_[1]->Charge);
    charge = selected_electrons_[0]->Charge*selected_electrons_[1]->Charge;
  } else if (num_selected_electrons == 1) {
    b_ch_ = 1;
    if (selected_muons_[0]->PT > selected_electrons_[0]->PT) {
      lep_p4_.push_back(selected_muons_[0]->P4());
      lep_p4_.push_back(selected_electrons_[0]->P4());
      lep_charge_.push_back(selected_muons_[0]->Charge);
      lep_charge_.push_back(selected_electrons_[0]->Charge);
      lep_pid_.push_back(-13*selected_muons_[0]->Charge);
      lep_pid_.push_back(-11*selected_electrons_[0]->Charge);
    } else {
      lep_p4_.push_back(selected_electrons_[0]->P4());
      lep_p4_.push_back(selected_muons_[0]->P4());
      lep_charge_.push_back(selected_electrons_[0]->Charge);
      lep_charge_.push_back(selected_muons_[0]->Charge);
      lep_pid_.push_back(-11*selected_electrons_[0]->Charge);
      lep_pid_.push_back(-13*selected_muons_[0]->Charge);
    }
    charge = selected_electrons_[0]->Charge*selected_muons_[0]->Charge;
  } else if (num_selected_electrons == 0) {
    b_ch_ = 2;
    lep_p4_.push_back(selected_muons_[0]->P4());
    lep_p4_.push_back(selected_muons_[1]->P4());
    lep_charge_.push_back(selected_muons_[0]->Charge);
    lep_charge_.push_back(selected_muons_[1]->Charge);
    lep_pid_.push_back(-13*selected_muons_[0]->Charge);
    lep_pid_.push_back(-13*selected_muons_[1]->Charge);
    charge = selected_muons_[0]->Charge*selected_muons_[1]->Charge;
  }

  auto dilepton = lep_p4_[0] + lep_p4_[1];

  if ( (dilepton.M() < kMinDileptonMass_) or (charge > 0) ) return false;
  if ( lep_p4_[0].Pt() < kMinFirstLeptonPT_  &&  lep_p4_[1].Pt() < kMinFirstLeptonPT_ )return false;
  cutflow_->Fill(2);

  // veto Z boson
  //if ( (b_ch_ != 1) && (std::abs(kZbosonMass_ - dilepton.M()) < 15 ) ) return false;
  if ( (b_ch_ != 1) && ( dilepton.M() > 76. && dilepton.M() < 106. ) ) return false;
  cutflow_->Fill(3);

  // MET cut
  auto missing_et = getMissingET(0);
  b_met_     = missing_et->MET;
  b_met_eta_ = missing_et->Eta;
  b_met_phi_ = missing_et->Phi;
  if ( b_ch_ != 1 && b_met_ < kMinMET_ ) return false;
  cutflow_->Fill(4);

  // jet selection
  for (Int_t i = 0; i < jets_->GetEntries(); i++) {
    auto jet = dynamic_cast<const Jet*>(jets_->At(i));
    if (jet->PT < kMinJetPT_) continue;
    if (std::fabs(jet->Eta) > kMaxJetEta_) continue;
    //if (jet->Mass < 0.0f) continue;
    if (jet->P4().DeltaR(lep_p4_[0]) < kJetOverlap_) continue;
    if (jet->P4().DeltaR(lep_p4_[1]) < kJetOverlap_) continue;
    selected_jets_.push_back(jet);
  }

  UInt_t  num_selected_jets  = selected_jets_.size();
  if (num_selected_jets < kHLTMinNumJets_) return false;
  cutflow_->Fill(5);

  // bjet selection
  std::vector<const Jet*> selected_b_jets;
  std::copy_if(selected_jets_.begin(),
               selected_jets_.begin() + 2,
               std::back_inserter(selected_b_jets),
               [](const Jet* each) {return each->BTag == 1;});

  UInt_t  num_b_jets = selected_b_jets.size();
  if (num_b_jets > kHLTMinNumBJets_) return false;
  cutflow_->Fill(6);

  //////////////////////////////////////////////////////////////////////////////
  // NOTE Jet-Parton Matching
  //////////////////////////////////////////////////////////////////////////////
  if (kIsTT_) {
    b_is_matched_ = matchPartonsWithJets();
    if (b_is_matched_)
      cutflow_->Fill(7);
  }
  b_is_tt_ = kIsTT_;

  // NOTE return HLT;
  return true;
}


Int_t TTDileptonAnalyser::getPFIndex(Int_t pid, Int_t charge) {
  if      (pid == pdgid::kElectron)     return pfindex::kElectron;
  else if (pid == pdgid::kAntiElectron) return pfindex::kAntiElectron;
  else if (pid == pdgid::kMuon)         return pfindex::kMuon;
  else if (pid == pdgid::kAntiMuon)     return pfindex::kAntiMuon;
  else if (charge > 0)                  return pfindex::kPositivelyChargedHadron;
  else                                  return pfindex::kNegativelyChargedHadron;
}


std::pair<Int_t, Int_t>
TTDileptonAnalyser::getJetPartonMatchingCode(const GenParticle* parton) {
  Int_t code1, code2;
  if (parton == q1_) {
    code1 = parton->PID;
  } else if (parton == q2_) {
    code1 = parton->PID;     // t or tbar or other
  } else {
    code1 = -1;
  }

  if (std::abs(parton->PID) == 5) {
    code2 = 1;
  } else if (std::abs(parton->PID) == 3) {
    code2 = 2;
  } else {
    code2 = -1;
  }
  return std::make_pair(code1, code2);
}

void TTDileptonAnalyser::HadronReconstruction() {
  for (Int_t i = 0; i < tracks_->GetEntries(); i++) {
    auto hadCand1 = (Track*) tracks_->At(i);
    if (hadCand1->Charge != 1) continue;
    if (abs(hadCand1->PID) == 11 || abs(hadCand1->PID) == 13) continue;
    if (hadCand1->PT < tkPTCut_) continue;
    if (fabs(hadCand1->Eta) > tkEtaCut_) continue;
    if (fabs(hadCand1->D0/hadCand1->ErrorD0) < tkIPSigXYCut_) continue;
    for (Int_t j = 0; j < tracks_->GetEntries(); j++) {
      auto hadCand2 = (Track*) tracks_->At(j);
      if (hadCand2->Charge != 1) continue;
      if (abs(hadCand2->PID) == 11 || abs(hadCand2->PID) == 13) continue;
      if (hadCand1->PT < tkPTCut_) continue;
      if (fabs(hadCand2->Eta) > tkEtaCut_) continue;
      if (fabs(hadCand2->D0/hadCand2->ErrorD0) < tkIPSigXYCut_) continue;

      TLorentzVector dauCand1_pion_tlv;
      TLorentzVector dauCand1_proton_tlv;
      TLorentzVector dauCand2_pion_tlv;
      TLorentzVector dauCand2_proton_tlv;

      auto genDau1 = (GenParticle*) hadCand1->Particle.GetObject();
      auto genDau2 = (GenParticle*) hadCand2->Particle.GetObject();
      dauCand1_pion_tlv.SetPtEtaPhiM(hadCand1->PT, hadCand1->Eta, genDau1->Phi, pionMass_);
      dauCand1_proton_tlv.SetPtEtaPhiM(hadCand1->PT, hadCand1->Eta, genDau1->Phi, protonMass_);
      dauCand2_pion_tlv.SetPtEtaPhiM(hadCand2->PT, hadCand2->Eta, genDau2->Phi, pionMass_);
      dauCand2_proton_tlv.SetPtEtaPhiM(hadCand2->PT, hadCand2->Eta, genDau2->Phi, protonMass_);

      TLorentzVector hadCand_pion_pion_tlv = dauCand1_pion_tlv + dauCand2_pion_tlv;
      TLorentzVector hadCand_pion_proton_tlv = dauCand1_pion_tlv + dauCand2_proton_tlv;
      TLorentzVector hadCand_proton_pion_tlv = dauCand1_proton_tlv + dauCand2_pion_tlv;

      if ((fabs(hadCand_pion_pion_tlv.M() - ksMass_)/ksMass_) < 0.1) { // Invariant mass cut
        m_recoHad.push_back({hadCand_pion_pion_tlv, dauCand1_pion_tlv, dauCand2_pion_tlv, -1, i, j, ksPdgId_, pionPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_pion_proton_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) { // Invariant mass cut
        m_recoHad.push_back({hadCand_pion_proton_tlv, dauCand1_pion_tlv, dauCand2_proton_tlv, -1, i, j, lambda0PdgId_, pionPdgId_, (-1)*protonPdgId_, -99, false, false});
      }
      if ((fabs(hadCand_proton_pion_tlv.M() - lambda0Mass_)/lambda0Mass_) < 0.3) { // Invariant mass cut
        m_recoHad.push_back({hadCand_proton_pion_tlv, dauCand1_proton_tlv, dauCand2_pion_tlv, -1, i, j, lambda0PdgId_, protonPdgId_, (-1)*pionPdgId_, -99, false, false});
      }
    }
  }
}




void TTDileptonAnalyser::analyseJet() {
  Int_t jet_count = 0;
  b_num_jets_ = selected_jets_.size();

  b_num_b_jets_ = 0;
  b_num_b_jets_2_ = 0;

  for (const auto & jet : selected_jets_) {
    jet_count++;

    TLorentzVector jet_p4 = jet->P4();
    float j_pt = jet->PT;
    float j_eta = jet->Eta;
    float j_phi = jet->Phi;
    //float j_e = jet->P4().E();   #### fix
    float j_e = jet_p4.E();
    float j_charge = 0;

    if (jet->BTag == 1) {
      b_num_b_jets_++;

      if(jet_count <= 2) b_num_b_jets_2_++;
    }

    b_jet_.push_back(jet->P4());
    b_jet_pt_.push_back(j_pt);
    b_jet_eta_.push_back(j_eta);
    b_jet_phi_.push_back(j_phi);
    b_jet_mass_.push_back(jet->Mass);
    b_jet_energy_.push_back(j_e);

    float jet_lep1_dr = jet_p4.DeltaR(lep_p4_[0]);
    float jet_lep2_dr = jet_p4.DeltaR(lep_p4_[1]);
    b_jet_lep1_dr_.push_back(jet_lep1_dr);
    b_jet_lep2_dr_.push_back(jet_lep2_dr);

    // Loop over cons

    std::vector<float> con_pt;
    std::vector<float> con_deta;
    std::vector<float> con_dphi;

    std::vector<float> track_pt;
    std::vector<float> track_eta;
    std::vector<float> track_phi;
    std::vector<float> track_deta;
    std::vector<float> track_dphi;
    std::vector<float> track_x;
    std::vector<float> track_ptrel;
    std::vector<float> track_pzrel;
    std::vector<float> track_d0;
    std::vector<float> track_dz;
    std::vector<float> track_delta_x;
    std::vector<float> track_delta_y;
    std::vector<float> track_delta_z;
    std::vector<float> track_pid;

    std::vector<float> tower_pt;
    std::vector<float> tower_eta;
    std::vector<float> tower_deta;
    std::vector<float> tower_dphi;
    std::vector<float> tower_x;
    std::vector<float> tower_ptrel;
    std::vector<float> tower_pzrel;
    std::vector<float> tower_d0;
    std::vector<float> tower_dz;
    std::vector<float> tower_E;
    std::vector<float> tower_Eem;
    std::vector<float> tower_Ehad;


    Int_t num_chad = 0, num_nhad = 0;
    Int_t num_electron = 0, num_muon = 0, num_photon = 0;

    int num_constituents = jet->Constituents.GetEntries();
    for (Int_t con_idx = 0; con_idx < num_constituents; con_idx++) {
      auto con = jet->Constituents.At(con_idx);

      float c_pt = -1;
      float c_eta = 100;
      float c_phi = 100;
      float c_charge = 100;
      float c_e = -1;

      if (auto track = dynamic_cast<const Track*>(con)) {
        c_pt = track->PT;
        c_eta = track->Eta;
        c_phi = track->Phi;
        c_charge = track->Charge;
        c_e = track->P4().E();

        Int_t abs_con_pid = std::abs(track->PID);
        if      (abs_con_pid == pdgid::kElectron)  num_electron++;
        else if (abs_con_pid == pdgid::kMuon)      num_muon++;
        else                                       num_chad++;

        if (kAlgorithm_ == "MomentumWeighted") {
          j_charge += c_charge * std::pow(c_pt/j_pt, kExponent_);
        } else if (kAlgorithm_ == "EnergyWeighted") {
          j_charge += c_charge * std::pow(c_e/j_e, kExponent_);
        }
        TLorentzVector trk_p4 = track->P4();
        float trk_pt = track->PT;
        float trk_eta = track->Eta;
        float trk_phi = track->Phi;
        float trk_d0 = track->D0;
        float trk_dz = track->DZ;
        float trk_deta = trk_eta - j_eta;
        float trk_dphi = TVector2::Phi_mpi_pi(trk_phi - j_phi);
        float trk_x = trk_pt / j_pt;
        float trk_pzrel = trk_p4.Dot(jet_p4) / jet_p4.P();
        float trk_ptrel = std::pow(trk_p4.P()*trk_p4.P() - trk_pzrel*trk_pzrel, 0.5);
        float trk_deltax = track->Xd;
        float trk_deltay = track->Yd;
        float trk_deltaz = track->Zd;
        float trk_pid = track->PID;

        track_pt.push_back(trk_pt);
        track_eta.push_back(trk_eta);
        track_phi.push_back(trk_phi);
        track_deta.push_back(trk_deta);
        track_dphi.push_back(trk_dphi);
        track_x.push_back(trk_x);
        track_ptrel.push_back(trk_ptrel);
        track_pzrel.push_back(trk_pzrel);
        track_d0.push_back(trk_d0);
        track_dz.push_back(trk_dz);
        track_delta_x.push_back(trk_deltax);
        track_delta_y.push_back(trk_deltay);
        track_delta_z.push_back(trk_deltaz);
        track_pid.push_back(trk_pid);
      } else if (auto tower = dynamic_cast<const Tower*>(con)) {
        c_pt = tower->ET;
        c_eta = tower->Eta;
        c_phi = tower->Phi;

        if (c_pt < 1) continue;
        if (tower->Eem == 0.0) {
          num_nhad++;
        } else if (tower->Ehad == 0.0) {
          num_photon++;
        } else {
          std::cout << "[[ERROR]]: Tower with Had " << tower->Ehad
                    << " and EM " << tower->Eem << " energy" << std::endl;
          }
        TLorentzVector tow_p4 = tower->P4();
        float tow_pt = tower->ET;
        float tow_eta = tower->Eta;
        float tow_phi = tower->Phi;
        float tow_E = tower->E;
        float tow_Eem = tower->Eem;
        float tow_Ehad = tower->Ehad;
        float tow_deta = tow_eta - j_eta;
        float tow_dphi = TVector2::Phi_mpi_pi(tow_phi - j_phi);
        float tow_x = tow_pt / j_pt;
        float tow_pzrel = tow_p4.Dot(jet_p4) / jet_p4.P();
        float tow_ptrel = std::sqrt(tow_p4.P()*tow_p4.P() - tow_pzrel*tow_pzrel);
        // float tow_pid = tower->Particles->PID;

        tower_pt.push_back(tow_pt);
        tower_eta.push_back(tow_eta);
        tower_deta.push_back(tow_deta);
        tower_dphi.push_back(tow_dphi);
        tower_x.push_back(tow_x);
        tower_ptrel.push_back(tow_ptrel);
        tower_pzrel.push_back(tow_pzrel);
        tower_E.push_back(tow_E);
        tower_Eem.push_back(tow_Eem);
        tower_Ehad.push_back(tow_Ehad);
      } else {
        std::cout << "[[ERROR]]: BAD DAUGHTER! " << con << std::endl;
      }

      float c_deta = c_eta - j_eta;
      float c_dphi = TVector2::Phi_mpi_pi(c_phi - j_phi);
      con_pt.push_back(c_pt);
      con_deta.push_back(c_deta);
      con_dphi.push_back(c_dphi);
    } // end loop over cons

    b_track_pt_.push_back(track_pt);
    b_track_eta_.push_back(track_eta);
    b_track_phi_.push_back(track_phi);
    b_track_deta_.push_back(track_deta);
    b_track_dphi_.push_back(track_dphi);
    b_track_x_.push_back(track_x);
    b_track_ptrel_.push_back(track_ptrel);
    b_track_pzrel_.push_back(track_pzrel);
    b_track_d0_.push_back(track_d0);
    b_track_dz_.push_back(track_dz);
    b_track_delta_x_.push_back(track_delta_x);
    b_track_delta_y_.push_back(track_delta_y);
    b_track_delta_z_.push_back(track_delta_z);
    b_track_pid_.push_back(track_pid);

    b_tower_pt_.push_back(track_pt);
    b_tower_eta_.push_back(track_eta);
    b_tower_deta_.push_back(track_deta);
    b_tower_dphi_.push_back(track_dphi);
    b_tower_x_.push_back(track_x);
    b_tower_ptrel_.push_back(track_ptrel);
    b_tower_pzrel_.push_back(track_pzrel);
    b_tower_E_.push_back(tower_E);
    b_tower_Eem_.push_back(tower_Eem);
    b_tower_Ehad_.push_back(tower_Ehad);


    ////////////////////////////////////////////////////////////////////////////
    // NOTE CMS Variables
    ////////////////////////////////////////////////////////////////////////////
    b_jet_num_chad_.push_back(num_chad);
    b_jet_num_nhad_.push_back(num_nhad);
    b_jet_num_electron_.push_back(num_electron);
    b_jet_num_muon_.push_back(num_muon);
    b_jet_num_photon_.push_back(num_photon);
    // Move CMSSW to higher version and use structured bindings
    float major_axis, minor_axis;
    std::tie(major_axis, minor_axis) = computeAngularSpread(
        con_deta, con_dphi, con_pt);

    b_jet_major_axis_.push_back(major_axis);
    b_jet_minor_axis_.push_back(minor_axis);

    // jet energy sharing variable
    float ptd = computeJetEnergySharingVariable(con_pt);
    b_jet_ptd_.push_back(ptd);

    b_jet_charge_.push_back(j_charge);

    ////////////////////////////////////////////////////////////////////////////
    // NOTE b-jet properties
    ////////////////////////////////////////////////////////////////////////////
    b_jet_b_tag_.push_back(jet->BTag);
    b_jet_flavor_.push_back(jet->Flavor);
    ////////////////////////////////////////////////////////////////////////////
    // NOTE jet-parton matching
    ////////////////////////////////////////////////////////////////////////////

    if (kIsTT_) {
      if (b_is_matched_) {
        Int_t matching_code = 999;
        Int_t detail_code = 999;

        Int_t parton_pid = 999; // detail
        Int_t is_t_decay = 999;
        Int_t is_b_jet = 999;

        if (auto parton = jet2parton_[jet]) {
          std::tie(matching_code, detail_code) = getJetPartonMatchingCode(parton);

          parton_pid = parton->PID;
          is_t_decay = (matching_code == 1) ? 1 : -1;
          is_b_jet = (detail_code == 1) ? 1 : -1;
        } else {
          matching_code = 0;
          detail_code = 0;

          parton_pid = 0;
          is_t_decay = 0;
          is_b_jet = 0;
        }

        b_jet_parton_match_.push_back(matching_code);
        b_jet_parton_match_detail_.push_back(detail_code);
        b_jet_parton_match_pid_.push_back(parton_pid);
        b_jet_parton_match_t_.push_back(is_t_decay);
        b_jet_parton_match_b_.push_back(is_b_jet);
      } else {
        b_jet_parton_match_.push_back(0);

        Int_t matching_code = 999;
        Int_t detail_code = 999;
        if (auto parton = jet2parton_[jet]) {
          std::tie(matching_code, detail_code) = getJetPartonMatchingCode(parton);
        } else {
          detail_code = 0;
        }
        b_jet_parton_match_detail_.push_back(detail_code);
      }
    }
  } // end loop over jets

}


void TTDileptonAnalyser::analyseLepton() {

  //b_lep_pid_.push_back(lep1_pid_);
  //b_lep_pid_.push_back(lep2_pid_);
  //b_lep_charge_.push_back(lep1_charge_);
  //b_lep_charge_.push_back(lep2_charge_);
  //for (const auto & lep : lep_p4_) {
  //  b_lep_.push_back(*lep);
  //  b_lep_pt_.push_back(lep->Pt());
  //  b_lep_eta_.push_back(lep->Eta());
  //  b_lep_phi_.push_back(lep->Phi());
  //  b_lep_mass_.push_back(lep->M());
  //  b_lep_energy_.push_back(lep->E());
  //} // end loop over leptons

  for (Int_t i = 0; i < 2; i++) {
    auto lep = lep_p4_[i];
    b_lep_.push_back(lep);
    b_lep_pt_.push_back(lep.Pt());
    b_lep_eta_.push_back(lep.Eta());
    b_lep_phi_.push_back(lep.Phi());
    b_lep_mass_.push_back(lep.M());
    b_lep_energy_.push_back(lep.E());

    b_lep_pid_.push_back(lep_pid_[i]);
    b_lep_charge_.push_back(lep_charge_[i]);
    b_lep_isMuon_.push_back((std::abs(lep_pid_[i])==13));
    //b_jet_parton_match_.push_back(lep_charge_[i] == 1? 3: 4);
  }
}


void TTDileptonAnalyser::analyse() {
  analyseLepton();
  analyseJet();
  // for (Int_t i = 0; i < events_->GetEntries(); i++) {
  //   auto event = (HepMCEvent*) events_->At(i);
  //   if (event->Weight == 0) b_weight_ = 0;
  //   else b_weight_ = event->Weight > 0 ? 1 : -1;
  // }

  if (kIsTT_) {
    if (b_is_matched_) {
      // NOTE partons
      b_t1_ = t1_->P4();
      b_q1_ = q1_->P4();
      b_w1_ = w1_->P4();
      b_l1_ = l1_->P4();
      b_nu1_ = nu1_->P4();

      b_t2_ = t2_->P4();
      b_q2_ = q2_->P4();
      b_w2_ = w2_->P4();
      b_l2_ = l2_->P4();
      b_nu2_ = nu2_->P4();

      b_l1_pid_  = l1_->PID;
      b_nu1_pid_  = nu1_->PID;
      b_l2_pid_  = l2_->PID;
      b_nu2_pid_  = nu2_->PID;

      auto nu1_p4 = nu1_->P4();
      auto nu2_p4 = nu2_->P4();
      b_nu1_pt_.push_back(nu1_p4.Pt());
      b_nu1_eta_.push_back(nu1_p4.Eta());
      b_nu1_phi_.push_back(nu1_p4.Phi());
      b_nu2_pt_.push_back(nu2_p4.Pt());
      b_nu2_eta_.push_back(nu2_p4.Eta());
      b_nu2_phi_.push_back(nu2_p4.Phi());

      for (const auto& [parton, jet] : parton2jet_) {
        float delta_energy = std::fabs(parton->E - jet->P4().Energy());
        float delta_pt = std::fabs(parton->PT - jet->PT);
        float delta_r = parton->P4().DeltaR(jet->P4());

        b_parton_jet_distance_.push_back(delta_r);

        b_sum_delta_energy_ += delta_energy;
        b_sum_delta_pt_ += delta_pt;
        b_sum_delta_r_ += delta_r;

      }

    }
  }

  if (kIsTT_) {
    if (b_is_matched_) {
      b_label_ = 1;
      out_tree_->Fill();
    } else {
      b_label_ = 2;
      ambig_tree_->Fill();
      //out_tree_->Fill();
    }
  } else {
    b_label_ = 0;
    out_tree_->Fill();
  }

}


void TTDileptonAnalyser::loop() {
  const Int_t kNumTotal = in_tree_->GetEntries();
  const Int_t kPrintFreq = kNumTotal / 20;
  TString msg_fmt = TString::Format("[%s/%d (%s %%)]", "%d", kNumTotal, "%.0f");

  for (Long64_t entry=0; entry < in_tree_->GetEntries(); entry++) {
    in_tree_->GetEntry(entry);

    if (entry == 0) std::cout << "start resetBranch" << std::endl;
    resetBranch();

    if (entry == 0) std::cout << "start resetMemberData" << std::endl;
    resetMemberData();

    if (entry == 0) std::cout << "start selectEvent" << std::endl;

    if (selectEvent()) {
      if (entry == 0) std::cout << "start analyse" << std::endl;
      analyse();
    } else {
      if (entry == 0) std::cout << "end selectEvent" << std::endl;
    }

    if (entry % kPrintFreq == 0) {
      float progress = 100 * float(entry) / float(kNumTotal);
      std::cout << TString::Format(msg_fmt, entry, progress) << std::endl;
    }
  }
}
