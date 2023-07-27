#ifndef DELPHES_EXTERNAL_JETDISCRIMINATINGVARIABLE_H_
#define DELPHES_EXTERNAL_JETDISCRIMINATINGVARIABLE_H_

#include <tuple>
#include <vector>

namespace delphes {

std::tuple<float, float> computeAngularSpread(
    const std::vector<float> & constituents_deta,
    const std::vector<float> & constituents_dphi,
    const std::vector<float> & constituents_pt);

std::tuple<float, float> computeAngularSpread(
    float jet_eta, float jet_phi,
    const std::vector<float> & constituents_eta,
    const std::vector<float> & constituents_phi,
    const std::vector<float> & constituents_pt);

float computeJetEnergySharingVariable(
    const std::vector<float> & constituents_pt);

} // NOTE end of delphys namespace

#endif // DELPHES_EXTERNAL_JETDISCRIMINATINGVARIABLE_H_
