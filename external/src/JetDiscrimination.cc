#include "delphes/external/interface/JetDiscrimination.h"

#include "TMatrixTSym.h"
#include "TMatrixDSymfwd.h"
#include "TMatrixFSymfwd.h"
#include "TMatrixDfwd.h"
#include "TMatrixTBase.h" // TVectorT
#include "TVectorF.h"
#include "TVectorD.h"
#include "TVector2.h"

#include <iostream>
#include <numeric> // inner_product
#include <algorithm> // transform
#include <cmath>
#include <vector>


namespace delphes {

std::tuple<float, float> computeAngularSpread(
    const std::vector<float> & constituents_delta_eta,
    const std::vector<float> & constituents_delta_phi,
    const std::vector<float> & constituents_pt) {

  if (constituents_delta_eta.size() != constituents_delta_phi.size() or
      constituents_delta_phi.size() != constituents_pt.size() or
      constituents_pt.size() != constituents_delta_eta.size()) {

    std::cerr << "computeAngularSpread:: different size" << std::endl;
    return std::make_tuple(-1.0, -1.0);
  }

  Int_t multiplicity = constituents_delta_eta.size();

  // sum of weight weight 
  float w2_sum = 0.0;
  // components of a covariance matrix
  float m00 = 0.0, m01 = 0.0, m11 = 0.0;

  float x, y, w2;
  for (Int_t idx = 0; idx < multiplicity; idx++) {
    x = constituents_delta_eta.at(idx);
    y = constituents_delta_phi.at(idx);
    w2 = std::pow(constituents_pt.at(idx), 2);

    w2_sum += w2;

    m00 += w2 * std::pow(x, 2);
    m01 -= w2 * x * y;
    m11 += w2 * std::pow(y, 2);
  } 

  // TMatrixTSm
  // Note that in this implementation both matrix element m[i][j] and m[j][i]
  // are updated and stored in memory . However, when making the object
  // persistent only the upper right triangle is stored .
  TMatrixTSym<float> covariance_matrix(2);
  covariance_matrix(0, 0) = m00;
  covariance_matrix(0, 1) = m01;
  covariance_matrix(1, 1) = m11;

  TVectorT<float> eigen_values;
  covariance_matrix.EigenVectors(eigen_values);

  // length of
  float major_axis = std::sqrt(eigen_values[0] / w2_sum);
  float minor_axis = std::sqrt(eigen_values[1] / w2_sum);

  //float eccentricity = std::sqrt(1 - (eigen_values[1] / eigen_values[0]));
  // return std::make_tuple(major_axis, minor_axis, eccentricity);
  return std::make_tuple(major_axis, minor_axis);
}

std::tuple<float, float> computeAngularSpread(
    float jet_eta, float jet_phi,
    const std::vector<float> & constituents_eta,
    const std::vector<float> & constituents_phi,
    const std::vector<float> & constituents_pt) {

  std::vector<float> constituents_delta_eta;
  std::transform(constituents_eta.begin(), constituents_eta.end(),
                 std::back_inserter(constituents_delta_eta),
                 [jet_eta] (float eta) -> float {
                   return eta - jet_eta;});

  std::vector<float> constituents_delta_phi;
  std::transform(constituents_phi.begin(), constituents_phi.end(),
                 std::back_inserter(constituents_delta_phi),
                 [jet_phi] (float phi) -> float {
                   return phi - jet_phi;});

  return computeAngularSpread(constituents_delta_eta, constituents_delta_phi,
                              constituents_pt);
}

float computeJetEnergySharingVariable(
    const std::vector<float> & constituents_pt) {

  float l1_norm = std::accumulate(
      constituents_pt.begin(), constituents_pt.end(), 0.0f);

  float l2_norm_squared = std::inner_product(
      constituents_pt.begin(), constituents_pt.end(),
      constituents_pt.begin(), 0.0f);

  float ptd = std::sqrt(l2_norm_squared) / l1_norm;

  return ptd;
}

} // end delphes namespace
