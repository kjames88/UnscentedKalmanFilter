#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse = VectorXd(4);  // 4d state
  rmse << 0, 0, 0, 0;
  if (estimations.size() > 0 &&
      estimations.size() == ground_truth.size()) {
    VectorXd sum = VectorXd(4);
    sum << 0, 0, 0, 0;
    for (int i=0; i < estimations.size(); ++i) {
      VectorXd d = estimations[i] - ground_truth[i];
      VectorXd d_square = d.array() * d.array();
      sum += d_square;
    }
    VectorXd mean = sum / estimations.size();
    rmse = mean.array().sqrt();
  }
  return rmse;
}

VectorXd Tools::PolarToCartesian(VectorXd const& polar) {
  double rho = polar[0];
  double phi = polar[1];
  double rho_dot = polar[2];
  double px = rho * cos(phi);
  double py = rho * sin(phi);
  VectorXd cartesian = VectorXd(2);
  cartesian << px, py;
  return cartesian;
}
