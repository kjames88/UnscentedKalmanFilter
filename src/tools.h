#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  static Eigen::VectorXd PolarToCartesian(Eigen::VectorXd const& polar);
  static double normalize_angle(double a) {
    double norm = a;
    while (norm > M_PI)
      norm -= (2.0 * M_PI);
    while (norm < -M_PI)
      norm += (2.0 * M_PI);
    return norm;
  }
  
};

#endif /* TOOLS_H_ */
