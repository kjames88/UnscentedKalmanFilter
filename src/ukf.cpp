#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <assert.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.97;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  is_initialized_ = false;
  time_us_ = 0;
  NIS_radar_ = 0;
  NIS_laser_ = 0;
  n_x_ = 5;  // {px, py, v, psi, psi_dot}
  n_aug_ = 7; // augment state with noise {vx, vy}
  lambda_ = 3 - n_aug_;
  n_sigma_ = 2 * n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, n_sigma_);
  weights_ = VectorXd(n_sigma_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i < n_sigma_; i++) {
    weights_(i) = 0.5 / (lambda_ + n_aug_);
  }

  std::cout << "weights" << std::endl << weights_ << std::endl;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  
  if (is_initialized_ == false) {
    VectorXd meas_x_std(n_x_);
    meas_x_std.fill(0.0);
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      VectorXd polar = meas_package.raw_measurements_;
      VectorXd cartesian = Tools::PolarToCartesian(polar);
      meas_x_std << cartesian[0], cartesian[1], 0, 0, 0;
    } else {
      meas_x_std << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
    }
    x_ = meas_x_std;
    P_ << (std_laspx_ * std_laspx_), 0, 0, 0, 0,
      0, (std_laspy_ * std_laspy_), 0, 0, 0,
      0, 0, 1.0, 0, 0,
      0, 0, 0, (std_radphi_ * std_radphi_), 0,
      0, 0, 0, 0, (std_radrd_ * std_radrd_);
    timestamp_q_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }
  
  double delta_t = ((double) meas_package.timestamp_ - timestamp_q_) / 1000000.0;  // usec to sec
  timestamp_q_ = meas_package.timestamp_;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    VectorXd polar = meas_package.raw_measurements_;
    VectorXd cartesian = Tools::PolarToCartesian(polar);
    std::cout << "measurements radar" << std::endl << cartesian << std::endl;
    
    if (use_radar_) {
      Prediction(delta_t);
      UpdateRadar(meas_package);
    }
  } else {

    std::cout << "measurements laser" << std::endl << meas_package.raw_measurements_ << std::endl;

    if (use_laser_) {
      Prediction(delta_t);
      UpdateLidar(meas_package);
    }
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  // Generate augmented covariance matrix Pa
  MatrixXd Q(2,2);
  Q << (std_a_ * std_a_), 0, 0, (std_yawdd_ * std_yawdd_);
  MatrixXd Pa(n_aug_, n_aug_);
  Pa.fill(0.0);
  Pa.topLeftCorner(n_x_, n_x_) = P_;
  Pa.block<2,2>(n_x_, n_x_) = Q;

  std::cout << "Pa" << std::endl << Pa << std::endl;
  
  // Predict sigma points
  //   Generate sigma points
  MatrixXd Xsig(n_aug_, n_sigma_);
  VectorXd x_aug(n_aug_);
  x_aug << x_, 0, 0;
  Xsig.col(0) = x_aug;
  double sq = sqrt(lambda_ + n_aug_);  // sqrt(3)
  MatrixXd A = Pa.llt().matrixL();  // matrix sqrt
  
  for (int i=0; i<n_aug_; i++) {
    VectorXd m = sq * A.col(i);
    Xsig.col(1 + i) = x_aug + m;
    Xsig.col(1 + n_aug_ + i) = x_aug - m;
  }

  std::cout << "pred Xsig " << std::endl << Xsig << std::endl;
  
  //  Map sigma points to create prediction X_k+1|k
  //    Equations located at:  Sigma Point Prediction Assignment 1
  Xsig_pred_.fill(0.0);
  for (int i=0; i<n_sigma_; i++) {
    VectorXd col = Xsig.col(i);
    VectorXd d(n_x_);
    double px = col(0);
    double py = col(1);
    double v = col(2);
    double psi = col(3);
    double psi_dot = col(4);
    double nu_a = col(5);
    double nu_psi = col(6);
    double dt2 = 0.5 * (delta_t * delta_t);
    double psi_dt = psi_dot * delta_t;
    if (psi_dot < 0.001) {
      d << v * cos(psi) * delta_t,
	v * sin(psi) * delta_t,
	0,
	psi_dt,
	0;
    } else {
      double m = v / psi_dot;
      d << m * (sin(psi + psi_dt) - sin(psi)),
	m * (-cos(psi + psi_dt) + cos(psi)),
	0,
	psi_dt,
	0;
    }
    VectorXd nu(n_x_);
    nu << dt2 * cos(psi) * nu_a,
      dt2 * sin(psi) * nu_a,
      delta_t * nu_a,
      dt2 * nu_psi,
      delta_t * nu_psi;
    Xsig_pred_.col(i) = col.head(n_x_) + d + nu;
  }

  std::cout << "pred Xsig_pred " << std::endl << Xsig_pred_ << std::endl;
  
  // Predict state
  x_.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }
  
  // Predict covariance
  P_.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    VectorXd col = Xsig_pred_.col(i) - x_;
    col(3) = Tools::normalize_angle(col(3));
    P_ += weights_(i) * col * col.transpose();
  }

  std::cout << "pred x_ " << std::endl << x_ << std::endl;
  std::cout << "pred P_ " << std::endl << P_ << std::endl;

  for (int i=0; i<5; i++)
    assert(P_(4,i) < 10000);
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  int n_z = 2;
  MatrixXd Zsig(n_z, n_sigma_);
  MatrixXd R(n_z, n_z);
  R << (std_laspx_ * std_laspx_), 0,
    0, (std_laspy_ * std_laspy_);

  for (int i=0; i < n_sigma_; i++) {
    VectorXd col = Xsig_pred_.col(i);
    double px = col[0];
    double py = col[1];
    VectorXd zcol(n_z);
    zcol << px, py;
    Zsig.col(i) = zcol;
  }

  // Predicted measurement mean
  VectorXd z_pred(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  // Predicted covariance
  MatrixXd S(n_z, n_z);
  S.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    VectorXd col = Zsig.col(i);
    col -= z_pred;
    S += weights_(i) * col * col.transpose();
  }
  S += R;
  MatrixXd S_inv = S.inverse();

  // Kalman gain
  MatrixXd T(n_x_, n_z);
  T.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::normalize_angle(x_diff(3));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    T += weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K(n_x_, n_z);
  K = T * S_inv;

  // Extract measurement
  VectorXd z_meas = meas_package.raw_measurements_;

  // Update state
  VectorXd z_diff = z_meas - z_pred;
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  std::cout << "lidar x_ " << std::endl << x_ << std::endl;
  std::cout << "lidar P_ " << std::endl << P_ << std::endl;

  // Compute NIS
  NIS_laser_ = z_diff.transpose() * S_inv * z_diff;

  std::cout << "NIS_laser " << NIS_laser_ << std::endl << std::endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  int n_z = 3;
  MatrixXd Zsig(n_z, n_sigma_);
  MatrixXd R(n_z, n_z);
  R << (std_radr_ * std_radr_), 0, 0,
    0, (std_radphi_ * std_radphi_), 0,
    0, 0, (std_radrd_ * std_radrd_);
  
  for (int i=0; i < n_sigma_; i++) {
    VectorXd col = Xsig_pred_.col(i);
    double px = col[0];
    double py = col[1];
    double v = col[2];
    double psi = col[3];
    double rho = sqrt((px * px) + (py * py));
    double phi = atan2(py, px);
    double rho_dot = rho > 0.001 ? ((px * cos(psi) * v) + (py * sin(psi) * v)) / rho : 0;
    VectorXd zcol(n_z);
    zcol << rho, phi, rho_dot;
    Zsig.col(i) = zcol;
  }

  //  std::cout << "Zsig" << std::endl << Zsig << std::endl;
  

  // Predicted measurement mean
  VectorXd z_pred(n_z);  // z_k+1|k
  z_pred.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //  std::cout << "z_pred" << std::endl << z_pred << std::endl;
  
  // Predicted covariance
  MatrixXd S(n_z, n_z);
  for (int i=0; i < n_sigma_; i++) {
    VectorXd d = Zsig.col(i) - z_pred;
    // angle normalization
    d(1) = Tools::normalize_angle(d(1));
    S += (weights_(i) * d * d.transpose());
  }
  S += R;
  MatrixXd S_inv = S.inverse();

  // Kalman gain
  MatrixXd T(n_x_, n_z);
  T.fill(0.0);
  for (int i=0; i < n_sigma_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Tools::normalize_angle(x_diff(3));
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Tools::normalize_angle(z_diff(1));
    T += weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K(n_x_, n_z);
  K = T * S_inv;

  // Extract measurement
  VectorXd z_meas = meas_package.raw_measurements_;

  //  std::cout << "z_meas" << std::endl << z_meas << std::endl;
  
  // Update state
  VectorXd z_diff = z_meas - z_pred;
  z_diff(1) = Tools::normalize_angle(z_diff(1));
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  //  std::cout << "radar z_diff" << std::endl << z_diff << std::endl;
  std::cout << "radar x_ " << std::endl << x_ << std::endl;
  std::cout << "radar P_ " << std::endl << P_ << std::endl;
  
  // Compute NIS
  NIS_radar_ = z_diff.transpose() * S_inv * z_diff;

  std::cout << "NIS_radar " << NIS_radar_ << std::endl << std::endl;
}
