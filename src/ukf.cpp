#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  
  is_initialized_ = false;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_.setZero();

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3; // set this to 3 initially

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  // time when the state is true?
  time_us_ = 0;
  Xsig_pred_.setZero();
  // number of arguments 
  n_x_ = 5;
  n_aug_ = 7;
  // sigma point spreading parameter
  lambda_ = 3 - n_x_; // default vaule
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*************
  Initialization 
  **************/
  if (!is_initialized_) {
    // set the timestamp
    time_us_ = meas_package.timestamp_;
    // initialize x 
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 1, 1, 1;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double ro = meas_package.raw_measurements_(0);
      double theta = meas_package.raw_measurements_(1);
      x_ << ro*cos(theta), ro*sin(theta), 1, 1, 1;
    }

    // initialize P
    P_ << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

    // done initialization
    is_initialized_ = true;
    return;
  }

  /*************
  Prediction 
  **************/
  double delta_t = meas_package.timestamp_ - time_us_;
  time_us_ = meas_package.timestamp_;
  Prediction(delta_t);

  /*************
  Update
  **************/
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*
  Generate augmented sigma points.
  */

  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(2, 2) = Q;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }

  /*
  Predict sigma points.
  */

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  for (int i = 0; i < 2*n_aug+1; i++) {
      VectorXd current_x = Xsig_aug.col(i);
      double p_x = current_x(0);
      double p_y = current_x(1);
      double v = current_x(2);
      double phi = current_x(3);
      double phi_dot = current_x(4);
      double noise_a = current_x(5);
      double noise_phi = current_x(6);
      
      // a and b are two vectors in the prediction equation
      VectorXd a = VectorXd::Zero(n_x);
      VectorXd b = VectorXd::Zero(n_x);
      
      if (fabs(phi_dot) < 0.001) {
        a(0) = v*cos(phi)*delta_t;
        a(1) = v*sin(phi)*delta_t;
      }
      else {
        a(0) = (v/phi_dot) * (sin(phi+phi_dot*delta_t) - sin(phi));
        a(1) = (v/phi_dot) * (-cos(phi+phi_dot*delta_t) + cos(phi));
        a(3) = phi_dot*delta_t;
      }
      
      b(0) = 0.5*delta_t*delta_t * cos(phi) * noise_a;
      b(1) = 0.5*delta_t*delta_t * sin(phi) * noise_a; 
      b(2) = delta_t * noise_a;
      b(3) = 0.5*delta_t*delta_t * noise_phi;
      b(4) = delta_t * noise_phi;
      
      Xsig_pred.col(i) = Xsig_aug.col(i).head(n_x) + a + b;
  }

  /*
  Calculate mean state and variance 
  */
 
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  int num_points = 2*n_aug+1;
  //set weights
  weights(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < num_points; i++) {
      weights(i) = 0.5 * 1 / (lambda_ + n_aug_);
  }
  //predict state mean
  x_.setZero();
  P_.setZero();
  for (int i = 0; i < num_points; i++) {
      x_ += weights(i) * Xsig_pred.col(i);
  }
  //predict state covariance matrix
  for (int i = 0; i < num_points; i++) {
      P_ += weights(i) * (Xsig_pred.col(i) - x_) * (Xsig_pred.col(i) - x_).transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}

/**
 * Generate sigma points
 */
void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

  //create augmented mean state
  x_aug.setZero();
  x_aug.head(n_x_) = x_;
  
  //create augmented covariance matrix
  P_aug.setZero();
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_*std_a_, 0,
       0, std_yawdd_*std_yawdd_;
  P_aug.bottomRightCorner(2, 2) = Q;
  
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug; i++) {
      Xsig_aug.col(i+1) = x_aug + sqrt(lambda_ + n_aug_) * A.col(i);
      Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda_ + n_aug_) * A.col(i);
  }
  //write result
  *Xsig_out = Xsig_aug;
}

