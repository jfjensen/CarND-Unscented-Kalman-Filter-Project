#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

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
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_; // ?

  weights_ = VectorXd(2 * n_aug_ + 1);


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

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
 
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float ro_dot = meas_package.raw_measurements_(2);
      float px = ro * cos(phi);
      float py = ro * sin(phi);
      float vx = ro_dot * cos(phi);
      float vy = ro_dot * sin(phi);
      x_ << px, py, 0, 0, 0;

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */

      float px = meas_package.raw_measurements_(0);
      float py = meas_package.raw_measurements_(1);
      x_ << px, py, 0, 0, 0;
      
    }
    time_us_ = meas_package.timestamp_;
    // done initializing, no need to predict or update

    P_.fill(0.0);

    int n_sig = 2 * n_aug_ + 1;
    double w1 = lambda_ / (lambda_ + n_aug_);
    double w2 = 1 / (2 * (lambda_ + n_aug_));
    weights_(0) = w1;
    for (unsigned int i = 1; i < n_sig; i++)
    {
      weights_(i) = w2;
    }

    is_initialized_ = true;
    return;
  }

  

  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  //cout << "So far so good..." << dt << endl;
  Prediction(dt);


  // if (use_radar_)
  // {
  //   // Radar updates
  //   UpdateRadar(meas_package);
   
  // }
  // if (use_laser_)
  // {
  //   // Laser updates
  //   UpdateLidar(meas_package);
    
  // }


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

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  int n_sig = 2 * n_aug_ + 1;
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig);

  x_aug.head(n_x_) = x_;
  x_aug(n_aug_ - 2) = 0;
  x_aug(n_aug_ - 1) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner( n_x_, n_x_ ) = P_;
  P_aug(n_aug_ - 2, n_aug_ - 2) = std_a_ * std_a_;
  P_aug(n_aug_ - 1, n_aug_ - 1) = std_yawdd_ * std_yawdd_;

  MatrixXd A = P_aug.llt().matrixL();

  Xsig_aug.col(0) = x_aug;

  float sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_);

  //std::cout << "so far so good..." << std::endl;

  for (unsigned int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i + 1)          = x_aug + sqrt_lambda_n_aug * A.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_lambda_n_aug * A.col(i);
  }

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig);

  for (unsigned int i = 0; i < n_sig; i++)
  {
    VectorXd x_k(n_x_);
    VectorXd sigm(n_aug_);
    sigm = Xsig_aug.col(i);

    float px      = sigm(0);
    float py      = sigm(1);
    float v       = sigm(2);
    float psi     = sigm(3);
    float psi_dot = sigm(4);
    float v_a     = sigm(5);
    float psi_dot_dot = sigm(6);

    x_k << px, py, v, psi, psi_dot;

    VectorXd term1(n_x_);
    VectorXd term2(n_x_);

    float dt = delta_t;

    if (v_a == 0 and psi_dot_dot == 0)
    {
      term2 << 0,0,0,0,0;
    }
    else
    {
      term2 << 0.5 * dt * dt * cos(psi) * v_a,
               0.5 * dt * dt * sin(psi) * v_a,
               dt * v_a,
               0.5 * dt * dt * psi_dot_dot,
               dt * psi_dot_dot;
    }

    if (fabs(psi_dot) > 0.001)
    {
      term1 << (v / psi_dot) * (sin(psi + (psi_dot * dt)) - sin(psi)),
               (v / psi_dot) * (-cos(psi + (psi_dot * dt)) + cos(psi)),
               0,
               psi_dot * dt,
               0;
    }
    else
    {
      term1 << v * cos(psi) * dt,
               v * sin(psi) * dt,
               0,
               psi_dot * dt,
               0;
    }

    Xsig_pred.col(i) = x_k + term1 + term2;


  }

  VectorXd x_new = VectorXd(n_x_);
  x_new.fill(0.0);
  for (unsigned int i = 0; i < n_sig; i++)
  {
    x_new = x_new + (weights_(i) * Xsig_pred.col(i));
  }

  x_ = x_new;

  MatrixXd P_new = MatrixXd(n_x_, n_x_);
  P_new.fill(0.0);
  for (unsigned int i = 0; i < n_sig; i++)
  {
    
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_new = P_new + weights_(i) * x_diff * x_diff.transpose();
  }

  P_ = P_new;

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
