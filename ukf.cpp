#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double PI = 3.1415926;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  //state of initialization
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = PI/8;

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

  // dimension of original state
  n_x_ = 5;

  // dimension of augmented state
  n_aug_ = 7;

  // set the value of lambda
  lambda_ = 3 - n_aug_;

  // Matrix of sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  
  // augmented state
  x_aug_ = VectorXd(n_aug_);

  // augmented convariance matrix
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // weight vector
  weights_ = VectorXd(2 * n_aug_ + 1);

  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2* n_aug_ + 1; i++){
	  double weight_i = 1/(2*(lambda_+n_aug_));
	  weights_(i) = weight_i;
  }

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
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
	if (!is_initialized_){
		// initialize the state and covariance
		x_.fill(0.0);
		P_<<1,0,0,0,0,
			0,1,0,0,0,
			0,0,1,0,0,
			0,0,0,1,0,
			0,0,0,0,1;
		P_aug_.fill(0.0);
		P_aug_.topLeftCorner(5,5) = P_;
		P_aug_(5,5) = std_a_ * std_a_;
		P_aug_(6,6) = std_yawdd_ * std_yawdd_;

		// read the first data
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
			double rho = meas_package.raw_measurements_(0);
			double phi = meas_package.raw_measurements_(1);
			double rho_dot = meas_package.raw_measurements_(2);
			double px = rho*cos(phi);
			double py = rho*sin(phi);
			double vx = rho_dot*cos(phi);
			double vy = rho_dot*sin(phi);
			double v = sqrt(vx*vx+vy*vy);

			x_ << px, py, v, 0, 0;
			cout<<"Radar: x_: "<<x_<<endl;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
			x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0, 0, 0;
		}

		// aviod division by zero
		if (fabs(x_(0)) < 0.0001 and fabs(x_(1)) < 0.0001) {
            x_(0) = 0.0001;
            x_(1) = 0.0001;
		}

		// augmented state
		x_aug_.head(5) = x_;
		x_aug_(5) = 0;
		x_aug_(6) = 0;
		
		time_us_ = meas_package.timestamp_;

		is_initialized_ = true;

		return;
	}

	double delta_t = (meas_package.timestamp_ - time_us_)/1000000.0;
	time_us_ = meas_package.timestamp_;

	Prediction(delta_t);

	if (meas_package.sensor_type_ == MeasurementPackage::RADAR){
		UpdateRadar(meas_package);
	}
	else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
		UpdateLidar(meas_package);
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
	// set the augmented state and covariance
	P_aug_.fill(0.0);
	P_aug_.topLeftCorner(5,5) = P_;
	P_aug_(5,5) = std_a_ * std_a_;
	P_aug_(6,6) = std_yawdd_ * std_yawdd_;

	x_aug_.head(5) = x_;
	x_aug_(5) = 0;
	x_aug_(6) = 0;

	// choose sigma points
	MatrixXd L = P_aug_.llt().matrixL();

	MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2*n_aug_ + 1);

	Xsig_aug_.col(0) = x_aug_;

	for (int i = 0; i < n_aug_; i++){
		Xsig_aug_.col(i+1) = x_aug_ + sqrt(lambda_ + n_aug_) * L.col(i);
		Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_ + n_aug_) * L.col(i);
	}

	// predict sigma points
        for (int i = 0; i< 2*n_aug_+1; i++)
    {
        //extract values for better readability
        double px = Xsig_aug_(0,i);
        double py = Xsig_aug_(1,i);
        double v = Xsig_aug_(2,i);
        double yaw = Xsig_aug_(3,i);
        double yawd = Xsig_aug_(4,i);
        double nu_a = Xsig_aug_(5,i);
        double nu_yawdd = Xsig_aug_(6,i);
        
        //Predicted values
        double px_p, py_p;
        
        //Avoid division by zero
        if (fabs(yawd) > 0.0001) {
            px_p = px + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        }
        else {
            px_p = px + v*delta_t*cos(yaw);
            py_p = py + v*delta_t*sin(yaw);
        }
        
        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;
        
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v + nu_a * delta_t;
        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;
        
        Xsig_pred_(0,i) = px_p;
        Xsig_pred_(1,i) = py_p;
        Xsig_pred_(2,i) = v_p;
        Xsig_pred_(3,i) = yaw_p;
        Xsig_pred_(4,i) = yawd_p;
        
        //Predicted state mean
        x_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            x_ = x_ + weights_(i) * Xsig_pred_.col(i);
        }
        
        // predicted state covariance matrix
        P_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {
            
            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;
            // angle normalization
            while (x_diff(3)> PI) x_diff(3)-=2.*PI;
            while (x_diff(3)< -PI) x_diff(3)+=2.*PI;
            
            P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
        }
    }
	cout << "Xsig_pred_: "<<Xsig_pred_<<endl;
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
	// dimension of measurement
	int n_z=2;
    VectorXd z = meas_package.raw_measurements_;
    cout<<"z (Lidar): "<<z<<endl;
    
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    
    // find sigma points
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        Zsig.col(i) << p_x, p_y;
    }
    cout<<"Zsig (Lidar): "<<Zsig<<endl;
    
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    z_pred=Zsig * weights_;
    
    cout<<"z_pred (Lidar): "<<z_pred<<endl;
    
    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    // measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;

    S = S + R;
    
    cout << "S (Lidar): "<<S<<endl;
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    //residual
    VectorXd z_diff = z - z_pred;
    
    cout << "z_diff (Lidar) "<<z_diff<<endl;
    
    //update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
    
    //NIS Update
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

	cout<< "Ladar_NIS: "<<NIS_laser_<<endl;
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
	// dimension of measurement
	int n_z=3;
    VectorXd z = meas_package.raw_measurements_;
    cout<<"z (Radar): "<<z<<endl;
    
	// sigma points
    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
    Zsig.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        
        // extract values for better readibility
        double p_x = Xsig_pred_(0,i);
        double p_y = Xsig_pred_(1,i);
        double v  = Xsig_pred_(2,i);
        double yaw = Xsig_pred_(3,i);
        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;
        
        if (fabs(p_x)<0.0001) {
            p_x=0.0001;
        }
        if (fabs(p_y) < 0.0001) {
            p_y = 0.0001;
        }
        
        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
        Zsig(1,i) = atan2(p_y,p_x);                                 //phi
        Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
    }
    cout<<"Zsig (Radar): "<<Zsig<<endl;
    
    VectorXd z_pred = VectorXd(n_z);
    z_pred.fill(0.0);
    z_pred=Zsig*weights_;
    
    cout<<"z_pred (Radar): "<<z_pred<<endl;
    
    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        //angle normalization
        while (z_diff(1)>PI)z_diff(1) -= 2.*PI;
		while (z_diff(1)<-PI)z_diff(1) += 2.*PI;
        
        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    
    // measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<    std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    
    S = S + R;
    cout << "S (Radar): "<<S<<endl;
    
    MatrixXd Tc = MatrixXd(n_x_, n_z);
    Tc.fill(0.0);
    
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        while (z_diff(1)>PI)z_diff(1) -= 2.*PI;
		while (z_diff(1)<-PI)z_diff(1) += 2.*PI;

        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // angle normalization
        while (x_diff(3)>PI)z_diff(1) -= 2.*PI;
		while (x_diff(3)<-PI)z_diff(1) += 2.*PI;
        
        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }
    
    // kalman gain K;
    MatrixXd K = Tc * S.inverse();
    
    VectorXd z_diff=z-z_pred;
    // angle normalization
	while (z_diff(1)>PI)z_diff(1) -= 2.*PI;
	while (z_diff(1)<-PI)z_diff(1) += 2.*PI;
    
    // update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();
    
    // NIS Update
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

	cout<<"Radar_NIS: "<<NIS_radar_<<endl;
}
