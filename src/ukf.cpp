#include "ukf.h"
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


  is_initialized_ = false;
  //time_us_= 0;
  time_us_=0;
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);



  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.6;//30,3,1,0.8,0.7,0.8,0.8,0.8,0.8

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;//30,3,1,0.8,0.7,0.7,0.6,0.5,0.8

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;//0.15 0.12

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;//0.15 0.12

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;//0.3  0.2 0.25 0.26

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;//0.03 0.02 0.025 0.026

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;//0.3 0.2 0.25 0.26

  n_x_=5;
  n_aug_=7;
  lambda_=3-n_aug_;

  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  weights_=VectorXd(2*n_aug_+1);

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
  if (!is_initialized_) {
     // first measurement

    //cout << "UKF Initialization!! " << endl;
    time_us_ = meas_package.timestamp_;
    //cout<<"time_us_ is :"<<time_us_<<endl;

    //cout << "UKF: " << endl;

    //Eigen::VectorXd x_in;
    //VectorXd x_ = VectorXd(5);
    //cout<<x_ << 0, 0, 0, 0, 0;

      P_ = MatrixXd::Identity(5,5);
      P_(0,0)=0.15;
      P_(1,1)=0.15;
      //cout<<"P_\n"<<P_<<endl;
      /*
      P_.fill(0.0);
    P_(0,0)=1.0;
    P_(1,1)=1.0;
    P_(2,2)=1.0;
    P_(3,3)=1.0;
    P_(4,4)=1.0;
       */
      //cout<<"P_ is "<<P_<<endl;

    //Eigen::MatrixXd H_in;
    //cout<<"meas_package.sensor_type_ is:"<<meas_package.sensor_type_<<endl;

    double weight_0 = lambda_/(lambda_+n_aug_);
      //cout<<"weight_0 is "<<weight_0<<endl;

    weights_(0) = weight_0;
      //cout<<"weights_(0) is "<<weights_(0)<<endl;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights_
      double weight = 0.5 / (n_aug_ + lambda_);
        //cout<<"weight is "<<weight<<endl;
      weights_(i) = weight;
    }
    //cout<<"weights is "<<weights_<<endl;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //cout<<"RADAR Initialization!!"<<endl;
      float ro=meas_package.raw_measurements_(0);
      float phi=meas_package.raw_measurements_(1);
      //float ro_dot=meas_package.raw_measurements_(2);
      x_ << ro*cos(phi),ro*sin(phi), 0.0, 0.0,0.0;
      //time_us_ = meas_package.timestamp_;
      //is_initialized_ = true;
      //std::cout<<"efk_.x_:"<<ekf_.x_<<endl;


    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //cout<<"LIDAR Initialization!!"<<endl;
      //cout<<"meas_package.raw_measurements_ is:"<<meas_package.raw_measurements_<<endl;
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0,0.0;
      //ekf_.x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0;
      //cout<<"x_ is :"<<x_<<endl;

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


  //ekf_.Predict();
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
  time_us_ = meas_package.timestamp_;
  //cout<<"dt is:"<<dt<<endl;

  Prediction(dt);


  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    //UpdateEKF(meas_package.raw_measurements_);

      UpdateRadar(meas_package);
      //std::cout<<"after RADAR UKF Update"<<endl;


  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
    // Laser updates


    UpdateLidar(meas_package);
    //std::cout<<"After Laser UKF Update"<<endl;
  }

  // print the output
  //cout << "x_ = " << x_ << endl;
  //cout << "P_ = " << P_ << endl;


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
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0.0;
  x_aug(6) = 0.0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //cout<<"Xsig_aug is "<<Xsig_aug<<endl;

  //MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug__ + 1);

  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    const double p_x       = Xsig_aug(0,i);
    const double p_y       = Xsig_aug(1,i);
    const double v         = Xsig_aug(2,i);
    const double yaw       = Xsig_aug(3,i);
    const double yawd      = Xsig_aug(4,i);
    const double nu_a      = Xsig_aug(5,i);
    const double nu_yawdd  = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero


    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column

    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  //cout<<"Xsig_pred_ is "<<Xsig_pred_<<endl;

  /*
  //set weights_
  for(int i=0;i<2*n_aug_+1;i++){
    if(i==0){
      weights_(i)=lambda/(lambda+n_aug_);
    }
    else{
      weights_(i)=1.0/(2.0*(lambda+n_aug_));
    }
  }
   */
  //predict state mean
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    x_pred +=weights_(i)*Xsig_pred_.col(i);
      //x_pred+=Xsig_pred_*weights_;
  }
  //cout<<"x_pred is "<<x_pred<<endl;

  //update x_ for x_pred
  x_=x_pred;


  MatrixXd P_pred = MatrixXd(n_x_, n_x_);

  //predict state covariance matrix
  P_pred.fill(0.0);
  VectorXd x_diff=VectorXd(n_x_);
  for(int i=0;i<2*n_aug_+1;i++){
    x_diff=Xsig_pred_.col(i)- x_pred;
    P_pred+=weights_(i)*x_diff*x_diff.transpose();
  }

  //cout<<"P_pred is "<<P_pred<<endl;

  //update P_ for P_pred
  P_=P_pred;
    //cout<<"Prediction x_ is \n"<<x_<<endl;
    //cout<<"Xsig_pred is \n"<<Xsig_pred_<<endl;
    //cout<<"Prediction P_ is \n"<<P_<<endl;


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
  //cout<<"enter Lidar Update\n";
  const VectorXd z=meas_package.raw_measurements_;
  const int n_z_=2;




  //Calculat measurement sigment point
  MatrixXd Zsig;
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);

    // measurement model
    Zsig(0, i) = p_x;
    Zsig(1, i) = p_y;
  }
    //cout<<"Zsig is "<<Zsig<<endl;

  // Calculate mean predicted measurement
    //mean predicted measurement
  VectorXd z_pred;
  z_pred = VectorXd(n_z_);

  z_pred=Zsig*weights_;
    //cout<<"weights_ is "<<weights_<<endl;

  /*
  z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
  z_pred=z_pred+weights_(i)*Zsig.col(i);
  }
   */


  //Calculate measurement covariance matrix

  MatrixXd S = MatrixXd(n_z_, n_z_);
  //S.fill(0.0);
  S << std_laspx_*std_laspx_, 0,
          0, std_laspy_*std_laspy_;
  //S(2,2)=std_radrd_*std_radrd_;

  for(int i=0;i<2*n_aug_+1;i++){
  VectorXd z_diff=Zsig.col(i)-z_pred;
  //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  S +=weights_(i)*z_diff*z_diff.transpose();
  }
  //cout<<"S is "<<S<<endl;
  // creat Matrix Tc for Cross-correlation Matrix
  MatrixXd Tc;
  Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){

    VectorXd z_diff = Zsig.col(i) - z_pred;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
    //cout<<"Tc is "<<Tc<<endl;

  //calculate Kalman gain K;
  MatrixXd K =Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff=z-z_pred;
  //while(z_diff(1)>M_PI) z_diff(1)-=2.0*M_PI;
  //while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();


  //Calculate NIS
    //cout<<"Lidar x_ is "<<x_<<endl;
  //cout<<"Lidar z is "<<z<<endl;
  //cout<<"Lidar z_pred is "<<z_pred<<endl;
  double NIS=(z-z_pred).transpose()*S.inverse()*(z-z_pred);
  //cout<<"Lidar NIS is \n"<<NIS<<endl;
    cout<<NIS<<endl;
  //  cout<<"Finish Lidar Update"<<endl;
}
  //Calculate NIS
    //MatrixXd H_= MatrixXd(2, 5);
    //H_ << 1, 0, 0, 0, 0,
      //      0, 1, 0, 0, 0;
    //cout<<"H_ is "<<H_<<endl;
    //measurement covariance matrix S

  /*
    MatrixXd S;
    S = MatrixXd(n_z_, n_z_);

    S.fill(0.0);
    //cout<<"S zero is "<<S<<endl;
    S(0,0)=std_laspx_*std_laspx_;
    S(1,1)=std_laspy_*std_laspy_;
    //cout<<"S is "<<S<<endl;
    //cout<<"x_ is "<<x_<<endl;
    VectorXd z_pred = H_ * x_;
    //cout<<"z_pred is "<<z_pred<<endl;
    VectorXd y = z - z_pred;
    cout<<"y is "<<y<<endl;
    MatrixXd Ht = H_.transpose();
    MatrixXd S_ = H_ * P_ * Ht + S;
    MatrixXd Si = S_.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;
    //cout<<"K is "<<K<<endl;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
    cout<<"x_ is "<<x_<<endl;
    cout<<"P_ is "<<P_<<endl;

    //Calculate NIS
    cout<<"Lidar z is "<<z<<endl;
    cout<<"Lidar z_pred is "<<z_pred<<endl;
    //double NIS=(z-z_pred).transpose()*S.inverse()*(z-z_pred);
    //cout<<NIS<<endl;

    */


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
  //cout<<"Enter Radar Update!\n";
  //define const value for Radar measurement value
  const VectorXd z=meas_package.raw_measurements_;


  const int n_z_=3;

  //mean predicted measurement
  VectorXd z_pred;
  z_pred = VectorXd(n_z_);

  //measurement covariance matrix S
  MatrixXd S;
  S = MatrixXd(n_z_, n_z_);


  //Calculat measurement sigment point
  MatrixXd Zsig;
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

  for(int i=0;i<2*n_aug_+1;i++){
    double px=Xsig_pred_(0,i);
    double py=Xsig_pred_(1,i);
    double v_=Xsig_pred_(2,i);
    double fai_=Xsig_pred_(3,i);
    Zsig(0,i)=sqrt(px*px+py*py);
    Zsig(1,i)=atan2(py,px);
    if(Zsig(0,i)<0.0001){
      std::cout<<"Error rho_dot is too big\n";
      Zsig(2,i)=0.0;
    }
    else{
      Zsig(2,i)=(px*cos(fai_)*v_+py*sin(fai_)*v_)/Zsig(0,i);
    }
  }

  // Calculate mean predicted measurement
  z_pred.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){
    z_pred=z_pred+weights_(i)*Zsig.col(i);
  }
    //cout<<"Radar z is \n"<<z<<endl;
    //cout<<"Radar z_pred is \n"<<z_pred<<endl;

  //Calculate measurement covariance matrix
  S.fill(0.0);
  S(0,0)=std_radr_*std_radr_;
  S(1,1)=std_radphi_*std_radphi_;
  S(2,2)=std_radrd_*std_radrd_;

  for(int i=0;i<2*n_aug_+1;i++){
    VectorXd z_diff=Zsig.col(i)-z_pred;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    S +=weights_(i)*z_diff*z_diff.transpose();
  }

  // creat Matrix Tc for Cross-correlation Matrix
  MatrixXd Tc;
  Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);
  for(int i=0;i<2*n_aug_+1;i++){

    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K =Tc*S.inverse();

  //update state mean and covariance matrix
  VectorXd z_diff=z-z_pred;
  while(z_diff(1)>M_PI) z_diff(1)-=2.0*M_PI;
  while(z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;


  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  //cout<<"Radar x_ is "<<x_<<endl;


  //Calculate NIS
  //cout<<"Radar z is "<<z<<endl;
  //cout<<"Radar z_pred is "<<z_pred<<endl;
  double NIS=(z-z_pred).transpose()*S.inverse()*(z-z_pred);
  cout<<NIS<<endl;


}
