#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  //check the validatity of the following inputs:

  if(estimations.size() != ground_truth.size() || estimations.size() == 0){
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }
  
  //Accumulate squared residuals
  for(unsigned int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    //coefficeient-wise multiplication
    residual = residual.array()*residual.array();
    rmse += residual;
  }

  //Calculate teh mean
  rmse = rmse/estimations.size();

  //Calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}