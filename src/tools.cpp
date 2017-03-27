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

    // TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	int est_size = estimations.size();
	int gt_size  = ground_truth.size();

	if (est_size == 0 or gt_size == 0)
	{
		std::cout << "Error vectors should not be zero!" << std::endl;
		return rmse;
	}
	//  * the estimation vector size should equal ground truth vector size
	if (est_size != gt_size)
	{
		std::cout << "Vector sizes should be the same!" << std::endl;
		return rmse;
	}

	//std::cout << "rmse 1: " << rmse << std::endl;
	//accumulate squared residuals
	
	for(int i=0; i < estimations.size(); ++i){
		//std::cout << "============================================================" << std::endl;
		//std::cout << "estimations : " << estimations[i] << " i: " << i << std::endl;
		//std::cout << "ground_truth : " << ground_truth[i] << " i: " << i << std::endl;
    	VectorXd temp(4);
    	temp = (estimations[i].array() - ground_truth[i].array());
    	//std::cout << "temp : " << temp << " i: " << i << std::endl;
    	rmse = rmse.array() + (temp.array() * temp.array());
    	//std::cout << "============================================================" << std::endl;
	}

	//std::cout << "rmse 2: " << rmse << std::endl;
	//calculate the mean
	rmse = rmse.array()/est_size;
	//std::cout << "rmse 3: " << rmse << std::endl;

	//calculate the squared root
	rmse  = rmse.array().sqrt();
	//std::cout << "rmse 4: " << rmse << std::endl;

	//return the result
	return rmse;
}
