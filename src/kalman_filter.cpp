#include "kalman_filter.h"
#include <iostream>
#include "tools.h"
#include <cmath>
#include <stdlib.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;
using std::vector;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
    x_ = x_in;
    P_ = P_in;
    F_ = F_in;
    H_ = H_in;
    R_ = R_in;
    Q_ = Q_in;
    
}

void KalmanFilter::Predict() {
    /**
     * predict the state
     */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;

    
}

void KalmanFilter::Update(const VectorXd &z) {
    /**
     * update the state by using Kalman Filter equations
     */
    
    VectorXd y = z - H_ * x_;
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    I_ = Eigen::MatrixXd::Identity(x_.size(), x_.size());
    //new state
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
    
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    /**
     * update the state by using Extended Kalman Filter equations
     */
    
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);
    VectorXd h_p  = VectorXd(3);
    
    float rho = sqrt( px*px+py*py );
    float phi  = 0;
    float rho_dot = 0;
    
    //check division by zero
    
    // avoid division by zero
    if(fabs(px) < 0.0001){
        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
    }else{
        phi = atan2(py,px);
    }
    
    
    if (rho < 0.0001) {
        cout << "Error while converting vector x_ to polar coordinates: Division by Zero" << endl;
    }else{
        rho_dot = (px*vx + py*vy) / rho;
    }
    
    h_p << rho, phi, rho_dot;
    
    VectorXd y = z - h_p;
    
    while ( y(1) < -M_PI) {
        y(1)  +=  2*M_PI;
    }
    while (y(1) > M_PI) {
        y(1)  -= 2*M_PI;
    }
    
    if (y(1)  > M_PI || y(1) < -M_PI) {
        cout << "Error" << endl;
    }
    
    MatrixXd Ht = H_.transpose();
    MatrixXd S = H_ * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd K =  P_ * Ht * Si;
    
    I_ = Eigen::MatrixXd::Identity(x_.size(), x_.size());
    //new state
    x_ = x_ + (K * y);
    P_ = (I_ - K * H_) * P_;
}
