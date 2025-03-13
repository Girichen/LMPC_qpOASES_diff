/*
 * ==================================================
 * Project:    MPC
 * File:       mpc_qpOASES.cpp
 * Author:     Girichen
 * Email:      845031364@qq.com
 * License:    MIT License 
 * Created on: 2025-03-13
 * Description: MPC
 * ==================================================
 */

// Copyright (c) 2025, Girichen
// This code is licensed under the MIT License.
// See the LICENSE file for details.
#ifndef MPC_qp_H_
#define MPC_qp_H_

#include <qpOASES.hpp>
#include <Eigen/Eigen>
#include <vector>
#include "math.h"
#include <iostream>
#include <ctime>
using namespace Eigen;
using namespace std;
USING_NAMESPACE_QPOASES

class MPC_qp
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    MPC_qp(int N,double dt) : N_(N), dt_(dt) {}

    ~MPC_qp() {}

    Vector2d solveMPC(Vector3d& x_0,vector<Vector3d>& X_r,Eigen::Vector2d& odom_data) const {
        if(X_r.size() != static_cast<size_t>(N_)){
            return Vector2d(0,0);
        }

        double last_theta = x_0(2);
        vector<MatrixXd> A_r(N_), B_r(N_), A_multiply2;
        MatrixXd A_multiply1;
        MatrixXd A_multiply2_temp;
        MatrixXd O_r(3*N_,1);
        MatrixXd C_bar = MatrixXd::Identity(3*N_,3*N_);
        MatrixXd Q = MatrixXd::Zero(3 * N_, 3 * N_);
        MatrixXd R = MatrixXd::Zero(2 * N_, 2 * N_);
        MatrixXd A_bar = MatrixXd::Zero(3 * N_, 3);
        MatrixXd B_bar = MatrixXd::Zero(3 * N_, 2 * N_);
        MatrixXd A_r_init = MatrixXd::Identity(3,3);
        MatrixXd B_r_init = MatrixXd::Zero(3,2);

        MatrixXd X_ref(3*N_,1);
        Eigen::MatrixXd O_r_temp(3, 1);
        O_r_temp <<  odom_data(0) * sin(x_0(2)) * dt_ * x_0(2),
                    -odom_data(0) * cos(x_0(2)) * dt_ * x_0(2),
                    0;
        A_r_init(0,2) = -odom_data(0) * sin(x_0(2)) * dt_;
        A_r_init(1,2) = odom_data(0) * cos(x_0(2)) * dt_; 
        B_r_init(0,0) = cos(x_0(2)) * dt_;
        B_r_init(1,0) = sin(x_0(2)) * dt_;
        B_r_init(2,1) = dt_;           
        for(int k = 0; k < N_; ++k){
            Q(k * 3, k * 3) = x_weight_;  
            Q(k * 3 + 1, k * 3 + 1) = y_weight_;  
            Q(k * 3 + 2, k * 3 + 2) = yaw_weight_;  
            R(k * 2, k * 2) = v_weight_; 
            R(k * 2 + 1, k * 2 + 1) = w_weight_;
            if (X_r[k](2) - last_theta > M_PI) {
                X_r[k](2) -= 2 * M_PI;
            } else if (X_r[k](2) - last_theta < -M_PI) {
                X_r[k](2) += 2 * M_PI;
            }  
            last_theta = X_r[k](2);
            X_ref.block<3,1>(k*3,0) = X_r[k];
            
            O_r.block<3,1>(k * 3, 0) = O_r_temp;
            A_r[k] = A_r_init;
            B_r[k] = B_r_init;

            if (k == 0) {
                A_multiply1 = A_r[k];//A(k)
            } else { //k >= 1
                A_multiply1 = A_r[k] * A_multiply1;//A_multiply1 is A_r(k - 1)
                if(k == 1){
                    A_multiply2_temp = A_r[k];//A(k + 1)
                } else {
                    A_multiply2_temp = A_r[k] * A_multiply2_temp;
                }
                A_multiply2.push_back(A_multiply2_temp);//A(k+1) \ A(k+2)A(k+1) \ A(k+3)A(k+2)A(k+1) 
            } 

            int A_multiply2_size = A_multiply2.size();
            for(int i = 0; i < A_multiply2_size; ++i){
                B_bar.block<3,2>(3 * k , 2 * i) = A_multiply2[A_multiply2_size - 1 - i] * B_r[i];//B_bar not finish in this step
                C_bar.block<3,3>(3 * k , 3 * i) = A_multiply2[A_multiply2_size - 1 - i];//B_bar not finish in this step
            }
            
            A_bar.block<3,3>(3 * k, 0) = A_multiply1;//A_bar finished
            B_bar.block<3,2>(3 * k, 2 * k) = B_r[k];//B_bar finished
        }
        
        MatrixXd E = A_bar*x_0+C_bar*O_r-X_ref;
        MatrixXd Hesse = 2 * (B_bar.transpose() * Q * B_bar + R);
        VectorXd gradient = 2 * B_bar.transpose() * Q * E;
       
        real_t A[2*N_*2*N_],lb[2*N_],ub[2*N_];
        QProblem mpc_qp_solver(2 * N_, 1);
        mpc_qp_solver.setPrintLevel(PrintLevel::PL_NONE);

        for(int i=0;i<2*N_;i++)
        {
            if (i % 2 == 0) {
                 lb[i] = v_min_ ;
                 ub[i] = v_max_ ;
            } else {
                 lb[i] = w_min_ ;
                 ub[i] = w_max_ ;
            }
        }

        int_t nWSR = 200;
        mpc_qp_solver.init(Hesse.data(),gradient.data(),A,lb,ub,nullptr,nullptr,nWSR);
        real_t x_solution[2 * N_];
        mpc_qp_solver.getPrimalSolution(x_solution);
        Vector2d u_k;
        u_k(0) = x_solution[0];
        u_k(1) = x_solution[1];
        //std::cout << "Program execution time: " << duration << " seconds" << std::endl;
        return u_k;
    }

private:
    int N_ = 10;
    double dt_ = 0.1;
    double v_min_ = -0.6;
    double v_max_ = 0.6;
    double w_min_ = -0.6;
    double w_max_ = 0.6;
    double x_weight_ = 30;
    double y_weight_ = 30;
    double yaw_weight_ = 10;
    double v_weight_ = 5;
    double w_weight_ = 5;
};

#endif
