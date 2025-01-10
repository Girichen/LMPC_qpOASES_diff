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
    MPC_qp(int N,double dt) : N_(N), dt_(dt){

    }

    ~MPC_qp() {}

    Vector2d solveMPC(Vector3d& x_0,vector<Vector3d>& X_r,vector<Vector2d>& U_r){
        clock_t start = std::clock();
        vector<MatrixXd> A_r(N_), B_r(N_), A_multiply2;
        MatrixXd A_multiply1;
        MatrixXd A_multiply2_temp;
        Vector3d X_0 = x_0 - X_r[0];
        if(X_0(2) < -M_PI){
            X_0(2) += 2 * M_PI;
        } else if (X_0(2) > M_PI){
            X_0(2) -= 2 * M_PI;
        } else {
            
        }
        MatrixXd Q = MatrixXd::Identity(3 * N_, 3 * N_) * omega0_;
        MatrixXd R = MatrixXd::Identity(2 * N_, 2*N_) * omega1_;
        MatrixXd A_bar = MatrixXd::Zero(3 * N_, 3);
        MatrixXd B_bar = MatrixXd::Zero(3 * N_, 2 * N_);
        MatrixXd A_r_init = MatrixXd::Identity(3,3);
        MatrixXd B_r_init = MatrixXd::Zero(3,2);

        for(int k = 0; k < N_; ++k){
            A_r[k] = A_r_init;
            B_r[k] = B_r_init;

            A_r[k](0,2) = -U_r[k](0) * sin(X_r[k](2));
            A_r[k](1,2) = U_r[k](0) * cos(X_r[k](2));

            B_r[k](0,0) = cos(X_r[k](2)) * dt_;
            B_r[k](1,0) = sin(X_r[k](2)) * dt_;
            B_r[k](2,1) = dt_;

            if (k == 0) {
                A_multiply1 = A_r[k];//A(k)
            } else {
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
            }
            
            A_bar.block<3,3>(3 * k, 0) = A_multiply1;//A_bar finished
            B_bar.block<3,2>(3 * k, 2 * k) = B_r[k];//B_bar finished
        }

        MatrixXd Hesse = 2 * (B_bar.transpose() * Q * B_bar + R);
        VectorXd gradient = 2 * B_bar.transpose() * Q * A_bar * X_0;
       
        real_t H[2*N_*2*N_],g[2*N_],A[2*N_*2*N_],lb[2*N_],ub[2*N_];
        QProblem mpc_qp_solver(2 * N_, 1);
        mpc_qp_solver.setPrintLevel(PrintLevel::PL_NONE);

        for(int i=0;i<2*N_;i++)
        {
            g[i] = gradient(i);
            if(i % 2 == 0){
                lb[i] = v_min_ - U_r[0](0);
                ub[i] = v_max_ - U_r[0](0);
            } else {
                lb[i] = w_min_ - U_r[0](1);
                ub[i] = w_max_ - U_r[0](1);
            }
            for(int j = 0; j < 2 * N_; j++)
            {
                H[i*2*N_+j] = Hesse(i,j);
            }
        }

        int_t nWSR = 800;
        mpc_qp_solver.init(H,g,A,lb,ub,nullptr,nullptr,nWSR);
        real_t x_solution[2 * N_];
        mpc_qp_solver.getPrimalSolution(x_solution);
        Vector2d u_k;
        u_k(0) = x_solution[0];
        u_k(1) = x_solution[1];
        clock_t end = std::clock();
        double duration = double(end - start) / CLOCKS_PER_SEC;
        //std::cout << "Program execution time: " << duration << " seconds" << std::endl;
        return u_k;
    }
private:
    int N_ = 10;
    double dt_ = 0.1;
    double omega0_ = 0.01;
    double omega1_ = 5.0;
    double v_min_ = -0.6;
    double v_max_ = 0.6;
    double w_min_ = -0.6;
    double w_max_ = 0.6;
};

#endif
