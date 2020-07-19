/*********************************************************************
* Copyright 2020 Yuki Onishi
* 
* This software is released under the MIT license.
* http://opensource.org/licenses/mit-license.php
*********************************************************************/

#include "lmik/lmik.hpp"

namespace lmik {

LMIK::LMIK(/* args */)
    : df_(10e-1)
    , W_E_(MatrixXd::Identity(1, 1)) {
}

LMIK::~LMIK() {
}

bool LMIK::solve(VectorXd& q_init, VectorXd& x_d, uint32_t iter, double threshold) {
    uint32_t n = q_init.size();
    uint32_t m = x_d.size();
    VectorXd q = q_init;

    // allocates and checks sizes of vectors and matrices
    VectorXd e = - calcForwardKinematics(q);  // error vector
    if(e.size() != m) {
        return false;
    } else {
        e += x_d;  // initalization
    }

    MatrixXd J = calcBasicJacobian(q);  // Jacobian matrix
    if(J.rows() != m || J.cols() != n) {
        return false;
    }

    if(W_E_.rows() != m || W_E_.cols() != m) {
        return false;
    }

    // solves once
    VectorXd g   = J.transpose() * (W_E_ * e).eval();
    MatrixXd W_N = (e.norm() + df_) * MatrixXd::Identity(n, n);
    MatrixXd H   = J.transpose() * W_E_ * J + W_N;
    
    Eigen::LDLT<MatrixXd> H_ldlt(H);
    VectorXd dq = H_ldlt.solve(g);
    q += dq;

    // iteration
    for(uint32_t i = 0; i < iter - 1; i++) {
        if(dq.norm() < threshold) {
            break;
        }

        e = x_d - calcForwardKinematics(q);
        J = calcBasicJacobian(q);

        g   = J.transpose() * (W_E_ * e).eval();
        W_N = (e.norm() + df_) * MatrixXd::Identity(n, n);
        H   = J.transpose() * W_E_ * J + W_N;
        H_ldlt.compute(H);
        dq = H_ldlt.solve(g);
        q += dq;
    }

    // assings the solution to the arg.
    q_init = q;

    return true;
}

}  // namespace lmik
