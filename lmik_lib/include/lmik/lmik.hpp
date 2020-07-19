/*********************************************************************
* Copyright 2020 Yuki Onishi
* 
* This software is released under the MIT license.
* http://opensource.org/licenses/mit-license.php
*********************************************************************/

#ifndef LMIK_LMIK_HPP_
#define LMIK_LMIK_HPP_

#include <functional>

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Cholesky>

namespace lmik {

using Eigen::MatrixXd;
using Eigen::VectorXd;

class LMIK {
public:
    LMIK();
    ~LMIK();

    void setBasicJacobianFunc(const std::function<MatrixXd&(const VectorXd&)> func) {
        calcBasicJacobian = func;
    }
    void setForwardKinematicsFunc(const std::function<VectorXd&(const VectorXd&)> func) {
        calcForwardKinematics = func;
    }

    void setDampingFactor(const double df) {
        df_ = df;
    }
    void setWeightMatrix(const MatrixXd& W) {
        W_E_ = W;
    }

    bool solve(VectorXd& q_init, VectorXd& x_d, uint32_t iter = 500, double threshold = 10e-5);

protected:
    std::function<MatrixXd&(const VectorXd&)> calcBasicJacobian;
    std::function<VectorXd&(const VectorXd&)> calcForwardKinematics;

private:
    double   df_;   // damping factor
    MatrixXd W_E_;  // weights
};

}  // namespace lmik

#endif  // LMIK_LMIK_HPP_
