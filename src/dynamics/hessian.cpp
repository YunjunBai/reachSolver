#include "NonlinearSys.h"
#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include "derivatives.h"

namespace reachSolver{

// template class NonlinearSys<double>;

// template <typename Number>
// std::vector<IntervalMatrix> NonlinearSys<Number>::hessian_(IntervalMatrix im1, IntervalMatrix im2){
std::vector<IntervalMatrix> hessian(mFile_f_type mFile_f_, IntervalMatrix im1, IntervalMatrix im2){
    typedef autodiff::dual2nd (*func_p)(autodiff::ArrayXdual2nd& vector1, autodiff::ArrayXdual2nd& vector2);
    func_p mFile_f[6];
    for(int i=0; i<6; i++){
        mFile_f[i] = (func_p)mFile_f_[i];
    }
    // using Eigen::MatrixXd;
    // using Eigen::VectorXd;

    // Seprate inf and sup
    Matrix_t<double>im1_inf = im1.inf;
    Matrix_t<double>im1_sup = im1.sup;
    Matrix_t<double>im2_inf = im2.inf;
    Matrix_t<double>im2_sup = im2.sup;

    autodiff::ArrayXdual2nd tmp_im1_inf, tmp_im1_sup, tmp_im2_inf, tmp_im2_sup;
    tmp_im1_inf = autodiff::ArrayXdual2nd(im1_inf.rows());
    tmp_im1_sup = autodiff::ArrayXdual2nd(im1_sup.rows());
    tmp_im2_inf = autodiff::ArrayXdual2nd(im2_inf.rows());
    tmp_im2_sup = autodiff::ArrayXdual2nd(im2_sup.rows());
    for(int i=0; i<im1_inf.rows(); i++){
        tmp_im1_inf[i] = im1_inf(i, 0);
        tmp_im1_sup[i] = im1_sup(i, 0);
    }
    for(int i=0; i<im2_inf.rows(); i++){
        tmp_im2_inf[i] = im2_inf(i, 0);
        tmp_im2_sup[i] = im2_sup(i, 0);
    }

    // MatrixXd Jdyn_comb_inf = jacobian(this->mfile_, wrt(im1_inf, im2_inf), at(im1_inf, im2_inf), F);
    // MatrixXd Jdyn_comb_sup = jacobian(this->mfile_, wrt(im1_sup, im2_sup), at(im1_sup, im2_sup), F);
    // assert(Jdyn_comb_inf.rows() == Jdyn_comb_sup.rows());
    std::vector<IntervalMatrix> result = std::vector<IntervalMatrix>(im1_inf.rows());
    IntervalMatrix tmp_result;
    autodiff::dual2nd u; // the output scalar u = f(x) evaluated together with Hessian below
    autodiff::VectorXdual g; // gradient of f(x) evaluated together with Hessian below
    for(int k = 0; k<im1_inf.rows(); k++){
        tmp_result.inf = autodiff::hessian(mFile_f[k], autodiff::wrt(tmp_im1_inf, tmp_im2_inf), autodiff::at(tmp_im1_inf, tmp_im2_inf), u, g); 
        tmp_result.sup = autodiff::hessian(mFile_f[k], autodiff::wrt(tmp_im1_sup, tmp_im2_sup), autodiff::at(tmp_im1_sup, tmp_im2_sup), u, g); 
        result[k] = tmp_result;
    }
    return result;

}

}