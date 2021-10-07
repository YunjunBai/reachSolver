#include "NonlinearSys.h"
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include "derivatives.h"

namespace reachSolver{

// template class NonlinearSys<double>;

// template <typename Number>
// void NonlinearSys<Number>::jacobian_(Vector_t<Number> vector1, Vector_t<Number> vector2, Matrix_t<Number> matrix1, Matrix_t<Number> matrix2){
void jacobian(mFile_type mFile_, Vector_t<double> vector1, Vector_t<double> vector2, Matrix_t<double> matrix1, Matrix_t<double> matrix2){

    autodiff::VectorXreal F;       // the output vector F = f(x, p, q) evaluated together with Jacobian below

    matrix1  = autodiff::jacobian(mFile_, autodiff::wrt(vector1), autodiff::at(vector1, vector2), F);       // evaluate the function and the Jacobian matrix Jx = dF/dx
    matrix2  = autodiff::jacobian(mFile_, autodiff::wrt(vector2), autodiff::at(vector1, vector2), F);       // evaluate the function and the Jacobian matrix Jx = dF/du

    // std::cout << "F = \n" << F << std::endl;     // print the evaluated output vector F
    // std::cout << "Jx = \n" << Jx << std::endl;   // print the evaluated Jacobian matrix dF/dx
    // std::cout << "Jp = \n" << Jp << std::endl;   // print the evaluated Jacobian matrix dF/dp
    // std::cout << "Jq = \n" << Jq << std::endl;   // print the evaluated Jacobian matrix dF/dq
    // std::cout << "Jqpx = \n" << Jqpx << std::endl; // print the evaluated Jacobian matrix [dF/dq, dF/dp, dF/dx]
}

}