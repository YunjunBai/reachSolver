#include "NonlinearSys.h"

namespace reachSolver{

typedef Vector_t<double> (*mFile_type)(Vector_t<double> vector1, Vector_t<double> vector2);
typedef double (**mFile_f_type)(Vector_t<double> vector1, Vector_t<double> vector2);

void jacobian(mFile_type mFile_, Vector_t<double> vector1, Vector_t<double> vector2, Matrix_t<double> matrix1, Matrix_t<double> matrix2);
std::vector<IntervalMatrix> hessian(mFile_f_type mFile_f_, IntervalMatrix im1, IntervalMatrix im2);

}