#ifndef LINEAR_SOLVER
#define LINEAR_SOLVER

#include <iostream>
#include <Eigen/Dense>
#include <iomanip>
#include <vector>
using namespace std;
using namespace Eigen;

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
//typedef Eigen::PermutationMatrix<-1, -1, uint> permutation;


class LSolver {
public:
	vec forward_subst(const mat & L, const vec & b);
	vec backward_subst(const mat & U, const vec & b);

	tuple<mat, mat> lu_no_pivoting(mat A);
	tuple<mat, mat> tridiag_lu_no_pivoting(mat A);
	mat cholesky(mat A);
};

#endif // !LINEAR_SOLVER
