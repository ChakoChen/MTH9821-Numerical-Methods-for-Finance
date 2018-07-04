#include <cmath>
#include <functional>
#include <algorithm>
#include <tuple>
#include <Eigen>

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix <-1, -1, uint32_t> permutation;
enum StoppingCriterion { consecutive, residual };


vec forward_subst(const mat & L, const vec & b);

vec backward_subst(const mat & U, const vec & b);

std::tuple <mat, mat> lu_no_pivoting(mat A);

std::tuple <permutation, mat, mat> lu_row_pivoting(mat A);

std::tuple<mat, mat> lu_tridiag(mat A);

vec lu_solver(mat A, vec b);

vec lu_tridiag_solver(mat A, vec b);

mat cholesky(mat A);

std::tuple<vec, uint32_t, double> gauss_seidel(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion);

std::tuple<vec, uint32_t, double> gauss_seidel(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion);


std::tuple<vec, uint32_t, double> jacobi(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion);

std::tuple<vec, uint32_t, double> jacobi(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion);

std::tuple<vec, uint32_t, double> sor(const mat & A, const vec & b, const vec & x_0, const double tolerance, const double omega, const StoppingCriterion criterion);

std::tuple<vec, uint32_t, double> sor(const mat & A, const vec & b, const double tolerance, const double omega, const StoppingCriterion criterion);

vec psor_tridiag(const mat & A, const vec & b, const vec & x_0, const double tolerance, const double omega, const bool early_ex);

vec regressor_solver(mat A, vec b);
