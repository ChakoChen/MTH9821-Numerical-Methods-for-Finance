#include "LinearSolver.hpp"
using namespace std;

vec LSolver::forward_subst(const mat & L, const vec & b) {
	/* Solve x for Lx = b;
	L: nonsingular lower triangular matrix of size n
	b: column vector of size n
	*/
	long n = L.rows();
	vec x(n);
	x.setZero();
	double sum;

	if (L.rows() != L.cols() || L.cols() != b.size()){
		cout << "forward_subst: dimension mismatch\n";
		return x;
	}

	x(0) = b(0) / L(0, 0);
	for (long j = 1; j <= n - 1; j++){
		sum = 0.;
		for (long k = 0; k<j; k++){
			sum += L(j, k) * x(k);
		}
		x(j) = (b(j) - sum) / L(j, j);
	}
	return x;
}

vec LSolver::backward_subst(const mat & U, const vec & b) {
	/* Solve x for Ux = b;
	U: nonsingular upper triangular matrix of size n
	b: column vector of size n
	*/

	long n = U.rows();
	vec x(n);
	x.setZero();
	double sum;

	if (U.rows() != U.cols() || U.cols() != b.size()){
		cout << "backward_subst: dimension mismatch\n";
		return x;
	}

	x(n - 1) = b(n - 1) / U(n - 1, n - 1);
	for (long j = n - 2; j >= 0; j--){
		sum = 0.;
		for (long k = j + 1; k <= n - 1; k++){
			sum += U(j, k)*x(k);
		}
		x(j) = (b(j) - sum) / U(j, j);
	}

	return x;
}

tuple<mat, mat> LSolver::lu_no_pivoting(mat A) { 
  // [L, U] = lu_no_pivoting(A)
  // A: nonsingular matrix with size n
  // L: lower triangular matrix with entries 1 on main diagonal
  // U: upper triangular matrix
  // A = LU

	long n = A.rows();
	mat L(n, n), U(n, n);

	L.setZero(); U.setZero();
	for (int i = 0; i <= n - 2; i++) {
		for (int k = i; k <= n - 1; k++) {
			U(i, k) = A(i, k);
			L(k, i) = A(k, i) / U(i, i);
		}
		for (int j = i + 1; j <= n - 1; j++) {
			for (int k = i + 1; k <= n - 1; k++) {
				A(j, k) -= L(j, i)*U(i, k);
			}
		}
	}
	L(n - 1, n - 1) = 1; U(n - 1, n - 1) = A(n - 1, n - 1);

	return make_tuple(L, U);
}

tuple<mat, mat> LSolver::tridiag_lu_no_pivoting(mat A) {
	// [L, U] = tridiag_lu_no_pivoting(A)
	// A: nonsingular tridiagonal matrix with size n
	// L: lower triangular matrix with entries 1 on main diagonal
	// U: upper triangular matrix
	// A = LU

	long n = A.rows();
	mat L(n, n), U(n, n);

	L.setZero(); U.setZero();
	for (int i = 0; i < n - 1; i++) {
		L(i, i) = 1;
		L(i + 1, i) = A(i + 1, i) / A(i, i);
		U(i, i) = A(i, i);
		U(i, i + 1) = A(i, i + 1);
		A(i + 1, i + 1) = A(i + 1, i + 1) - L(i + 1, i) * U(i, i + 1);
	}
	L(n - 1, n - 1) = 1;
	U(n - 1, n - 1) = A(n - 1, n - 1);

	return make_tuple(L, U);
}

mat LSolver::cholesky(mat A){ 
	// U = cholesky(A)
	// A: symmetric positive definite matrix of size n
	// U: upper triangular matrix s.t. (U^T)U=A

	long i, j, k;
	long n = A.rows();
	mat U(n, n);

	U.setZero();
	for (i = 0; i <= n - 2; i++) {
		U(i, i) = sqrt(A(i, i));
		for (k = i + 1; k <= n - 1; k++) {
			U(i, k) = A(i, k) / U(i, i);
		}
		for (j = i + 1; j <= n - 1; j++) {
			for (k = j; k <= n - 1; k++) {
				A(j, k) -= U(i, j)*U(i, k);
			}
		}
	}
	U(n - 1, n - 1) = sqrt(A(n - 1, n - 1));

	return U;
}