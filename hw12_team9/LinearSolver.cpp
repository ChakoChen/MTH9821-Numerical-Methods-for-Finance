#include "LinearSolver.hpp"

vec forward_subst(const mat & L, const vec & b)
{
	int n = (int)L.rows();
	vec x(n);

	for (int i = 0; i < n; i++)
    {
		double sum = 0;
		for (int j = 0; j <= i - 1; j++)
			sum += L(i, j)*x(j);

		x(i) = (b(i)-sum)/L(i, i);
	}
	return x;
}

vec backward_subst(const mat & U, const vec & b)
{
    int n = (int)U.rows();
	vec x(n);
	for (int i=n-1; i>=0; i--)
    {
		double sum = 0;
		for (int j = i+1; j<n; j++)
            sum += U(i, j)*x(j);

        x(i) = (b(i) - sum) / U(i, i);
	}

	return x;
}

std::tuple <mat, mat> lu_no_pivoting(mat A)
{
	int n = (int)A.rows();
	mat L(n, n), U(n, n);
	L.setIdentity(); U.setIdentity();
	for (int i = 0; i < n; i++)
    {
		for (int j=i; j<n; j++)
        {
            U(i, j) = A(i, j);
			L(j, i) = A(j, i)/U(i, i);
        }

		for (int j=i+1; j<n; j++)
			for (int k=i+1; k<n; k++)
                A(j, k) = A(j, k)-L(j, i)*U(i, k);
	}

	return std::make_tuple(L, U);
}

std::tuple <permutation, mat, mat> lu_row_pivoting(mat A)
{
	int n = (int)A.rows();
	mat L(n, n), U(n, n);

	L.setIdentity(); U.setIdentity();

	//permutation matrix
	permutation P(n);

	P.setIdentity();

	for (int i = 0; i < n; i++) {

		int i_max = i;

		double temp = 0;

		permutation perm(n);

		perm.setIdentity();

		for (int j = i; j < n; j++) {

			if (abs(A(j, i)) > temp) {

				i_max = j;

				temp = abs(A(j, i));

			}
		}

		perm.applyTranspositionOnTheLeft(i, i_max);

		perm.applyThisOnTheLeft(A);

		perm.applyThisOnTheLeft(P);

		if (i > 0) {

			perm.applyThisOnTheLeft(L);

			perm.applyThisOnTheRight(L);

		}

		for (int j = i; j < n; j++) {

			U(i, j) = A(i, j);

			L(j, i) = A(j, i) / U(i, i);

		}
		//updating
		for (int j = i + 1; j < n; j++) {

			for (int k = i + 1; k < n; k++) {

				A(j, k) = A(j, k) - L(j, i) * U(i, k);

			}

		}

	}

	return std::make_tuple(P, L, U);

}


std::tuple<mat, mat> lu_tridiag(mat A) {

	int n = A.rows();

	mat L(n, n);

	mat U(n, n);

	L.setIdentity(); U.setIdentity();

	for (int i = 0; i < n - 1; i++) {

		L(i, i) = 1;

		L(i + 1, i) = A(i + 1, i) / A(i, i);

		U(i, i) = A(i, i);

		U(i, i + 1) = A(i, i + 1);

		A(i + 1, i + 1) = A(i + 1, i + 1) - L(i + 1, i) * U(i, i + 1);

	}

	L(n - 1, n - 1) = 1;

	U(n - 1, n - 1) = A(n - 1, n - 1);

	return std::make_tuple(L, U);
}

vec lu_solver(mat A, vec b) {

	int n = A.rows();

	mat L(n, n);

	mat U(n, n);

	vec y(n), x(n);

	std::tuple<mat, mat> decomp = lu_no_pivoting(A);

	L = std::get<0>(decomp);

	U = std::get<1>(decomp);

	//forward substitution
	y = forward_subst(L, b);

	//backward substitution
	x = backward_subst(U, y);

	return x;

}


vec lu_tridiag_solver(mat A, vec b) {

	int n = A.rows();

	mat L(n, n);

	mat U(n, n);

	vec y(n), x(n);

	std::tuple<mat, mat> decomp = lu_tridiag(A);

	L = std::get<0>(decomp);

	U = std::get<1>(decomp);

	//forward substitution
	y(0) = b(0) / L(0, 0);

	for (int j = 1; j < n; j++) {
		y(j) = (b(j) - L(j, j - 1) * y(j - 1)) / L(j, j);
	}

	//backward substitution
	x(n - 1) = y(n - 1) / U(n - 1, n - 1);

	for (int j = (n - 2); j >= 0; j--)

		x(j) = (y(j) - U(j, j + 1) * x(j + 1)) / U(j, j);

	return x;

}

mat cholesky(mat A) {

	int n = A.rows();

	mat U(n, n);

	U.setIdentity();

	for (int i = 0; i < n; i++) {

		U(i, i) = sqrt(A(i, i));

		for (int j = i + 1; j < n; j++) {

			U(i, j) = A(i, j) / U(i, i);

		}

		//updating
		for (int j = i + 1; j < n; j++) {

			for (int k = j; k < n; k++) {

				A(j, k) = A(j, k) - U(i, j) * U(i, k);

			}

		}

	}

	return U;

}


std::tuple<vec, uint32_t, double> gauss_seidel(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion) {

	int n = A.rows();

	//A = L_A + U_A + D_A
	mat L_A(n, n), U_A(n, n), D_A(n, n);

	L_A.setZero(), U_A.setZero(), D_A.setZero();

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			if (j < i) {

				L_A(i, j) = A(i, j);

			}

			else if (i == j) {

				D_A(i, j) = A(i, j);

			}

			else {

				U_A(i, j) = A(i, j);

			}

		}

	}

	vec x_old = x_0, x_new = x_0;

	vec r = b - (A * x_0);

	vec b_new = forward_subst(D_A + L_A, b);

	//number of iteration
	int n_iter = 0;

	double diff = 0, tol = 0;

	switch (criterion) {

	case residual:

		tol = tolerance * r.norm();

		diff = r.norm();

		while (diff > tol) {

			x_new = b_new - forward_subst(D_A + L_A, U_A * x_old);

			//update residue
			diff = (b - (A * x_new)).norm();

			x_old = x_new;

			n_iter += 1;

		}

	case consecutive:

		tol = tolerance;

		diff = 10 * tolerance;

		while (diff > tol) {

			x_new = b_new - forward_subst(D_A + L_A, U_A * x_old);

			//update diff
			diff = (x_new - x_old).norm();

			x_old = x_new;

			n_iter += 1;

		}

	}

	return std::make_tuple(x_new, n_iter, diff);

}

std::tuple<vec, uint32_t, double> gauss_seidel(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion) {

	int n = A.rows();

	vec x_0(n);

	x_0.setZero();

	return gauss_seidel(A, b, x_0, tolerance, criterion);

}


std::tuple<vec, uint32_t, double> jacobi(const mat & A, const vec & b, const vec & x_0, const double tolerance, const StoppingCriterion criterion) {

	int n = A.rows();

	//A = L_A + U_A + D_A
	mat L_A(n, n), U_A(n, n), D(n, n);

	L_A.setZero(), U_A.setZero(), D.setZero();

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			if (j < i) {

				L_A(i, j) = A(i, j);

			}

			else if (i == j) {

				D(i, j) = 1 / A(i, j);

			}

			else {

				U_A(i, j) = A(i, j);

			}

		}

	}

	vec x_old = x_0, x_new = x_0;

	vec r = b - (A * x_0);

	vec b_new = D * b;

	//number of iteration
	int n_iter = 0;

	double diff = 0, tol = 0;

	switch (criterion) {

	case residual:

		tol = tolerance * r.norm();

		diff = r.norm();

		while (diff > tol) {

			x_new = b_new - D * (L_A + U_A) * x_old;

			//update residue
			diff = (b - (A * x_new)).norm();

			x_old = x_new;

			n_iter += 1;

		}

	case consecutive:

		tol = tolerance;

		diff = 10 * tolerance;

		while (diff > tol) {

			x_new = b_new - D * (L_A + U_A) * x_old;

			//update diff
			diff = (x_new - x_old).norm();

			x_old = x_new;

			n_iter += 1;

		}

	}

	return std::make_tuple(x_new, n_iter, diff);

}

std::tuple<vec, uint32_t, double> jacobi(const mat & A, const vec & b, const double tolerance, const StoppingCriterion criterion) {

	int n = A.rows();

	vec x_0(n);

	x_0.setZero();

	return jacobi(A, b, x_0, tolerance, criterion);

}

std::tuple<vec, uint32_t, double> sor(const mat & A, const vec & b, const vec & x_0, const double tolerance, const double omega, const StoppingCriterion criterion) {

	int n = A.rows();

	//A = L_A + U_A + D_A
	mat L_A(n, n), U_A(n, n), D_A(n, n);

	L_A.setZero(), U_A.setZero(), D_A.setZero();

	for (int i = 0; i < n; i++) {

		for (int j = 0; j < n; j++) {

			if (j < i) {

				L_A(i, j) = A(i, j);

			}

			else if (i == j) {

				D_A(i, j) = A(i, j);

			}

			else {

				U_A(i, j) = A(i, j);

			}

		}

	}

	vec x_old = x_0, x_new = x_0;

	vec r = b - (A * x_0);

	vec b_new = omega * forward_subst(D_A + omega * L_A, b);

	//number of iteration
	int n_iter = 0;

	double diff = 0, tol = 0;

	switch (criterion) {

	case residual:

		tol = tolerance * r.norm();

		diff = r.norm();

		while (diff > tol) {

			x_new = b_new + forward_subst(D_A + (omega * L_A), (1 - omega) * (D_A * x_old) - (omega * U_A * x_old));

			//update residue
			diff = (b - (A * x_new)).norm();

			x_old = x_new;

			n_iter += 1;

		}

	case consecutive:

		tol = tolerance;

		diff = 10 * tolerance;

		while (diff > tol) {

			x_new = b_new + forward_subst(D_A + (omega * L_A), (1 - omega) * (D_A * x_old) - (omega * U_A * x_old));

			//update diff
			diff = (x_new - x_old).norm();

			x_old = x_new;

			n_iter += 1;

		}

	}

	return std::make_tuple(x_new, n_iter, diff);

}

std::tuple<vec, uint32_t, double> sor(const mat & A, const vec & b, const double tolerance, const double omega, const StoppingCriterion criterion) {

	int n = A.rows();

	vec x_0(n);

	x_0.setZero();

	return sor(A, b, x_0, tolerance, omega, criterion);

}


vec psor_tridiag(const mat & A, const vec & b, const vec & x_0,  const double tolerance, const double omega, const bool early_ex) {

	int n = A.rows();

	Eigen::VectorXd x_old = x_0;

	Eigen::VectorXd x_new(n);

	double diff = 10 * tolerance;

	while (diff > tolerance) {

		for (int j = 0; j < n; j++) {

			double temp1 = 0;

			double temp2 = 0;

			if (j >= 1) { temp1 += A(j, j - 1) * x_new(j - 1); }

			if (j < n - 1) { temp2 += A(j, j + 1) * x_old(j + 1); }

			x_new(j) = (1 - omega) * x_old(j) - omega * (temp1 + temp2) / A(j, j) + omega * b(j) / A(j, j);

			if (early_ex) { x_new(j) = std::max(x_new(j), x_0(j)); }

		}

		diff = (x_new - x_old).norm();

		x_old = x_new;

	}

	return x_new;
}

vec regressor_solver(mat A, vec b) {

	return lu_solver(A.transpose() * A, A.transpose() * b);

}
