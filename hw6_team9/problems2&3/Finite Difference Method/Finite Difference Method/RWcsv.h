#ifndef RWcsv_HPP
#define RWcsv_HPP
#include <fstream>
#include <Eigen/Dense>

using namespace Eigen;

const IOFormat matrix_format(FullPrecision, 0, ", ", "\n", "", "", "", "");

//transform STL vector into eigen vector for output into csv file
template<typename t>
VectorXd VecTransform(std::vector<t> v)
{
	std::size_t sz = v.size();
	double* ptr = &v[0];
	Eigen::Map<Eigen::VectorXd> eigenvector(ptr, sz);
	return eigenvector;
}

//write into csv
template<typename M>
M load_csv(const std::string & path) {
	std::ifstream indata;
	indata.open(path);
	std::string line;
	std::vector<double> values;
	uint rows = 0;
	while (std::getline(indata, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		while (std::getline(lineStream, cell, ',')) {
			values.push_back(std::stod(cell));
		}
		++rows;
	}
	return Map<const Matrix<typename M::Scalar, Dynamic, Dynamic, RowMajor>>(values.data(), rows, values.size() / rows);
}

template<typename M>
static void save_csv(const M & A, const std::string & path) {
	std::ofstream csv(path);
	csv << A.format(matrix_format);
	csv.close();
}

#endif
