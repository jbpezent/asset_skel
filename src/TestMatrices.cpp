#include "TestMatrices.h"

Eigen::SparseMatrix<double, Eigen::RowMajor> collocation_kkt_matrix(int nx, int nu, int cardstates, int numsegs)
{

	Eigen::SparseMatrix<double, Eigen::RowMajor> m;

	int numstates = (cardstates - 1) * numsegs + 1;
	int numVars = (nx + nu + 1) * numstates;
	int numCons = (cardstates - 1) * numsegs * nx;


	int hsize = cardstates * (nx + nu + 1);
	int jsize = (cardstates - 1) * nx;

	int size = numVars + numCons;


	m.resize(size, size);

	std::vector<Eigen::Triplet<double>> trips;

	for (int i = 0; i < size; i++) {
		if (i < numVars) {
			trips.push_back(Eigen::Triplet<double>(i, i, 10.0));
		}
		else {
			trips.push_back(Eigen::Triplet<double>(i, i, 0.00));
		}
	}

	Eigen::MatrixXd jrand(hsize, jsize);
	jrand.setZero();

	Eigen::MatrixXd hrand(hsize, hsize);
	hrand.setZero();

	for (int ns = 0; ns < numsegs; ns++) {

		int hstart = (nx + nu + 1) * (cardstates - 1) * ns;
		int jstart = numVars + jsize * ns;

		jrand.setRandom();
		hrand.setRandom();

		jrand = jrand.cwiseAbs();
		hrand = (hrand.cwiseAbs() + hrand.cwiseAbs().transpose()) / 10.0;

		for (int i = 0; i < hsize; i++) {
			for (int j = i; j < hsize; j++) {
				trips.push_back(Eigen::Triplet<double>(hstart + i, hstart + j, hrand(i, j)));
			}

			for (int j = 0; j < jsize; j++) {
				trips.push_back(Eigen::Triplet<double>(hstart + i, jstart + j, jrand(i, j)));
			}
		}

	}


	m.setFromTriplets(trips.begin(), trips.end());

	m.makeCompressed();

	return m;
}
