#include "staticlib.h"
#include "TestMatrices.h"

#ifndef USE_ACCELERATE_SPARSE

#include "PardisoInterface.h"
#include "mkl.h"

bool test_sparselib() {

    using namespace Eigen;
    
    int num_threads = 8;

    int nx = 7;         
    int nu = 3;  
    int cardstates = 2; 

    SparseMatrix<double,EIGEN_STORAGE_ORDER> kkt = collocation_kkt_matrix(nx, nu, cardstates, 128);

    PardisoLDLT<SparseMatrix<double, EIGEN_STORAGE_ORDER>> kktsol;
    kktsol.m_ord = 2;
    kktsol.m_pivotstrat = 1;
    kktsol.m_pivotpert = 8;
    kktsol.m_matching = 1;
    kktsol.m_scaling = 0;
    kktsol.m_iterref = 0;
    kktsol.m_alg = 0;
    kktsol.m_msglvl = 0;
    kktsol.m_parsolve = 0;
    kktsol.setParams();

    double rtol = 1.0e-12; 

    std::vector<int> nthreads{ 1,2,3,4,5,6,7,8 };
    std::vector<int> nsegs{ 32,64,128,256 };

    for (auto threads : nthreads) {
        mkl_set_num_threads(threads);

        for (auto segs : nsegs) {

            SparseMatrix<double, EIGEN_STORAGE_ORDER> kkt = collocation_kkt_matrix(nx, nu, cardstates, segs);

            Eigen::VectorXd x(kkt.cols());
            Eigen::VectorXd b(kkt.cols());
            Eigen::VectorXd r(kkt.cols());

            b.setOnes();

            kktsol.compute(kkt);

            x = kktsol.solve(b);

            r = kkt.selfadjointView<Upper>() * x - b;

            if (r.norm() > rtol) return false;
        }
    }

    return true;
}

#else

#include "AccelerateInterface.h"
#include "Accelerate/Accelerate.h"
#include <iostream>
#include <cstdlib>

void accelerate_set_num_threads(int num_threads) {
    // Respect user-defined number of threads for Accelerate and set if unset
    const char* env_p = std::getenv("VECLIB_MAXIMUM_THREADS");
    if (!env_p)
        setenv("VECLIB_MAXIMUM_THREADS", std::to_string(num_threads).c_str(), 1);
}

bool test_sparselib() {

    using namespace Eigen;

    int num_threads = 8;

    int nx = 7;
    int nu = 3;
    int cardstates = 2;

    AccelerateLDLT<SparseMatrix<double, EIGEN_STORAGE_ORDER>, Upper> kktsol;
    kktsol.setOrder(SparseOrderMetis);
    kktsol.setIterativeRefinement(true);

    double rtol = 1.0e-12;

    std::vector<int> nthreads{ 1,2,3,4,5,6,7,8 };
    std::vector<int> nsegs{ 32,64,128,256 };
    for (auto threads : nthreads) {
        accelerate_set_num_threads(threads);

        for (auto segs : nsegs) {

            SparseMatrix<double, EIGEN_STORAGE_ORDER> kkt = 
                collocation_kkt_matrix(nx, nu, cardstates, segs);

            Eigen::VectorXd x(kkt.cols());
            Eigen::VectorXd b(kkt.cols());
            Eigen::VectorXd r(kkt.cols());

            b.setOnes();

            kktsol.compute(kkt);

            if (kktsol.info() != Success) {
                std::cout << "Decomposition failed" << std::endl;
                return false;
            }

            x = kktsol.solve(b);

            if (kktsol.info() != Success) {
                std::cout << "Solve failed" << std::endl;
                return false;
            }

            r = kkt.selfadjointView<Upper>() * x - b;

            std::cout << r.norm() << std::endl;
            if (r.norm() > rtol) return false;
        }
    }

    return true;
}

#endif