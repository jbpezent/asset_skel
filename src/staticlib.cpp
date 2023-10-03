#include "staticlib.h"
#include "TestMatrices.h"

// This is more than likely not how we want to handle MLK vs Accelerate
// but its simple and works for now...
#define USE_ACCELERATE

#ifndef USE_ACCELERATE
    #include "PardisoInterface.h"
    #include "mkl.h"

    bool test_sparselib() {

        using namespace Eigen;
        
        int num_threads = 8;

        int nx = 7;         
        int nu = 3;  
        int cardstates = 2; 


        SparseMatrix<double,RowMajor> kkt = collocation_kkt_matrix(nx, nu, cardstates, 128);


        

        PardisoLDLT<SparseMatrix<double, RowMajor>> kktsol;
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

                SparseMatrix<double, RowMajor> kkt = collocation_kkt_matrix(nx, nu, cardstates, segs);

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
    bool test_sparselib() {
        return true;
    }
#endif