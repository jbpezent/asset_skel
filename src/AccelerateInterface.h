#ifndef EIGEN_ACCELERATESUPPORT_H
#define EIGEN_ACCELERATESUPPORT_H

#include <Eigen/src/Core/util/DisableStupidWarnings.h>
#include <Accelerate/Accelerate.h>

#include <Eigen/Sparse>

#include <cmath>

namespace Eigen {

template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
class AccelerateImpl;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateLLT
 * \brief A direct Cholesky (LLT) factorization and solver based on Accelerate
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo_ additional information about the matrix structure. Default is Lower.
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateLLT
 */
template <typename MatrixType, int UpLo = Lower>
using AccelerateLLT = AccelerateImpl<MatrixType, UpLo | Symmetric, SparseFactorizationCholesky, true>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateLDLT
 * \brief The default Cholesky (LDLT) factorization and solver based on Accelerate
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo_ additional information about the matrix structure. Default is Lower.
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateLDLT
 */
template <typename MatrixType, int UpLo = Lower>
using AccelerateLDLT = AccelerateImpl<MatrixType, UpLo | Symmetric, SparseFactorizationLDLT, true>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateLDLTUnpivoted
 * \brief A direct Cholesky-like LDL^T factorization and solver based on Accelerate with only 1x1 pivots and no pivoting
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo_ additional information about the matrix structure. Default is Lower.
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateLDLTUnpivoted
 */
template <typename MatrixType, int UpLo = Lower>
using AccelerateLDLTUnpivoted = AccelerateImpl<MatrixType, UpLo | Symmetric, SparseFactorizationLDLTUnpivoted, true>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateLDLTSBK
 * \brief A direct Cholesky (LDLT) factorization and solver based on Accelerate with Supernode Bunch-Kaufman and static
 * pivoting
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo_ additional information about the matrix structure. Default is Lower.
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateLDLTSBK
 */
template <typename MatrixType, int UpLo = Lower>
using AccelerateLDLTSBK = AccelerateImpl<MatrixType, UpLo | Symmetric, SparseFactorizationLDLTSBK, true>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateLDLTTPP
 * \brief A direct Cholesky (LDLT) factorization and solver based on Accelerate with full threshold partial pivoting
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 * \tparam UpLo_ additional information about the matrix structure. Default is Lower.
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateLDLTTPP
 */
template <typename MatrixType, int UpLo = Lower>
using AccelerateLDLTTPP = AccelerateImpl<MatrixType, UpLo | Symmetric, SparseFactorizationLDLTTPP, true>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateQR
 * \brief A QR factorization and solver based on Accelerate
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateQR
 */
template <typename MatrixType>
using AccelerateQR = AccelerateImpl<MatrixType, 0, SparseFactorizationQR, false>;

/** \ingroup AccelerateSupport_Module
 * \typedef AccelerateCholeskyAtA
 * \brief A QR factorization and solver based on Accelerate without storing Q (equivalent to A^TA = R^T R)
 *
 * \warning Only single and double precision real scalar types are supported by Accelerate
 *
 * \tparam MatrixType_ the type of the sparse matrix A, it must be a SparseMatrix<>
 *
 * \sa \ref TutorialSparseSolverConcept, class AccelerateCholeskyAtA
 */
template <typename MatrixType>
using AccelerateCholeskyAtA = AccelerateImpl<MatrixType, 0, SparseFactorizationCholeskyAtA, false>;

namespace internal {
template <typename T>
struct AccelFactorizationDeleter {
  void operator()(T* sym) {
    if (sym) {
      SparseCleanup(*sym);
      delete sym;
      sym = nullptr;
    }
  }
};

template <typename DenseVecT, typename DenseMatT, typename SparseMatT, typename NumFactT>
struct SparseTypesTraitBase {
  typedef DenseVecT AccelDenseVector;
  typedef DenseMatT AccelDenseMatrix;
  typedef SparseMatT AccelSparseMatrix;

  typedef SparseOpaqueSymbolicFactorization SymbolicFactorization;
  typedef NumFactT NumericFactorization;

  typedef AccelFactorizationDeleter<SymbolicFactorization> SymbolicFactorizationDeleter;
  typedef AccelFactorizationDeleter<NumericFactorization> NumericFactorizationDeleter;
};

template <typename Scalar>
struct SparseTypesTrait {};

template <>
struct SparseTypesTrait<double> : SparseTypesTraitBase<DenseVector_Double, DenseMatrix_Double, SparseMatrix_Double,
                                                       SparseOpaqueFactorization_Double> {};

template <>
struct SparseTypesTrait<float>
    : SparseTypesTraitBase<DenseVector_Float, DenseMatrix_Float, SparseMatrix_Float, SparseOpaqueFactorization_Float> {
};

// Taken from https://github.com/ceres-solver/ceres-solver/blob/master/internal/ceres/accelerate_sparse.cc
void* resizeForAccelerateAlignment(const size_t required_size, 
                                         std::vector<uint8_t>* mem_vec) 
{
  // As per the Accelerate documentation, all workspace memory passed to the
  // sparse solver functions must be 16-byte aligned.
  constexpr int kAccelerateRequiredAlignment = 16; 
  // Although malloc() on macOS should always be 16-byte aligned, it is unclear
  // if this holds for new(), or on other Apple OSs (phoneOS, watchOS etc).
  // As such we assume it is not and use std::align() to create a (potentially
  // offset) 16-byte aligned sub-buffer of the specified size within workspace.
  mem_vec->resize(required_size + kAccelerateRequiredAlignment);
  size_t size_from_aligned_start = mem_vec->size();
  void* aligned_solve_workspace_start =
      reinterpret_cast<void*>(mem_vec->data());
  aligned_solve_workspace_start = std::align(kAccelerateRequiredAlignment,
                                             required_size,
                                             aligned_solve_workspace_start,
                                             size_from_aligned_start);
  return aligned_solve_workspace_start;
}

}  // end namespace internal

template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
class AccelerateImpl : public SparseSolverBase<AccelerateImpl<MatrixType_, UpLo_, Solver_, EnforceSquare_> > {
 protected:
  using Base = SparseSolverBase<AccelerateImpl>;
  using Base::derived;
  using Base::m_isInitialized;

 public:
  using Base::_solve_impl;

  typedef MatrixType_ MatrixType;
  typedef typename MatrixType::Scalar Scalar;
  typedef typename MatrixType::StorageIndex StorageIndex;
  enum { ColsAtCompileTime = Dynamic, MaxColsAtCompileTime = Dynamic };
  enum { UpLo = UpLo_ };

  using AccelDenseVector = typename internal::SparseTypesTrait<Scalar>::AccelDenseVector;
  using AccelDenseMatrix = typename internal::SparseTypesTrait<Scalar>::AccelDenseMatrix;
  using AccelSparseMatrix = typename internal::SparseTypesTrait<Scalar>::AccelSparseMatrix;
  using SymbolicFactorization = typename internal::SparseTypesTrait<Scalar>::SymbolicFactorization;
  using NumericFactorization = typename internal::SparseTypesTrait<Scalar>::NumericFactorization;
  using SymbolicFactorizationDeleter = typename internal::SparseTypesTrait<Scalar>::SymbolicFactorizationDeleter;
  using NumericFactorizationDeleter = typename internal::SparseTypesTrait<Scalar>::NumericFactorizationDeleter;

  AccelerateImpl() {
    m_isInitialized = false;

    auto check_flag_set = [](int value, int flag) { return ((value & flag) == flag); };

    if (check_flag_set(UpLo_, Symmetric)) {
      m_sparseKind = SparseSymmetric;
      m_triType = (UpLo_ & Lower) ? SparseLowerTriangle : SparseUpperTriangle;
    } else if (check_flag_set(UpLo_, UnitLower)) {
      m_sparseKind = SparseUnitTriangular;
      m_triType = SparseLowerTriangle;
    } else if (check_flag_set(UpLo_, UnitUpper)) {
      m_sparseKind = SparseUnitTriangular;
      m_triType = SparseUpperTriangle;
    } else if (check_flag_set(UpLo_, StrictlyLower)) {
      m_sparseKind = SparseTriangular;
      m_triType = SparseLowerTriangle;
    } else if (check_flag_set(UpLo_, StrictlyUpper)) {
      m_sparseKind = SparseTriangular;
      m_triType = SparseUpperTriangle;
    } else if (check_flag_set(UpLo_, Lower)) {
      m_sparseKind = SparseTriangular;
      m_triType = SparseLowerTriangle;
    } else if (check_flag_set(UpLo_, Upper)) {
      m_sparseKind = SparseTriangular;
      m_triType = SparseUpperTriangle;
    } else {
      m_sparseKind = SparseOrdinary;
      m_triType = (UpLo_ & Lower) ? SparseLowerTriangle : SparseUpperTriangle;
    }

    m_order = SparseOrderDefault;
    //m_order = SparseOrderMetis;
    m_doIterativeRefinement = false;
    m_iterativeRefinementIterations = 2;
  }

  explicit AccelerateImpl(const MatrixType& matrix) : AccelerateImpl() { compute(matrix); }

  ~AccelerateImpl() {}

  inline Index cols() const { return m_nCols; }
  inline Index rows() const { return m_nRows; }

  ComputationInfo info() const {
    eigen_assert(m_isInitialized && "Decomposition is not initialized.");
    return m_info;
  }

  void analyzePattern(const MatrixType& matrix);

  void factorize(const MatrixType& matrix);

  void compute(const MatrixType& matrix);

  template <typename Rhs, typename Dest>
  void _solve_impl(const MatrixBase<Rhs>& b, MatrixBase<Dest>& dest) const;

  /** Sets the ordering algorithm to use. */
  void setOrder(SparseOrder_t order) { m_order = order; }

  void setIterativeRefinement(bool iterativeRefinement) { 
    m_doIterativeRefinement = iterativeRefinement; 
  }

  void setIterativeRefinementIterations(int iterations) {
    eigen_assert(iterations >= 0 && "Number of iterations must be non-negative.");
    m_iterativeRefinementIterations = iterations;
  }

 private:
  template <typename T>
  void buildAccelSparseMatrix(const SparseMatrix<T>& a) {
    const Index nColumnsStarts = a.cols() + 1;

    m_columnStarts.resize(nColumnsStarts);

    for (Index i = 0; i < nColumnsStarts; i++) m_columnStarts[i] = a.outerIndexPtr()[i];

    SparseAttributes_t attributes{};
    attributes.transpose = false;
    attributes.triangle = m_triType;
    attributes.kind = m_sparseKind;

    SparseMatrixStructure structure{};
    structure.attributes = attributes;
    structure.rowCount = static_cast<int>(a.rows());
    structure.columnCount = static_cast<int>(a.cols());
    structure.blockSize = 1;
    structure.columnStarts = m_columnStarts.data();
    structure.rowIndices = const_cast<int*>(a.innerIndexPtr());

    m_matrix.structure = structure;
    m_matrix.data = const_cast<T*>(a.valuePtr());
  }

  void doAnalysis() {
    m_numericFactorization.reset(nullptr);

    SparseSymbolicFactorOptions fopts{};
    fopts.control = SparseDefaultControl;
    fopts.orderMethod = m_order;
    fopts.order = nullptr;
    fopts.ignoreRowsAndColumns = nullptr;
    fopts.malloc = malloc;
    fopts.free = free;
    fopts.reportError = nullptr;

    m_symbolicFactorization.reset(new SymbolicFactorization(SparseFactor(Solver_, m_matrix.structure, fopts)));

    SparseStatus_t status = m_symbolicFactorization->status;

    updateInfoStatus(status);

    if (status != SparseStatusOK) m_symbolicFactorization.reset(nullptr);
  }

  void doFactorization() {
    SparseStatus_t status = SparseStatusReleased;

    if (m_symbolicFactorization) {

      SparseNumericFactorOptions nopts{};
      nopts.control = SparseDefaultControl;
      nopts.scalingMethod = SparseScalingDefault;
      nopts.scaling = nullptr;
      // Default values set by Apple
      nopts.pivotTolerance = 0.01;                   // Recommended value for difficult matrices in double
      nopts.zeroTolerance = 1e-4 * __DBL_EPSILON__;  // "A few" orders of magnitude below epsilon.

      // Get factor and workspace size
      const int factorSize = 
        std::is_same<Scalar, double>::value 
            ? m_symbolicFactorization->factorSize_Double 
            : m_symbolicFactorization->factorSize_Float;
      const int workspaceSize =
        std::is_same<Scalar, double>::value 
            ? m_symbolicFactorization->workspaceSize_Double 
            : m_symbolicFactorization->workspaceSize_Float;

      m_numericFactorization.reset(new NumericFactorization(SparseFactor(
        *m_symbolicFactorization, m_matrix, nopts, 
        internal::resizeForAccelerateAlignment(factorSize, &m_factorStorage), 
        internal::resizeForAccelerateAlignment(workspaceSize, &m_workspace))));

      status = m_numericFactorization->status;

      if (status != SparseStatusOK) m_numericFactorization.reset(nullptr);
    }

    updateInfoStatus(status);
  }

 protected:
  void updateInfoStatus(SparseStatus_t status) const {
    switch (status) {
      case SparseStatusOK:
        m_info = Success;
        break;
      case SparseFactorizationFailed:
      case SparseMatrixIsSingular:
        m_info = NumericalIssue;
        break;
      case SparseInternalError:
      case SparseParameterError:
      case SparseStatusReleased:
      default:
        m_info = InvalidInput;
        break;
    }
  }

  std::vector<long> m_columnStarts;
  mutable AccelSparseMatrix m_matrix;
  mutable ComputationInfo m_info;
  mutable std::vector<uint8_t> m_factorStorage;
  mutable std::vector<uint8_t> m_workspace;
  Index m_nRows, m_nCols;
  std::unique_ptr<SymbolicFactorization, SymbolicFactorizationDeleter> m_symbolicFactorization;
  std::unique_ptr<NumericFactorization, NumericFactorizationDeleter> m_numericFactorization;
  SparseKind_t m_sparseKind;
  SparseTriangle_t m_triType;
  SparseOrder_t m_order;
  bool m_doIterativeRefinement;
  int m_iterativeRefinementIterations;
};

/** Computes the symbolic and numeric decomposition of matrix \a a */
template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
void AccelerateImpl<MatrixType_, UpLo_, Solver_, EnforceSquare_>::compute(const MatrixType& a) {
  if (EnforceSquare_) eigen_assert(a.rows() == a.cols());

  m_nRows = a.rows();
  m_nCols = a.cols();

  buildAccelSparseMatrix(a);

  doAnalysis();

  if (m_symbolicFactorization) doFactorization();

  m_isInitialized = true;
}

/** Performs a symbolic decomposition on the sparsity pattern of matrix \a a.
 *
 * This function is particularly useful when solving for several problems having the same structure.
 *
 * \sa factorize()
 */
template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
void AccelerateImpl<MatrixType_, UpLo_, Solver_, EnforceSquare_>::analyzePattern(const MatrixType& a) {
  if (EnforceSquare_) eigen_assert(a.rows() == a.cols());

  m_nRows = a.rows();
  m_nCols = a.cols();

  std::vector<long> columnStarts;
  buildAccelSparseMatrix(a, columnStarts);

  doAnalysis();

  m_isInitialized = true;
}

/** Performs a numeric decomposition of matrix \a a.
 *
 * The given matrix must have the same sparsity pattern as the matrix on which the symbolic decomposition has been
 * performed.
 *
 * \sa analyzePattern()
 */
template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
void AccelerateImpl<MatrixType_, UpLo_, Solver_, EnforceSquare_>::factorize(const MatrixType& a) {
  eigen_assert(m_symbolicFactorization && "You must first call analyzePattern()");
  eigen_assert(m_nRows == a.rows() && m_nCols == a.cols());

  if (EnforceSquare_) eigen_assert(a.rows() == a.cols());

  std::vector<long> columnStarts;

  buildAccelSparseMatrix(a, columnStarts);

  doFactorization();
}

template <typename MatrixType_, int UpLo_, SparseFactorization_t Solver_, bool EnforceSquare_>
template <typename Rhs, typename Dest>
void AccelerateImpl<MatrixType_, UpLo_, Solver_, EnforceSquare_>::_solve_impl(const MatrixBase<Rhs>& b,
                                                                              MatrixBase<Dest>& x) const {
  if (!m_numericFactorization) {
    m_info = InvalidInput;
    return;
  }

  eigen_assert(m_nRows == b.rows());
  eigen_assert(((b.cols() == 1) || b.outerStride() == b.rows()));

  SparseStatus_t status = SparseStatusOK;

  Scalar* b_ptr = const_cast<Scalar*>(b.derived().data());
  Scalar* x_ptr = const_cast<Scalar*>(x.derived().data());

  AccelDenseMatrix xmat{};
  xmat.attributes = SparseAttributes_t();
  xmat.columnCount = static_cast<int>(x.cols());
  xmat.rowCount = static_cast<int>(x.rows());
  xmat.columnStride = xmat.rowCount;
  xmat.data = x_ptr;

  AccelDenseMatrix bmat{};
  bmat.attributes = SparseAttributes_t();
  bmat.columnCount = static_cast<int>(b.cols());
  bmat.rowCount = static_cast<int>(b.rows());
  bmat.columnStride = bmat.rowCount;
  bmat.data = b_ptr;

  const int nrhs = (bmat.attributes.transpose) ? bmat.rowCount : bmat.columnCount;
  const int workspaceSize = m_numericFactorization->solveWorkspaceRequiredStatic + 
    nrhs*m_numericFactorization->solveWorkspaceRequiredPerRHS;

  void* ws = internal::resizeForAccelerateAlignment(workspaceSize, &m_workspace);
  assert(ws != nullptr && "Accelerate workspace alignment failed");

  SparseSolve(*m_numericFactorization, bmat, xmat, ws);

  updateInfoStatus(status);

  if (m_doIterativeRefinement)
  {
    auto n = vDSP_Length(x.rows() * x.cols());
    auto r_mem = std::vector<Scalar>(x.rows() * x.cols(), Scalar(0));

    AccelDenseMatrix ref_mat{};
    ref_mat.attributes = SparseAttributes_t();
    ref_mat.columnCount = static_cast<int>(x.cols());
    ref_mat.rowCount = static_cast<int>(x.rows());
    ref_mat.columnStride = ref_mat.rowCount;
    ref_mat.data = r_mem.data();

    for (int i = 0; i < m_iterativeRefinementIterations; ++i) {
        // Calculate residual and store in ref_mat
        vDSP_vnegD(
            bmat.data, 1,
            ref_mat.data, 1, n
        );
        SparseMultiplyAdd(m_matrix, xmat, ref_mat);

        // Solve for correction and store in ref_mat
        SparseSolve(*m_numericFactorization, ref_mat, ws);

        // vDSP operation that calculates x -= correction
        vDSP_vsubD(
            ref_mat.data, 1,
            xmat.data, 1,
            xmat.data, 1,
            n
        );
    }
  }
}

}  // end namespace Eigen

#endif  // EIGEN_ACCELERATESUPPORT_H

