/**
 * @file slate_eigensolvers.cpp
 * @brief C-linkage wrappers around SLATE 2025.05.28 eigensolvers.
 *
 * These four functions are designed to be called from Fortran via
 * ISO C binding (bind(C)).  They are intended as replacements for
 * ScaLAPACK eigensolver calls inside ABINIT's m_slk Fortran module,
 * specifically in:
 *   - subroutine compute_eigen_problem       (standard: A*Z = Z*diag(w))
 *   - subroutine compute_generalized_eigen_problem (generalised: A*Z = B*Z*diag(w))
 *
 * Functions provided
 * ------------------
 *   slate_zheev_c   complex*16 Hermitian standard eigenproblem
 *                   Replaces PZHEEVD / PZHEEVX
 *   slate_dsyev_c   real*8 symmetric standard eigenproblem
 *                   Replaces PDSYEVD / PDSYEVX
 *   slate_zhegv_c   complex*16 Hermitian generalised eigenproblem (itype=1)
 *                   Replaces PZHEGVD / PZHEGVX
 *   slate_dsygv_c   real*8 symmetric generalised eigenproblem (itype=1)
 *                   Replaces PDSYGVD / PDSYGVX
 *
 *
 * -----------------------------------------------------------------------
 * Design notes
 * -----------------------------------------------------------------------
 * - All matrix data resides on CPU host buffers throughout.  When
 *   use_gpu != 0, SLATE stages data to GPU devices internally for the
 *   compute kernels and returns results to host on exit.
 * - SLATE heev / hegv always compute all n eigenvalues and eigenvectors.
 *   The nev argument is accepted for API compatibility; when nev < n only
 *   the first nev entries of w and the first nev columns of Z are
 *   meaningful to the caller.
 * - The process grid is assumed to follow the row-major ordering used by
 *   ABINIT (BLACS_GRIDINIT 'R'), which matches SLATE's fromScaLAPACK
 *   convention: rank r -> row r/npcol, column r%npcol.
 * - Only itype = 1 (A x = lambda B x) is implemented for the generalised
 *   problems.  B is overwritten by its Cholesky factor on exit.
 * - The block size nb must satisfy mb == nb (square blocks), consistent
 *   with SLATE's requirement for eigensolver routines.
 * - Errors inside SLATE (std::exception derivatives) are caught, reported
 *   to stderr, and escalated via MPI_Abort.
 */

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/* This file should only be processed if SLATE was enabled at configuration */
#ifdef HAVE_LINALG_SLATE

#include <slate/slate.hh>   // slate::HermitianMatrix, heev, hegv, Options, ...
#include <blas.hh>           // blas::Uplo

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

// ---------------------------------------------------------------------------
// Module-private helpers (anonymous namespace)
// ---------------------------------------------------------------------------
namespace {

/**
 * @brief Build slate::Options selecting the desired compute target.
 *
 * @param use_gpu  Non-zero → offload compute to GPU devices.
 *                 Zero      → run entirely on CPU (HostTask).
 */
inline slate::Options make_opts(int use_gpu) noexcept
{
    slate::Options opts;
    opts[slate::Option::Target] = (use_gpu != 0)
        ? slate::Target::Devices
        : slate::Target::HostTask;
    return opts;
}

/**
 * @brief Print an error message to stderr, then abort all MPI processes.
 *
 * Declared [[noreturn]] so that callers inside catch blocks do not need
 * a dummy return value.
 */
[[noreturn]] void fatal(MPI_Comm comm, const std::string& msg)
{
    std::fprintf(stderr, "\n[slate_eigensolvers] FATAL: %s\n", msg.c_str());
    std::fflush(stderr);
    MPI_Abort(comm, EXIT_FAILURE);
    std::terminate();   // unreachable; silences compiler [[noreturn]] warnings
}

} // anonymous namespace

// ---------------------------------------------------------------------------
extern "C" {
// ---------------------------------------------------------------------------

/**
 * @brief Standard complex Hermitian eigenvalue problem:  A * Z = Z * diag(w)
 *
 * SLATE wraps the existing ScaLAPACK-distributed Fortran buffers via
 * fromScaLAPACK (no data copy) and calls slate::heev.
 *
 * @param[in]     n        Global matrix order.
 * @param[in]     nb       ScaLAPACK block size (mb == nb required).
 * @param[in]     nprow    Number of process-grid rows.
 * @param[in]     npcol    Number of process-grid columns.
 * @param[in]     comm_f   Fortran MPI communicator handle (c_int / integer(4)).
 * @param[in]     lda      Local leading dimension of @p a_data on this rank.
 * @param[in,out] a_data   Local buffer of matrix A (Hermitian, upper triangle).
 *                         Destroyed on exit.
 * @param[in]     ldz      Local leading dimension of @p z_data.
 * @param[out]    z_data   Local buffer receiving the distributed eigenvectors.
 * @param[out]    w        Length-n array of eigenvalues in ascending order.
 *                         Identical on every process after the call.
 * @param[in]     nev      Number of eigenpairs requested (0 = all).
 *                         SLATE always computes all n; read only the first
 *                         @p nev entries of @p w and columns of @p z_data.
 * @param[in]     use_gpu  Non-zero to offload computation to GPU.
 * @param[out]    info     0 on success.  Errors are currently fatal.
 */
void slate_zheev_c(
    int  n,    int  nb,
    int  nprow, int npcol,
    int  comm_f,
    int  lda,  void* a_data,
    int  ldz,  void* z_data,
    double* w,
    int  nev,
    int  use_gpu,
    int* info)
{
    using scalar_t = std::complex<double>;
    *info = 0;

    MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));

    try {
        // Wrap existing Fortran/ScaLAPACK host buffers.
        // SLATE operates directly on these arrays; no allocation or copy occurs.
        auto A = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(a_data),
            static_cast<int64_t>(lda),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto Z = slate::Matrix<scalar_t>::fromScaLAPACK(
            static_cast<int64_t>(n),
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(z_data),
            static_cast<int64_t>(ldz),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        std::vector<double> Lambda(n);

        slate::heev(A, Lambda, Z, make_opts(use_gpu));

        // Eigenvalues are already identical on every rank; just copy the vector.
        std::copy(Lambda.cbegin(), Lambda.cend(), w);
    }
    catch (const std::exception& ex) {
        fatal(comm, std::string("slate_zheev_c: ") + ex.what());
    }
}

// ---------------------------------------------------------------------------

/**
 * @brief Standard real symmetric eigenvalue problem:  A * Z = Z * diag(w)
 *
 * Identical interface to slate_zheev_c; @p a_data and @p z_data are double*.
 * Uses HermitianMatrix<double> (real Hermitian ≡ real symmetric).
 */
void slate_dsyev_c(
    int  n,    int  nb,
    int  nprow, int npcol,
    int  comm_f,
    int  lda,  void* a_data,
    int  ldz,  void* z_data,
    double* w,
    int  nev,
    int  use_gpu,
    int* info)
{
    using scalar_t = double;
    *info = 0;

    MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));

    try {
        auto A = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(a_data),
            static_cast<int64_t>(lda),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto Z = slate::Matrix<scalar_t>::fromScaLAPACK(
            static_cast<int64_t>(n),
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(z_data),
            static_cast<int64_t>(ldz),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        std::vector<double> Lambda(n);

        // slate::heev is templated; for scalar_t=double it calls the
        // symmetric (dsyev) path internally.
        slate::heev(A, Lambda, Z, make_opts(use_gpu));

        std::copy(Lambda.cbegin(), Lambda.cend(), w);
    }
    catch (const std::exception& ex) {
        fatal(comm, std::string("slate_dsyev_c: ") + ex.what());
    }
}

// ---------------------------------------------------------------------------

/**
 * @brief Generalised complex Hermitian eigenvalue problem (itype=1):
 *        A * Z = B * Z * diag(w)
 *
 * @param[in]     ldb      Local leading dimension of @p b_data on this rank.
 * @param[in,out] b_data   Local buffer of the Hermitian positive-definite
 *                         matrix B (upper triangle convention).
 *                         Overwritten by the Cholesky factor of B on exit.
 *
 * All other parameters have the same meaning as in slate_zheev_c.
 */
void slate_zhegv_c(
    int  n,    int  nb,
    int  nprow, int npcol,
    int  comm_f,
    int  lda,  void* a_data,
    int  ldb,  void* b_data,
    int  ldz,  void* z_data,
    double* w,
    int  nev,
    int  use_gpu,
    int* info)
{
    using scalar_t = std::complex<double>;
    *info = 0;

    MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));

    //std::cout << "nprow: " << nprow << " npcol: " << npcol << std::endl;
    //std::cout.flush();
    try {
        auto A = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(a_data),
            static_cast<int64_t>(lda),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto B = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(b_data),
            static_cast<int64_t>(ldb),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto Z = slate::Matrix<scalar_t>::fromScaLAPACK(
            static_cast<int64_t>(n),
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(z_data),
            static_cast<int64_t>(ldz),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        std::vector<double> Lambda(n);

        // itype = 1  →  A x = lambda B x
        // On exit: A is destroyed, B holds the Cholesky factor, Z holds eigenvectors.
        slate::hegv(static_cast<int64_t>(1), A, B, Lambda, Z, make_opts(use_gpu));

        std::copy(Lambda.cbegin(), Lambda.cend(), w);
    }
    catch (const std::exception& ex) {
        fatal(comm, std::string("slate_zhegv_c: ") + ex.what());
    }
}

// ---------------------------------------------------------------------------

/**
 * @brief Generalised real symmetric eigenvalue problem (itype=1):
 *        A * Z = B * Z * diag(w)
 *
 * Identical interface to slate_zhegv_c; @p a_data, @p b_data, and @p z_data
 * are double*.
 */
void slate_dsygv_c(
    int  n,    int  nb,
    int  nprow, int npcol,
    int  comm_f,
    int  lda,  void* a_data,
    int  ldb,  void* b_data,
    int  ldz,  void* z_data,
    double* w,
    int  nev,
    int  use_gpu,
    int* info)
{
    using scalar_t = double;
    *info = 0;

    MPI_Comm comm = MPI_Comm_f2c(static_cast<MPI_Fint>(comm_f));

    try {
        auto A = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(a_data),
            static_cast<int64_t>(lda),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto B = slate::HermitianMatrix<scalar_t>::fromScaLAPACK(
            blas::Uplo::Lower,
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(b_data),
            static_cast<int64_t>(ldb),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        auto Z = slate::Matrix<scalar_t>::fromScaLAPACK(
            static_cast<int64_t>(n),
            static_cast<int64_t>(n),
            static_cast<scalar_t*>(z_data),
            static_cast<int64_t>(ldz),
            static_cast<int64_t>(nb),
            nprow, npcol, comm);

        std::vector<double> Lambda(n);

        // slate::hegv is templated; scalar_t=double dispatches to the
        // symmetric (dsygv) path internally.
        slate::hegv(static_cast<int64_t>(1), A, B, Lambda, Z, make_opts(use_gpu));

        std::copy(Lambda.cbegin(), Lambda.cend(), w);
    }
    catch (const std::exception& ex) {
        fatal(comm, std::string("slate_dsygv_c: ") + ex.what());
    }
}

// ---------------------------------------------------------------------------
} // extern "C"
// ---------------------------------------------------------------------------
#endif
