## ARPACK SOFTWARE

ARPACK is a collection of Fortran77 subroutines designed to solve large scale eigenvalue problems.

The package is designed to compute a few eigenvalues and corresponding eigenvectors of a general n by n matrix A. It is most appropriate for large sparse or structured matrices A where structured means that a matrix-vector product w <- Av requires order n rather than the usual order n2  floating point operations. This software is based upon an algorithmic variant of the Arnoldi process called the Implicitly Restarted Arnoldi Method (IRAM). When the matrix A is symmetric it reduces to a variant of the Lanczos process called the Implicitly Restarted Lanczos Method (IRLM). These variants may be viewed as a synthesis of the Arnoldi/Lanczos process with the Implicitly Shifted QR technique that is suitable for large scale problems. For many standard problems, a matrix factorization is not required. Only the action of the matrix on a vector is needed.

ARPACK software is capable of solving large scale symmetric, nonsymmetric, and generalized eigenproblems from significant application areas. The software is designed to compute a few (k) eigenvalues with user specified features such as those of largest real part or largest magnitude. Storage requirements are on the order of n*k locations. No auxiliary storage is required. A set of Schur basis vectors for the desired k-dimensional eigen-space is computed which is numerically orthogonal to working precision. Numerically accurate eigenvectors are available on request.

### Important Features:

-   Reverse Communication Interface.
-   Single and Double Precision Real Arithmetic Versions for Symmetric, Non-symmetric,
-   Standard or Generalized Problems.
-   Single and Double Precision Complex Arithmetic Versions for Standard or Generalized Problems.
-   Routines for Banded Matrices - Standard or Generalized Problems.
-   Routines for The Singular Value Decomposition.
-   Example driver routines that may be used as templates to implement numerous Shift-Invert
-   strategies for all problem types, data types and precision.

***
[https://www.caam.rice.edu/software/ARPACK/](https://www.caam.rice.edu/software/ARPACK/)
