ORRB:  Orthogonal Rank-structured Reduction to Banded form
==========================================================

The Fortran 2003 code works on various Givens-weight parameterizations
and is intended to accompany the paper *On the Stability of Orthogonal
Transformation of Rank Structured Matrices*.  Since I expect to
continue to use this code as a basis for further research into related
algorithms, I have put more effort into writing, packaging, and
testing than I would normally do for a proof of concept
implementation.  Nevertheless, the code is experimental and no
significant attempt has been made at optimization.  These disclaimers
apply most strongly to the decomposition routines for rank structured
matrices stored as general matrices, which are contained in
`./src/general`.  These decomposition routines are not part of the
paper and were included only because they are convenient to use as a
means for generating givens weight representations.  They represent an
experiment in applying QR updating techniques; they did not end up
being very efficient.

Build Process
=============

The code has been tested on a Debian GNU/Linux 7 system with
`gfortran` version 4.7.2.  The makefile is arranged as a single top
level makefile with includes in `./src/src.mk`, `./test/test.mk`, and
`./exp/exp.mk`.  There are a number of higher level targets:

|-----------------|--------------------------------------------|
| `all`           | build all modules                          |
| `all_tests`     | build all tests                            |
| `run_all_tests` | build and run all tests                    |
| `run_exp`       | build and run numerical experiments        |
| `view_notes`    | Convert this file to pdf and open a viewer |
|-----------------|--------------------------------------------|

There are also a number of more specific targets in `./src/src.mk`,
`./test/test.mk`, and `./exp/exp.mk`.

This document is in Markdown with limited use of LaTeX for some
equations.  The target `view_notes` will build a PDF on systems on
which both pandoc and LaTeX are available.

General Design and Naming Conventions
=====================================

The modules provide a number of data structures for representing real
double precision and complex double precision banded givens weight
decompositions of UB, BT, UBT, WB, BV, and WBV type.  Definitions of
these structures are in `./src/orth_band_types.f90`.  Such types
typically include allocatable components and the module in which each
type is defined includes special routines for allocation or
deallocation, e.g. `d_new_ub` in `mod_orth_band_types` defined in
`./src/types/orth_band_types.f90`.

Most of the computational routines have interfaces at three different
levels.  The lowest level subroutines minimize the use of derived types
and have minimal error handling.  Such routines have a prefix `f_` in
their names and are intended to be called by higher level Fortran
routines that validate input, detect and return errors, and possibly
allocate storage.  The lowest level routines are also intended to
provide a very basic interface that can be called from other
languages.  An example of such a routine is `f_d_qr_bv_to_ub` in
`mod_qr_factorization` defined in `./src/solve/qr_factorization.f90`.
Note that computational routines are tagged with `d_` or `z_`
depending on type.  There are generic interfaces for the routines that
allow these tags to be omitted.

At the middle level, most subroutines take already allocated derived
types and do some computation, returning the result in another
allocated derived type.  Such routines typically have `_to_` in their
names but no `_f_` prefix.  For example the call `call
d_qr_bv_to_ub(bv,ub,sw,error)` takes a matrix reresented by a BV
decomposition and stored in `bv`, computes a QR factorization using
the algorithm from the paper, and stores the resulting UB
decomposition in `ub`.  The parameter `sw`, which is used to return
the orthogonal factor Q, is a derived type for representing sweeps of
plane rotations and is defined in `./src/transforms/sweeps.f90`.
Functions at this level often overwrite their input, as the above call
to `d_qr_bv_to_ub` does to `bv`.  If needed, there is a generic
interface, `copy`, which can be used to explicitly copy givens weight
representations so that they can be saved prior to a destructive,
in-place computation.

At the top level are functions that take a derived typed, do some
computation to it, and allocate storage for the results.  These
functions allocate needed temporary storage as well and do not
overwrite their input.  Such functions are typically tagged with
`_of`.  For example the call `swub=d_qr_of(bv,error)` is a function
call that computes a QR factorization of the matrix represented by the
BV decomposition stored in `bv`.  The variable `swub` has type `d_qr`
which is a derived type that contains a sweep of rotations
representing Q and a UB decomposition representing R.  The BV
decomposition is not modified.

Error Handling
==============

Middle level and top level computational routines that implement error
handling use an error type `type(error_info)` defined in the module
`mod_error_id` in the file `./misc/error_id.f90`.  Each routine that
implements error handling is assigned an ID number.  A structure
`error` of `type(error_info)` includes a component `error%code` for the
error code, an array `error%routines` for storing ID numbers of
calling routines, an index `error%rix` into `error%routines`
indicating the position of the current routine, and a flag
`error%halt` that indicates whether a routine that is passed the
variable `error` should halt at an error or return an error code.

The components `error%routines` and `error%rix` serve as a simple call
stack so that the history of calls leading up to an error can be
printed out.  Each routine pushes its ID on the stack on entry and
pops it off the stack on successful exit.  If an error is already set
on entry into a routine, the routine immediately exits, preserving the
information in `error`.  Thus it is possible to get information about
the first error in a sequence of subroutine or function calls without
having to check the error parameter after each call.

The parameter `error` in each computational routine is optional.  If
it is absent, the routine will halt and print a message at the first
error.  If `error` is present and `error%halt` is `.true.` then after
encountering an error, the routine will print an error message and the
calling history that led to the error.  If `error%halt` is `.false.`
then instead of halting, `error%code` is set with an appropriate error
code and the routine exits.

Inventory and Description of Specific Modules
=============================================

The modules provided by the code are in ./src/ in the various
directories listed below.

`./src/types`
-------------

These modules contain derived types used to represent banded matrices
or Givens weight representations, along with routines for allocations,
deallocations, generation of random structures, and basic
manipulations of the structures.

* `mod_band_types` in `./src/types/band_types.f90`: Routines for
  working with compact representations of banded matrices.  A banded
  matrix $A=[a_{jk}]$ with $a_{jk}=0$ for $k-j>s$ and $k-j < -r$ is
  stored in one of two ways: In the first format, elements of $A$ are
  stored in a matrix $B$ with columns aligned with those of $A$ and
  $a_{j,k}=b_{j-k+s+1,k}$. Thus each row of $B$ corresponds to a
  diagonal of $A$ and $B$ can be chosen to be $(r+s+1)\times n$.
  Routines working with this format are given the tag `_bc`. In the
  second format, $B$ has rows aligned with those of $A$ and
  $a_{j,k}=b_{j,k-j+r+1}$.  In this case, each column of $B$
  corresponds to a diagonal of $A$ and $B$ can be chosen to be
  $n\times (r+s+1)$.  Routines working with this format are given the
  tag `_br`.  The module provides routines for getting and setting
  elements $a_{jk}$ of a matrix $A$ stored in this format; for
  applying applying rotations in a truncated manner corresponding to
  $j$-leading or $j$-trailing transformations; for printing band
  matrices in a readable format; and for converting between different
  storage formats.

* `mod_orth_band_types` in `./src/types/orth_band_types.f90`: This
  module contains derived types for UB, BV, WB, BT, UBT, and WBV
  decompositions.  The types include components to store a band
  matrix, information on the bandwidth and amount of available
  storage, sines and cosines for rotations, and information on the
  rows or columns on which the rotations act.  There are special
  allocation and deallocation routines for each type.
  
* `mod_assemble` in `./src/types/assemble.f90`: This module contains
  routines for assembling a Givens weight parameterized matrix into an
  unstructured representation in general $n\times n$ array.
  
* `mod_random` in `./src/types/random.f90`: Routines for generating
  random band matrices and Givens weight representations.  Included
  are routines both for fixed and variable bandwidth matrices.

`./src/misc`
------------

* `mod_prec` in `./src/misc/prec.f90`: Defines kinds `dp` and `int32`
  and the parameter `eps` representing the unit round-off.

* `mod_error_id` in `./src/misc/error_id.f90`: Defines custom error
  types and structures holding information about routines and their
  associated error codes.

* `mod_utility` in `./src/misc/utility.f90`.  Various utilities,
  including: array deallocation functions that take optional
  parameters; norms; routines for printing matrices; routines for
  generating random matrices; and in-place transposition.
  
* `mod_misc` in `./src/misc/misc.f90`: A module that exists
  solely to collect together the other modules in `./src/misc`.
  
`./src/orth`
------------

Routines related to QR decomposition and rank determination:

* `mod_triangular` in `./src/orth/triangular.f90`: Solve and multiply
  routines for triangular matrices, including routines implementing
  the back-substitution from the LINPACK condition estimator.
* `mod_cond_triangular` in `./src/orth/cond_triangular.f90`: Routines
  for condition estimation and computation of null vectors for
  triangular matrices.
* `gs` in `./src/orth/gs.f90`: Routines for updating modified
  Gram-Schmidt QR decompositions.
* `orth` in `./src/orth/orth.f90`: A module that exists solely to
  collect together the other modules in `./src/orth`.
  
`./src/general`
---------------

Routines for converting a general matrix stored in an array into
various Givens weight representations.  These routines use modified
Gram-Schmidt updating and are not very efficient.  Only the `ub` and
`bv` are implemented directly.  The others are implemented in terms of
the `ub` and `bv` decomposition using transposition.

* `mod_general_bt` in `./src/general/general_bt.f90`.
* `mod_general_bv` in `./src/general/general_bv.f90`.
* `mod_general_wb` in `./src/general/general_wb.f90`.
* `mod_general_wbv` in `./src/general/general_wbv.f90`.
* `mod_general_ub` in `./src/general/general_ub.f90`.
* `mod_general_bt` in `./src/general/general_bt.f90`.
* `mod_general_ubt` in `./src/general/general_ubt.f90`.
* `mod_general` in `./src/general/general.f90`: Collects all of the above.

`./src/convert`
---------------

Linear complexity routines for conversion between leading and trailing
Givens weight decompositions.

* `mod_convert_ub_to_bv` in `./src/convert/convert_ub_to_bv.f90`
* `mod_convert_bv_to_ub` in `./src/convert/convert_bv_to_ub.f90`
* `mod_convert_bt_to_wb` in `./src/convert/convert_bt_to_wb.f90`
* `mod_convert_wb_to_bt` in `./src/convert/convert_wb_to_bt.f90`
* `mod_convert_ubt_to_wbv` in `./src/convert/convert_ubt_to_wbv.f90`
* `mod_convert_wbv_to_ubt` in `./src/convert/convert_wbv_to_ubt.f90`
* `mod_convert` in `./src/convert/convert.f90`: Collects all of the above.

`./src/solve`
-------------

The stable, linear complexity system solver described in the paper.

* `mod_row_compress` in `./src/solve/row_compress.f90`: Routines that
  when given a rank structured matrix A parameterized by a UBT
  decomposition, compute the BV decomposition of a row compression A=QC
  where C has band structure in its lower triangular part and rank structure
  in its upper triangular part.  The module also defines derived
  types `d_rc` and `z_rc` for storing both Q and C in a row compression.
* `mod_qr_factorization` in `./src/solve/qr_factorization.f90`:
  Routines that compute a QR factorization of a row compressed matrix
  C, where C is represented by a BV decomposition and the computed R
  is represented by a UB decomposition.  The module also defines
  derived types `d_qr` and `z_qr` for storing both Q and R in a QR
  factorization.
* `mod_solve` in `./src/solve/solve.f90`: Various routines for fast
  backward and forward substitution for upper triangular matrices represented
  by BV or UB decompositions.  The backward substitution solves $Tx=b$ while
  the forward substitution solve $x^T T = b^T$.
* `mod_cond_orth_band` in `./src/solve/cond_orth_band.f90`: Routines
  for condition estimation for a triangular matrix represented as a UB
  decomposition.

`./src/transforms`
------------------

* `mod_rotation` in `./src/transforms/rotation.f90`: A simple derived
  type for rotations, routines for computing rotations that introduce
  zeros, and routines for applying the rotations to matrices and
  vectors.
* `mod_shift` in `./src/transforms/shift.f90`: Simple routines for applying shifts
  to vectors and matrices.
* `mod_sweeps` in `./src/transforms/sweeps.f90`: Derived types for
  sequences of sweeps of plane rotations that act on adjacent rows.
  The format is flexible with separate arrays for the cosines, sines,
  and the rows on which the rotations act.  See the module for more
  detail.  The primary use is in representing sequences of order r
  leading or trailing transformations used in the computation of a row
  compression or QR factorization.  There are also routines for
  applying such transformations to unstructured matrices (with or
  without transposition) and routines for generating random sweeps
  (used in test code).  The data structures include a flag to indicate
  whether they should be applied in transposed or non-transposed form
  and there is a high level interface `trp` that performs
  transposition by copying the structure and toggling the flag.
* `mod_transforms` in `./src/transforms/transforms.f90`: Collects all
  of the above.

An Example
==========

The test code in `./test/` can serve as a guide to how the code can be
called.  As an example, consider the following code fragment, which
generates a random right hand side, a matrix $A$ represented by a
random UBT decomposition, performs a row compression $A=QC$ and a
$QR$ factorization $C=\hat{Q}R$, and then solves the system $Ax=b$.

```
program example
    use mod_orrb
    implicit none
    !
    type(error_info) :: error
    integer(kind=int32) :: na, lbwa, ubwa
    real(kind=dp), dimension(:,:), allocatable :: a
    real(kind=dp), dimension(:), allocatable :: x, rhs, res, rhs0
    type(d_qr), allocatable :: swub
    type(d_rc), allocatable :: swbv
    type(d_ubt), allocatable :: ubta
    !
    ! initialize_errors sets up array mapping error codes to error messages.
    call initialize_errors
    na=1000; lbwa=20; ubwa=25
    rhs=d_random_vector(na)
    rhs0=rhs
    ubta=d_random_ubt(na,lbwa,ubwa,error=error)
    a = general(ubta,error)
    swbv=rc_of(ubta,error)
    rhs=trp(swbv%sw) * rhs
    swub=qr_of(swbv%bv,error)
    rhs=trp(swub%sw) * rhs
    x=solve(swub%ub, rhs, error)
    res=matmul(a,x)-rhs0
    print *, "Norm of the residual: ", norm2(res)
end program example
```
