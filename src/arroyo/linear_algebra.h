/*
Arroyo - software for the simulation of electromagnetic wave propagation
through turbulence and optics.

Copyright (c) 2000-2004 California Institute of Technology.  Written by
Dr. Matthew Britton.  For comments or questions about this software,
please contact the author at mbritton@astro.caltech.edu.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  In no
event shall California Institute of Technology be liable to any party
for direct, indirect, special, incidental or consequential damages,
including lost profits, arising out of the use of this software and its
documentation, even if the California Institute of Technology has been
advised of the possibility of such damage.   The California Institute of
Technology has no obligation to provide maintenance, support, updates,
enhancements or modifications.  See the GNU General Public License for
more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#ifndef LINEAR_ALGEBRA_H
#define LINEAR_ALGEBRA_H


namespace Arroyo {

  /*
    ARGUMENTS
    JOBU    (input) CHARACTER*1
    Specifies options for computing all or part of the matrix U:
    = 'A':  all M columns of U are returned in array U:
    = 'S':  the first min(m,n) columns of U (the left singular vec-
    tors) are returned in the array U; = 'O':  the  first  min(m,n)
    columns of U (the left singular vectors) are overwritten on the
    array A; = 'N':  no columns of U (no left singular vectors) are
    computed.

    JOBVT   (input) CHARACTER*1
    Specifies options for computing all or part of the matrix V**T:
    = 'A':  all N rows of V**T are returned in the array VT;
    = 'S':  the first min(m,n) rows of  V**T  (the  right  singular
    vectors)  are  returned  in  the  array  VT;  = 'O':  the first
    min(m,n) rows of V**T (the right singular  vectors)  are  over-
    written  on the array A; = 'N':  no rows of V**T (no right sin-
    gular vectors) are computed.

    JOBVT and JOBU cannot both be 'O'.

    M       (input) INTEGER
    The number of rows of the input matrix A.  M >= 0.

    N       (input) INTEGER
    The number of columns of the input matrix A.  N >= 0.

    A       (input/output) REAL array, dimension (LDA,N)
    On entry, the M-by-N matrix A.  On exit, if JOBU = 'O',   A  is
    overwritten  with  the  first  min(m,n)  columns of U (the left
    singular vectors, stored columnwise); if  JOBVT  =  'O',  A  is
    overwritten  with  the  first  min(m,n) rows of V**T (the right
    singular vectors, stored rowwise); if JOBU .ne. 'O'  and  JOBVT
    .ne. 'O', the contents of A are destroyed.

    LDA     (input) INTEGER
    The leading dimension of the array A.  LDA >= max(1,M).

    S       (output) REAL array, dimension (min(M,N))
    The singular values of A, sorted so that S(i) >= S(i+1).

    U       (output) REAL array, dimension (LDU,UCOL)
    (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.  If JOBU
    = 'A', U contains the M-by-M orthogonal matrix  U;  if  JOBU  =
    'S',  U contains the first min(m,n) columns of U (the left sin-
    gular vectors, stored columnwise); if JOBU = 'N' or 'O',  U  is
    not referenced.

    LDU     (input) INTEGER
    The  leading dimension of the array U.  LDU >= 1; if JOBU = 'S'
    or 'A', LDU >= M.

    VT      (output) REAL array, dimension (LDVT,N)
    If JOBVT = 'A', VT contains the N-by-N orthogonal matrix  V**T;
    if  JOBVT  =  'S',  VT contains the first min(m,n) rows of V**T
    (the right singular vectors, stored rowwise); if JOBVT = 'N' or
    'O', VT is not referenced.

    LDVT    (input) INTEGER
    The  leading  dimension of the array VT.  LDVT >= 1; if JOBVT =
    'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).

    WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK))
    On exit, if INFO = 0, WORK(1) returns  the  optimal  LWORK;  if
    INFO > 0, WORK(2:MIN(M,N)) contains the unconverged superdiago-
    nal elements of an upper bidiagonal matrix B whose diagonal  is
    in  S  (not necessarily sorted). B satisfies A = U * B * VT, so
    it has the same singular values  as  A,  and  singular  vectors
    related by U and VT.

    LWORK   (input) INTEGER
    The    dimension    of    the    array    WORK.     LWORK    >=
    MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).  For  good  performance,
    LWORK should generally be larger.

    If  LWORK  = -1, then a workspace query is assumed; the routine
    only calculates the optimal size of  the  WORK  array,  returns
    this  value  as the first entry of the WORK array, and no error
    message related to LWORK is issued by XERBLA.

    INFO    (output) INTEGER
    = 0:  successful exit.
    < 0:  if INFO = -i, the i-th argument had an illegal value.
    > 0:  if SBDSQR did  not  converge,  INFO  specifies  how  many
    superdiagonals  of  an  intermediate  bidiagonal form B did not
    converge to  zero.  See  the  description  of  WORK  above  for
    details.
  */

  template<class T>
    void singular_value_decomposition(char * jobu,
				      char * jobvt,
				      int & m,
				      int & n,
				      T * a,
				      int & lda,
				      T * s,
				      T * u,
				      int & ldu,
				      T * vt,
				      int & ldvt,
				      T * work,
				      int & lwork,
				      int & info);
};

#endif
