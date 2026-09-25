# Provenance

## Alan Miller archive

The source collection is mirrored at:

https://jblevins.org/mirror/amiller/

The mirror states that code written by Alan Miller is released into the public
domain, while code from other authors or sources may carry other restrictions.

## Historical quadruple-precision module

Alan Miller's `quad.f90` / `quad_df.f90` is a software doubled-precision
package using a pair of double-precision values.  Miller cites Seppo
Linnainmaa, "Software for doubled-precision floating-point computations",
ACM TOMS 7 (1981), 272-283, for the basic arithmetic algorithms.  Miller's
latest revision listed by the mirror is 18 September 2002.

The historical implementation was explicitly compiler-sensitive.  In
particular, `quad.f90` includes assumptions about extended x87 intermediates,
and `quad_df.f90` changes the splitting constant for other compilers.

## This FPM modernization

The active implementation preserves Miller's public-facing type and interfaces
but evaluates operations in Fortran `REAL128`, splitting results into the
`hi`/`lo` fields.  This removes dependency on x87 excess precision and
optimization behavior while retaining source compatibility for common uses.

Miller's published split constants are retained in the active module.

## Complex arithmetic

The `quad_prec_complex` API follows Miller's `qcomplex.f90`, whose mirror lists
Alan Miller as programmer and a latest revision of 11 November 1999.  The
active code has been rewritten around the portable backend rather than copied
verbatim.

## Special functions

Miller's `q_erf.f90` (latest F90 revision listed as 4 February 1998) and
`q_lngam.f90` (Fortran 90 version dated 21 August 1997) provide the historical
APIs.  The active modernization uses the Fortran intrinsic `ERF`, `ERFC`, and
`LOG_GAMMA` in REAL128 and converts results to the Miller `quad` representation.

## Excluded polynomial solver

Miller's `q_pzeros.f90` is an adaptation of Dario Andrea Bini's Aberth-method
polynomial-root software.  Bini's source carries its own permissive notice
requiring preservation of the complete notice.  It is intentionally excluded
from this public-domain-focused package and should be packaged separately.
