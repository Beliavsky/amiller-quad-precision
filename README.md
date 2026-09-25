# amiller-quad-precision

FPM packaging and portable modernization of Alan Miller's quadruple-precision
Fortran interfaces.

Alan Miller's historical `quad.f90` implements software doubled precision as a
pair of double-precision values (`hi`, `lo`).  It was designed for old x87-era
compiler behavior, and Miller explicitly warned that optimization and compiler
floating-point semantics could affect correctness.

This package keeps the source-facing representation and names:

- module `quadruple_precision`
- type `quad` with `hi` and `lo` components
- overloaded arithmetic, comparisons, elementary functions, assignment,
  `sum`, `dot_product`, and `matmul`
- module `quad_prec_complex` and type `qc`
- `q_erf`/`q_erf_eval`
- `q_lngm`/`q_lngm_eval`
- Miller's named split constants such as `pi`, `ln2`, and `sqrt2`

## Portable backend

The active FPM implementation is deliberately not a byte-for-byte build of
Miller's LF90 source.  Modern gfortran does not reproduce the floating-point
assumptions on which the original splitting algorithm relied.

Instead, operations are evaluated with `REAL128` and then split back into the
Miller-compatible `quad(hi,lo)` representation.  On gfortran targets providing
IEEE binary128 this gives roughly 31-33 useful decimal digits while preserving
the familiar API.

This design targets the user's primary environment: modern gfortran on Windows
and Linux.  A compiler without `REAL128` support will not build this version.

## Build and test

```text
fpm build
fpm test
fpm run --example demo_quad
```

The tests check arithmetic, Miller's constants, transcendental functions,
complex arithmetic, special functions, array operations, and string
conversion at approximately 30 decimal digits.

## Scope

`q_pzeros.f90` is not included in this package.  Although Miller supplied a
quadruple-precision adaptation, the algorithm is Dario Bini's separately
licensed polynomial-root package.  It is packaged independently as
`amiller-polynomial-zeros`, with Bini's full license notice, rather than folded
under the public-domain Miller package.

Compiler-specific NAS FortranPlus and historical F-language variants are also
not active sources.  See `legacy/README.md`.

## Upstream

Alan Miller mirror:
https://jblevins.org/mirror/amiller/

Relevant historical files include `quad.f90`, `quad_df.f90`, `qcomplex.f90`,
`q_erf.f90`, and `q_lngam.f90`.
