# Krylov.jl: A Julia basket of hand-picked Krylov methods

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3548984.svg)](https://doi.org/10.5281/zenodo.3548984)

| **Documentation** | **Travis, AppVeyor and Cirrus build statuses** | **Coverage** |
|:-----------------:|:----------------------------------------------:|:------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSmoothOptimizers.github.io/Krylov.jl/dev) | [![Build Status](https://travis-ci.org/JuliaSmoothOptimizers/Krylov.jl.svg?branch=master)](https://travis-ci.org/JuliaSmoothOptimizers/Krylov.jl) [![Build status](https://ci.appveyor.com/api/projects/status/3xt558lune9f5r2v?svg=true)](https://ci.appveyor.com/project/dpo/krylov-jl) [![Build Status](https://api.cirrus-ci.com/github/JuliaSmoothOptimizers/Krylov.jl.svg)](https://cirrus-ci.com/github/JuliaSmoothOptimizers/Krylov.jl) | [![Coverage Status](https://coveralls.io/repos/github/JuliaSmoothOptimizers/Krylov.jl/badge.svg?branch=master)](https://coveralls.io/github/JuliaSmoothOptimizers/Krylov.jl?branch=master) [![codecov.io](https://codecov.io/github/JuliaSmoothOptimizers/Krylov.jl/coverage.svg?branch=master)](https://codecov.io/github/JuliaSmoothOptimizers/Krylov.jl?branch=master) |

## Purpose

This package provides implementations of certain of the most useful Krylov method for linear systems with special emphasis on methods for square systems, linear least-squares problems, linear least-norm problems, adjoint systems, saddle-point systems and symmetric quasi-definite (SQD) systems.

Square systems
<p align="center">
  <b><i>Ax = b</i></b>
</p>
should be solved when **_A_** has full column rank and **_b_** lies in the range space of **_A_**.
This situation occurs when
   * **_A_** is square and nonsingular, or
   * **_A_** is tall and has full column rank and **_b_** lies in the range of **_A_**.

Linear least-squares problems
<p align="center">
  minimize ‖<b><i>b</i></b> - <b><i>Ax</i></b>‖
</p>
should be solved when **_b_** is not in the range of **_A_**, regardless of the shape and rank of **_A_**.
If there are infinitely many such **_x_** (because **_A_** is rank deficient), identify the one with minimum norm.

Linear least-norm problems
<p align="center">
  minimize ‖<b><i>x</i></b>‖ &nbsp; subject to &nbsp; <b><i>Ax = b</i></b>
</p>
sould be solved when **_A_** is column-rank deficient but **_b_** is in the range of **_A_**.
This situation occurs when **_b_** is in the range of **_A_** and
   * **_A_** is square but singular, or
   * **_A_** is short and wide.

In some applications, adjoint systems
<p align="center">
  <b><i>Ax = b</i></b> &nbsp; and &nbsp; <b><i>Aᵀy = c</i></b>
</p>
saddle-point or symmetric quasi-definite systems
<p align="center">
  [<b><i>M </i></b>&nbsp;&nbsp;&nbsp;<b><i> A</i></b>]&nbsp; [<b><i>x</i></b>]            =           [<b><i>b</i></b>].
  <br>
  [<b><i>Aᵀ</i></b>&nbsp;&nbsp;      <b><i>-N</i></b>]&nbsp; [<b><i>y</i></b>]&nbsp;&nbsp;&nbsp;&nbsp;[<b><i>c</i></b>]&nbsp; 
</p>
are also encountered. Methods specialized for them have been developed too.


Krylov solvers are appropriate, in particular, in situations where such a problem must be solved but a factorization is not possible, either because:
* the operator is not available explicitly,
* the operator is dense, or
* factors would consume an excessive amount of memory and/or disk space.

Iterative methods are particularly appropriate in either of the following situations:
* the problem is sufficiently large that a factorization is not feasible or would be slower,
* an effective preconditioner is known in cases where the problem has unfavorable spectral structure,
* the operator can be represented efficiently as a sparse matrix,
* the operator is *fast*, i.e., can be applied with far better complexity than if it were materialized as a matrix. Often, fast operators would materialize as *dense* matrices.

## How to Install

Krylov can be installed and tested through the Julia package manager:

```julia
julia> ]
pkg> add Krylov
pkg> test Krylov
```
