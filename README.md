# Ewalder

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://SunnySuite.github.io/Ewalder.jl/stable/) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://SunnySuite.github.io/Ewalder.jl/dev/)
[![Build Status](https://github.com/SunnySuite/Ewalder.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/SunnySuite/Ewalder.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/SunnySuite/Ewalder.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/SunnySuite/Ewalder.jl)

A reference implementation for Ewald summation, which calculates the electrostatic energy of a periodic system.

- Support for charges and dipoles.
- Support for sheared volumes.
- Tin foil boundary conditions.
- Errors can be controlled to order $10^{-12}$.

## Example

Calculate the [Madelung
constant](https://en.wikipedia.org/wiki/Madelung_constant) for NaCl using its
primitive cell.

```julia
import Ewalder

latvecs = [[1,1,0], [1,0,1], [0,1,1]]
pos = [[0,0,0], [1,1,1]]
charges = [1, -1]
sys = Ewalder.System(; latvecs, pos)
E = Ewalder.energy(sys; charges)
@assert E ≈ -1.74756459
```

For more explanation, see the [package documentation](https://SunnySuite.github.io/Ewalder.jl/dev/).

Usage examples are contained in [`Ewalder/test`](https://github.com/SunnySuite/Ewalder.jl/tree/main/test).
