```@meta
CurrentModule = Ewalder
```

# Ewalder

Documentation for [Ewalder](https://github.com/SunnySuite/Ewalder.jl).

## Usage example

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

See the [`Ewalder/test`](https://github.com/SunnySuite/Ewalder.jl/tree/main/test) directory for more usage examples.

## What is Ewald summation?

The Coulomb energy between two point charges is $q_1 q_2 / 4\pi\epsilon_0 r$,
where $r$ is the separation distance, and $\epsilon_0$ is a physical constant
(the vacuum permittivity). For a system with multiple charges, the total Coulomb
energy is a sum of all such pair interactions. A frequent goal is to estimate
bulk material properties using computers simulations at much smaller scales.
Here, it is very effective to impose periodic boundary conditions on the finite
size simulation volume. Doing so effectively creates infinitely many copies
(image charges) of the charges present in the original, finite size system.
Ewald summation accounts for the interactions of all these charges through
long-range ($1/r$) pair interactions. In the _Ewalder_ example above, the output
`E` represents the dimensionless energy per ion of the NaCl crystal, i.e., table
salt.

There are some mathematical subtleties to Ewald summation. First, the system
must be net charge neutral, or the macroscopic Coulomb energy will diverge.
Second, if the original system has a nonzero net dipole moment (i.e., a
polarization), then the infinite sum over periodic images becomes only
conditionally convergent, and the result depends on the order of summation. This
conditional convergence reflects a true ambiguity---the physical conditions at
the surface of the material sample cannot be neglected. The simplest
possibility, employed here, is to impose so-called "tin foil boundary
conditions," which eliminate bound charge from the sample surface.

Ewald summation works by decomposing the Coulomb energy into two parts: (1) A
real-space part, evaluated as a sum over point charges at short range, and (2) a
Fourier space part, which captures the long-range interactions through a sum
over $\mathbf k$-vectors, which label frequencies in the decomposition of total
charge density. The `Ewalder` implementation uses a straightforward algorithm
that scales as $N^{3/2}$ for increasing system sizes $N$ (assuming efficient
calculation of the neighbor list). For very large scale simulations, it would be
preferable to use instead a method that scales near-linearly with $N$, such as
PPPM or PME.

## Mathematical details

For a full review of Ewald summation, please see our [mathematical writeup
](https://raw.githubusercontent.com/SunnySuite/Ewalder.jl/docs/math/ewald_review.pdf)
(PDF format).

## API

```@index
```

```@autodocs
Modules = [Ewalder]
```
