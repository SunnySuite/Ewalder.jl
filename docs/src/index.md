```@meta
CurrentModule = Ewalder
```

# Ewalder

Documentation for the Julia package
[Ewalder.jl](https://github.com/SunnySuite/Ewalder.jl) to perform Ewald
summation.

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

### Overview

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
charge density. The _Ewalder_ implementation uses a straightforward algorithm
that scales as $N^{3/2}$ for increasing system sizes $N$ (assuming efficient
calculation of the neighbor list). For very large scale simulations, it would be
preferable to use instead a method that scales near-linearly with $N$, such as
[Smooth Particle Mesh Ewald](https://doi.org/10.1021/ct900275y).

### Mathematical details

The full Ewald method is summarized below. See our [mathematical
writeup](https://raw.githubusercontent.com/SunnySuite/Ewalder.jl/main/docs/math/ewald_review.pdf)
(PDF format) for the derivation. 

The Ewald energy for a system of charges decomposes into three parts,

```math
E = E_S + E_L - E_\mathrm{self}.
```

The "short-range" energy is a sum of pairwise interactions real-space,
```math
E_{S} =\frac{1}{4\pi\epsilon_{0}}\frac{1}{2} \sum_{i,j,\mathbf{n}}{\vphantom{\sum}}' \frac{q_{i}q_{j}}{r_{ij\mathbf{n}}}\mathrm{erfc}\left(\frac{r_{ij\mathbf{n}}}{\sqrt{2}\,\sigma}\right).
```
Observe that ions indexed by $i$ and $j$ can interact across periodic system images, $\mathbf{n} \in \mathbb{R}^3$. The prime in $\sum'$ indicates that infinite self-interactions ($i=j$ and $\mathbf{n}=\mathbf{0}$) should be excluded from the sum. The periodic system may be sheared and the distance between image charges,
```math
\mathbf{r}_{ij\mathbf{n}} = \big|\mathbf{r}_j - \mathbf{r}_i + \sum_\alpha n_\alpha \mathbf{a}_\alpha \big|,
```
involves three (super-) lattice vectors $\mathbf{a}_\alpha$. The complementary error function $\mathrm{erfc}(x)$ decays rapidly in its argument $x$, such that the real-space interactions are localized to $r_{ij\mathbf{n}} \lesssim \sigma$. The length scale $\sigma$ is a tuneable parameter that can be adjusted for numerical efficiency.

The "long-range" energy is a sum over Fourier modes
```math
E_{L}=\frac{1}{2V\epsilon_{0}}\sum_{\mathbf{k}\neq\mathbf{0}}\frac{e^{-\sigma^{2}k^{2}/2}}{k^{2}}\left|\hat{\rho}(\mathbf{k})\right|^{2},
```
involving the Fourier transform of charge density,
```math
\hat{\rho}(\mathbf{k})=\sum_{i}q_{i}e^{-i\mathbf{k}\cdot\mathbf{r}_{i}},
```
and the system volume $V = (\mathbf{a}_1 \times \mathbf{a}_2) \cdot \mathbf{a}_3$. The sum runs overs Fourier modes consistent with the assumed periodicity, i.e. $\mathbf{k} = \sum_\alpha m_\alpha \mathbf{b}_\alpha$ for integer $m_\alpha$ and reciprocal (super-) lattice vectors defined to satisfy $\mathbf{a}_\alpha \cdot \mathbf{b}_\beta = 2\pi \delta_{\alpha \beta}$. The assumption of charge neutrality, $\hat{\rho}(\mathbf{0}) = 0$, justifies the removal of the mode $\mathbf{k}=\mathbf{0}$. Due to the Gaussian decay, only modes with $|\mathbf{k}| \lesssim 1/\sigma$ are relevant to the sum.

The energy $E_L$ can equivalently be viewed as an infinite real-space sum over Gaussian charge clouds. This periodic sum includes artificial "self" interactions,
```math
E_{\mathrm{self}}=\frac{1}{4\pi\epsilon_{0}}\frac{1}{\sqrt{2\pi}\,\sigma}\sum_{i}q_{i}^{2},
```
which must be removed from the final energy output $E$.

Both real- and Fourier- space sums are rapidly convergent. _Ewalder_ truncates these sums using the inequalities,
```math
r_{ij\mathbf{n}} \leq c_0 \sqrt{2} \sigma \\
|\mathbf{k}_{\mathbf{m}}| \leq c_0 \sqrt{2} / \sigma
```
where $c_0$ is a dimensionless parameter (defaulting to $c_0 = 6$ for about 13 digits of accuracy).

Apart from the small truncation error, the energy $E$ is mathematically independent of the parameter $\sigma$. A good choice of $\sigma$ should balance between the cost of the real-space sums in $E_S$ and the Fourier-space sums in $E_L$. Consider systems of fixed density, but with varying sizes controlled by the number of ions $N$. This package selects

```math
\sigma = L / c_1 N^{1/6},
```

where $c_1$ is a dimensionless parameter (defaulting to $c_1 = 2$) and $L = V^{1/3}$. This choice achieves, in principle, a total computational cost that scales as $O(N^{3/2})$, provided that the ion neighbor lists are calculated efficiently. For maximum efficiency, the user should provide a neighbor list. If one is not provided, this package will perform inefficient, brute-force calculation of the neighbor list.

In addition to point charges, _Ewalder_ also allows specification of point dipoles, each of which generates an electrostatic potential that decays like $1/r^2$, and can similarly be handled using Ewald summation. The final formulas can be found in our [PDF note](https://raw.githubusercontent.com/SunnySuite/Ewalder.jl/main/docs/math/ewald_review.pdf).

## API

```@index
```

```@autodocs
Modules = [Ewalder]
```
