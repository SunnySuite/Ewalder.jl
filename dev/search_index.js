var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Ewalder","category":"page"},{"location":"#Ewalder","page":"Home","title":"Ewalder","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Ewalder.","category":"page"},{"location":"#Usage-example","page":"Home","title":"Usage example","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Calculate the Madelung constant for NaCl using its primitive cell.","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Ewalder\n\nlatvecs = [[1,1,0], [1,0,1], [0,1,1]]\npos = [[0,0,0], [1,1,1]]\ncharges = [1, -1]\nsys = Ewalder.System(; latvecs, pos)\nE = Ewalder.energy(sys; charges)\n@assert E ≈ -1.74756459","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Ewalder/test directory for more usage examples.","category":"page"},{"location":"#What-is-Ewald-summation?","page":"Home","title":"What is Ewald summation?","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Coulomb energy between two point charges is q_1 q_2  4piepsilon_0 r, where r is the separation distance, and epsilon_0 is a physical constant (the vacuum permittivity). For a system with multiple charges, the total Coulomb energy is a sum of all such pair interactions. A frequent goal is to estimate bulk material properties using computers simulations at much smaller scales. Here, it is very effective to impose periodic boundary conditions on the finite size simulation volume. Doing so effectively creates infinitely many copies (image charges) of the charges present in the original, finite size system. Ewald summation accounts for the interactions of all these charges through long-range (1r) pair interactions. In the Ewalder example above, the output E represents the dimensionless energy per ion of the NaCl crystal, i.e., table salt.","category":"page"},{"location":"","page":"Home","title":"Home","text":"There are some mathematical subtleties to Ewald summation. First, the system must be net charge neutral, or the macroscopic Coulomb energy will diverge. Second, if the original system has a nonzero net dipole moment (i.e., a polarization), then the infinite sum over periodic images becomes only conditionally convergent, and the result depends on the order of summation. This conditional convergence reflects a true ambiguity–-the physical conditions at the surface of the material sample cannot be neglected. The simplest possibility, employed here, is to impose so-called \"tin foil boundary conditions,\" which eliminate bound charge from the sample surface.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Ewald summation works by decomposing the Coulomb energy into two parts: (1) A real-space part, evaluated as a sum over point charges at short range, and (2) a Fourier space part, which captures the long-range interactions through a sum over mathbf k-vectors, which label frequencies in the decomposition of total charge density. The Ewalder implementation uses a straightforward algorithm that scales as N^32 for increasing system sizes N (assuming efficient calculation of the neighbor list). For very large scale simulations, it would be preferable to use instead a method that scales near-linearly with N, such as PPPM or PME.","category":"page"},{"location":"#Mathematical-details","page":"Home","title":"Mathematical details","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"For a full review of Ewald summation, please see our mathematical writeup  (PDF format).","category":"page"},{"location":"#API","page":"Home","title":"API","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Ewalder]","category":"page"},{"location":"#Ewalder.System","page":"Home","title":"Ewalder.System","text":"System(; latvecs::Vector{NTuple{3, Float}},\n         pos::Vector{Vec3}\n         c0::Float64 = 6.0\n         c1::Float64 = 2.0)\n\nCreate a system that is periodic in each of the three supercell lattice vectors latvecs and with atom positions pos. The optional parameter c0 controls accuracy of Ewald summation. The default is c0 = 6, which yields errors of order 1e-12 (smaller c0 will run faster). The optional parameter c1 controls the balance between real and Fourier space summation (larger c1 implies more real-space summation).\n\n\n\n\n\n","category":"type"},{"location":"#Ewalder.energy-Tuple{System}","page":"Home","title":"Ewalder.energy","text":"energy(sys::System; charges=Float64[], dipoles=Vec3[])\n\nCalculate the Ewald energy in units of 14piepsilon_0 for a periodic system of charges and dipoles. If either charge or dipole parameters are omitted, they will be assumed zero.\n\n\n\n\n\n","category":"method"}]
}
