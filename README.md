# ED Circuits

A julia package for the exact diagonalization (ED) simulation of the quantum circuits featuring

* Many choices for state initialization
* Iterative gate application, including with Trotterization schemes
* Random gate generation
* Mid-circuit measurements
* Multiple observables, such as entanglement and Rényi entropies, mutual informaiton and stabilizer Rényi entropies
* Projected ensemble generation and associated moments and frame potential

The package was created for the simulations in the following works:

**H Lóio, G Lami, L Leone, M McGinley, X Turkeshi, J De Nardis,
"Quantum State Designs via Magic Teleportation",
arXiv:2510.13950
https://doi.org/10.48550/arXiv.2510.13950**

**G Cecile, H Lóio, J De Nardis,
"Measurement-induced phase transitions by matrix product states scaling"
Phys. Rev. Research 6, 033220
DOI: https://doi.org/10.1103/PhysRevResearch.6.033220** 

If you use this code in your research please cite the relevant works.

## User guide

Install with `./install.sh` and uninstall with `./unistall.sh`.

The package is poorly documented but you can find example scripts in the `examples` directory.
