# SMT Experiment

**SMT** (Scattering Matrix Tomography) uses scattering matrices to achieve high-resolution deep imaging inside scattering media with confocal spatio-temporal gates, digital pulse compression, refractive-index mismatch correction, and additional wavefront corrections for bot input and output.

This repository includes codes used in the experimental implementation of SMT.
Using [<code>run_all.m</code>](./reconstruct_SMT/run_all.m), one can generate SMT images from the experimentally measured hyperspectral reflection matrix.
Detailed documentation is given in the comments of the function files.

For running the code, one will need input data of reflection matrices which are available on [Zenodo](https://doi.org/10.5281/zenodo.8231065).

To run this code, one needs to install Flatiron NUFFT and NLopt. The detailed installation is given here
* Flatiron NUFFT: https://finufft.readthedocs.io/en/latest/install.html
* NLopt: https://nlopt.readthedocs.io/en/latest/NLopt_Installation/
