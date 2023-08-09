# SMT Experiment

**SMT** (Scattering Matrix Tomography) uses scattering matrices to achieve high-resolution deep imaging inside scattering media with confocal spatio-temporal gates, digital pulse compression, refractive-index mismath correction, and additional wavefront corrections for bot input and output.

This repository includes codes used in the epxerimental implementation of SMT.
Using [<code>run_all.m</code>](./reconstruct_SMT/run_all.m), one can generate SMT images from the experimentally measured hyperspectral reflection matrix.
Detailed documentation is given in the comments of the function files.

For running the code, one will need input data of reflection matrices which can be downloaded via
* [doi.org/10.25452/figshare.plus.22109147](https://plus.figshare.com/articles/dataset/Dataset_for_deep_volumetric_imaging_by_scattering_matrix_tomography_-_hyperspectral_reflection_matrices_of_a_USAF_target_under_1-millimeter_mouse_brain_tissue/22109147/1)
(Reflection matrices of a 2D target under 1-millimeter biological tissue)
* [doi.org/10.25452/figshare.plus.22114424](https://plus.figshare.com/articles/dataset/Dataset_for_deep_volumetric_imaging_by_scattering_matrix_tomography_-_hyperspectral_reflection_matrices_of_a_dense_colloid/22114424/1)
(Reflection matrices of a dense colloid)

To run this code, one needs to install Flatiron NUFFT and NLopt. The detailed installation is given here
* Flatiron NUFFT: https://finufft.readthedocs.io/en/latest/install.html
* NLopt: https://nlopt.readthedocs.io/en/latest/NLopt_Installation/
