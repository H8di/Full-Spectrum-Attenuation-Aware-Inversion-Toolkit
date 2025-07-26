# Seismic Inversion & Rock Physics Tools (wide range of frequency)

## 🎯 Research Focus

This repository supports my research in **Seismic Inversion**, **Elastic Tensor Analysis**, **Attenuation** and **Geomechanical Modeling** with a focus on subsurface characterization at wide range of scales. The codebase contributes to numerical simulations and inversion workflows involving:

- Effective Medium Theory and Rock Physics
- Elastic wave attenuation in isotropic & Anisotropic Media
- Eshelby-based inclusions and T-matrix formulation
- Green’s function computations and PDE-based modeling
- Finite Difference approximations for forward simulations

These tools allow building synthetic datasets, transforming stiffness tensors, analyzing wave response, attenuation, and validating inversion techniques. They provide a versatile foundation for custom seismic inversion in both isotropic and anisotropic geological settings.


This repository contains a subset of custom MATLAB tools developed for seismic inversion and rock physics modeling at differenr frequency range from core, well-logs, and seismic scale. The full set of modules involves elastic tensor manipulation, Green's functions, Eshelby inclusions, transformation matrices, T-Matrix Model, and inversion algorithms.

## ✅ Publicly Available Modules

The following files are included in this public repository to demonstrate code quality, structure, and MATLAB proficiency:

| File | Description |
|------|-------------|
| `Hash_Sht.m` | MATLAB function forComputes effective elastic moduli using the Hashin-Shtrikman bounds for composite materials. |
| `KELVIN_c4dto2d.m` | MATLAB function for Converts a 4th-order stiffness tensor to 2D matrix form in the Kelvin notation. |
| `paper1_model.m` | Main modeling script used for simulations and analysis presented in the first research paper (Paper 1). |
| `s4dto2d.m` | Converts 4D stress/strain tensor to 2D matrix. |
| `tau_calc.m` | MATLAB function for Compute stress relaxation time  based on rock physics parameters. See full documentation for usage. |
| `Voigt_Reuss.m` | Calculates effective properties using Voigt and Reuss bounds. |

## 📂 Module Descriptions

| File | Description |
|------|-------------|
| `Chaotic.m` | MATLAB function for Computes statistical distribution of fracture orientations and magnitudes under chaotic regimes. |
| `Christoffel.m` | Solves Christoffel equation for wave propagation in anisotropic media. |
| `Eshelby.m` | Calculates Eshelby tensor for inclusions in elastic media. |
| `EulurAngle.m` | Generates rotation matrix using Euler angles. |
| `GSA.m` | MATLAB function for Implements the Generalized Self-Consistent Approximation to estimate the effective stiffness of composite. |
| `Green_4D.m` | Computes 4D Green’s function for elastic response modeling. |
| `Hash_Sht.m` | MATLAB function forComputes effective elastic moduli using the Hashin-Shtrikman bounds for composite materials. |
| `KELVIN_c4dto2d.m` | MATLAB function for Converts a 4th-order stiffness tensor to 2D matrix form in the Kelvin notation. |
| `K_from_C.m` | Extracts bulk modulus from stiffness tensor C. |
| `Maxwell_Boltzmanpdf.m` | Computes Maxwell-Boltzmann distribution function. |
| `Mu_from_C.m` | Extracts shear modulus from stiffness tensor C. |
| `Rotate_Euler2d.m` | MATLAB function for modeling or applies 2D rotation to tensors or vectors using Euler angles for coordinate transformation. |
| `Voigt_Reuss.m` | Calculates effective properties using Voigt and Reuss bounds. |
| `Woods.m` | MATLAB function for compute effective bulk modulus of fluid mixtures using Wood’s formula for compressibility averaging. |
| `ab_gampdf.m` | Evaluates the gamma probability density function (PDF) for given shape (a) and scale (b) parameters. |
| `c2dto4d.m` | Converts 2D compliance matrix to 4D tensor. |
| `c4dto2d.m` | Converts 4D stiffness tensor to 2D matrix. |
| `contract.m` | Performs tensor contraction operation. |
| `elastictotensor.m` | Converts elastic constants into tensor representation. |
| `equality_detect.m` | MATLAB function for modeling or inversion tasks. See full documentation for usage. |
| `green2d.m` | 2D Green’s function for layered media. |
| `green4d.m` | Alternate 4D Green’s function implementation. |
| `modified_communicative_t_check.m` | MATLAB function for evaluate the commutative T-matrix model based on the formulation by Jakobsen (2003). |
| `paper1_model.m` Main modeling script used for simulations and analysis presented in the first research paper (Paper 1). |
| `paper1_test_a1_a2_a3_sep_w.m` | Testing different input cases for Paper 1 model. |
| `rotate2d.m` | Rotates 2D tensors or vectors using a specified angle for coordinate frame transformation; compatible with all statistical models including Gaussian, Gamma, Maxwell, and chaotic orientation. |
| `rotate_check.m` | Validation script for tensor rotation consistency. |
| `s2dto4d.m` | Converts 2D stress/strain matrix to 4D tensor. |
| `s4dto2d.m` |  Converts 4D stress/strain tensor to 2D matrix. |
| `t_matrix.m` | Computes the T-matrix for elastic wave scattering and effective property estimation in media with inclusions. |
| `t_matrix_perm.m` | T-matrix variant for permittivity or permeability effects. |
| `tau_calc.m` | MATLAB function for Compute stress relaxation time  based on rock physics parameters. See full documentation for usage. |

## 📊 **Figures Folder:**  
The [`graphs/`](./graphs) directory contains visual outputs illustrating the impact of various effective medium parameters on seismic wave velocity and attenuation across different frequency ranges. These plots reflect frequency-dependent behaviors modeled using the provided inversion and rock physics tools.

## 🔒 Full Codebase Access (Private Repository)

To maintain research integrity and prevent misuse, only selected files are shared publicly.

The full MATLAB codebase (30+ custom modules) includes:
- 4D/2D tensor conversions
- Green's functions in isotropic/anisotropic media
- Eshelby tensors and inclusion mechanics
- Christoffel matrix solvers
- Seismic response simulators
- Inversion routines and optimization modules

### 📧 Request Full Access
If you're an admissions committee member, faculty, or researcher wishing to evaluate the full codebase for academic purposes, please email me:

**Email:** Alimohamadi2@gmail.com  
**GitHub Username:** Include your GitHub username in the request.

Access will be granted via GitHub's private repository system.

## ⚠️ Licensing & Terms

All code in this repository is provided for evaluation only.  
Reproduction, redistribution, or use for research without permission is strictly prohibited.
