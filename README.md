# ParaSiF_CHT_Development

This repository contains the ongoing development of ParaSiF new feature on **Conjugate Heat Transfer** between FEniCSx and OpenFOAM that to be integrated into the ParaSiF Parallel Partitioned Simulation Framework once mature.
It is maintained as a **submodule** of the main ParaSiF repository: [ParaSiF Main Repository](https://github.com/ParaSiF/ParaSiF). Please note that all dev-submodules will not be shown in the main ParaSiF repository for users.

---

## Overview

Ongoing development of ParaSiF new feature on **Conjugate Heat Transfer** between FEniCSx and OpenFOAM via the **[MUI coupling library](https://mxui.github.io/)**.

This codebase is currently comprised of two solvers which are git submodules to this repository. For the solid region a FEM based solver, [heatSolverFenicX](https://github.com/blairSmcc03/heatSolverFenicsX), is used. For the fluid region Openfoam is used with [Coupled custom boundary conditions](https://github.com/blairSmcc03/customCHTBoundaryConditions) which perform the coupling operations.
---

## Future Improvements

This code is still a work in progress, there several improvements that could be made such as:

- Adding Aitken/fixed relaxation to the MUI fetch to improve convergence during strong coupling.
- Adding RBF sampler for heat flux across interface to support non-conforming meshes
- Re-integrate support for temporal sampling when using Strong Coupling (this is not really well defined)
- Optimise FEniCSx Computational Efficiency.

## Repository Structure

```
ParaSiF/dev/CHT/
├── fluid/                # fluid OpenFOAM solvers folder
│ ├── src/                # ParaSiF-specific OpenFOAM source code folder
│ │ ├── solvers/          # ParaSiF-specific OpenFOAM solvers
│ │ └── libs/             # Libs used by OpenFOAM solvers
│ │ |  ├── customCHTBoundaryConditions/          # custom coupled boundary conditions for openfoam
│ └── test/               # OpenFOAM unit test folder
├── structure/            # structure FEniCSx solvers folder
│ ├── src/                # ParaSiF-specific FEniCSx source code folder
│ │ |  ├── heatSolverFenicsX/                   # Solver for the heat equation in solids using FenicsX 
│ └── test/               # FEniCSx unit test folder
└── heatSolverFEniCSx_OpenFOAM              # examples of FEniCSx-OpenFOAM coupling
```

---
## Dependencies

The project combines the dependencies of both submodules. 

**customCHTBoundaryConditions** depends on Openfoam v2506 and MUI. 

heatSolverFenicsX depends on the following packaages. It was tested using Python 3.12.

Required packages:
- `scipy v1.16.3`
- `dolfinx v0.9.0`
- `ufl v2024.2.0`
- `basix v0.9.0`
- `mpi4py v4.1.1`
- `petsc4py v3.24.3`
- `numpy v2.3.5`
- `mui4py @ master`
  
---

## Installation

1. First download the code via git:
   ```bash
   git clone https://github.com/blairSmcc03/ParaSiF_CHT_Development.git
   git submodule init
   git submodule update
   ```
2. Install the dependencies. The best way to do this is using [spack](https://spack.io/)
   ```bash
   spack env create <environment name> spack.yaml
   spack env activate -p <environment name>
   spack install   # this will take a long time (compiling openfoam, fenics and mui)
   ```
## Usage and Example Cases

Once you have an up-to-date spack environment it is easy to run the example cases to test your implementation. To run caseA for example:
```bash
cd heatSolverFEniCSx_OpenFOAM/caseA
./Allrun
```
You can then view the output in Paraview by opening "fluid/fluid.foam" and "solid/output/fenicsx_solid_data.xdmf".

To create a new case follow the existing format with the openfoam case in the "fluid" directory and a heatSolverFenicsX case (format specified in heatSolverFenicsX/README.md) in the "solid" directory.

For now the mesh structure is limited to a simple box.

## Contributing

ParaSiF, including this repository, is an **open-source project**, and contributions from the community are warmly welcomed.

There are many ways you can help improve this submodule, including:

- Adding new features, libs or solvers
- Improving documentation, tests and examples
- Fixing bugs or refining existing functionality
- Sharing feedback and suggestions for enhancements

Your contributions, whether large or small, are highly valued and help make ParaSiF a stronger resource for the research community.

For detailed guidance on contributing, please see the [CONTRIBUTING.md](https://github.com/ParaSiF/ParaSiF/blob/main/CONTRIBUTING.md) in the main ParaSiF repository.

## License

Copyright (C) 2021–2025 The ParaSiF Development Team.  
Licensed under the **GNU General Public License v3 (GPL-3.0)**.

## Contact

For questions or contributions, please contact the ParaSiF Development Team

