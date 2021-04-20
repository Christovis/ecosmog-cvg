# ECOSMOG - cubic vector Galileon

[![Build Status](https://travis-ci.com/Christovis/ecosmog-cvg.svg?branch=master)](https://travis-ci.com/Christovis/ecosmog-cvg)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![arXiv](https://img.shields.io/badge/arXiv-2007.03042%20-green.svg)](https://arxiv.org/abs/2007.03042)

A open-source N-body simulation code for dark-matter only cosmological structure formation for cubic vector Galileon model of the [generalised Proca theory](https://arxiv.org/abs/1402.7026) (GP) gravity theory. It is implemented in a modified version of the [ECOSMOG](https://arxiv.org/abs/1110.1379) which is based on RAMSES. It uses adaptive mesh refinement and adaptive time integration to simulate self-gravitating fluids and is massively parallelizable as it makes use of the MPI communication library.

The code has multiple gravitational solvers to choose from, for which one needs to set corresponding values for the four variables extradof, extradof2, extradof3, and extradof4 in namelist/cosmo.nml (default values for these are all .false.) using the following rule:

| Model            | extradof  | extradof2  | extradof3  | extradof4  |
| ---------------- | :-------: | :--------: | :--------: | ---------: |
| LCDM             | False     | False      | False      | False      |
| QCDM             | False     | True       | False      | False      |
| Linear     w/o B | False     | True       | True       | False      |
| Linear     w B   | False     | True       | True       | True       |
| Non-linear w/o B | True      | True       | False      | False      |
| Non-linear w B   | True      | True       | False      | True       |

### Build/Install
The code has ben tested and run on [COSMA](https://www.dur.ac.uk/icc/cosma/), which uses CentOS 7.6 with a 3.10 Linux kernel.

**Prerequisites:**
* A Fortran compiler
* CMake
* Intel MPI (2018)
* IntelComp (2018)

To build the Fortran BMI bindings from source with cmake, run

```
$ cd ./src/bin
$ make
```

You can quickly test your installation by executing:
```
$ cd bin
$ make
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

Initial condition can be generated with [2LPTic](https://arxiv.org/abs/astro-ph/0606505) from a power spectrum generated with [CAMB](https://github.com/cmbant/CAMB).

The simulation results can be analysed through [Astrild](https://github.com/Christovis/astrild).
