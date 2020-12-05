# [ECOSMOG](https://arxiv.org/abs/1110.1379) - [cubic vector Galileon](https://arxiv.org/abs/2007.03042)

A open-source N-body simulation code for cosmological structure formation for cubic vector Galileon model of the [generalised Proca theory](https://arxiv.org/abs/1402.7026) (GP) gravity theory. It is implemented in a modified version of the _ECOSMOG_ which is based on _RAMSES_. It uses Adaptive Mesh Refinement to simulate self-gravitating fluids and is massively parallelizable as it makes use of the MPI communication library.

You can quickly test your installation by executing:
```
$ cd bin
$ make
$ cd ..
$ bin/ramses1d namelist/tube1d.nml
```

This version supports different variants of quartic Galileon simulation, for which one needs to set corresponding values for the three variables extradof, extradof2 and extradof3 in namelist/cosmo.nml (default values for these are all .false.) using the following rule:

| Model            | extradof  | extradof2  | extradof3  | extradof4  |
| ---------------- | :-------: | :--------: | :--------: | ---------: |
| LCDM             | False     | False      | False      | False      |
| QCDM             | False     | True       | False      | False      |
| Linear     w/o B | False     | True       | True       | False      |
| Linear     w B   | False     | True       | True       | True       |
| Non-linear w/o B | True      | True       | False      | False      |
| Non-linear w B   | True      | True       | False      | True       |
