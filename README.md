# MSRMD-CGMM
This comprises the essential input files and scripts necessary for conducting constrained density functional theory (CDFT)-based fitting, as well as Multiscale Reactive Molecular Dynamics/Coarse-Grained Molecular Mechanics (MS-RMD/CG-MM) simulations.

Key steps:

1) Constrained Density Functional Theory (CDFT) calculations: These were computed using [CP2K](https://www.cp2k.org/)

2) Diabatic and force matching: Codes provided

3) Coarse-graining: Using [OpenMSCG](https://software.rcc.uchicago.edu/mscg/) for multiscale coarse-graining force-matching

4) MS-RMD/AA-MM and MS-RMD/CG-MM simulations: LAMMPS patched with RAPTOR and PLUMED
- [LAMMPS](https://www.lammps.org/#gsc.tab=0): main classical MD engine
- [RAPTOR](https://software.rcc.uchicago.edu/raptor/home.php): required for reactive molecular dynamics (i.e, where bonding topology can change)
- [PLUMED](https://www.plumed.org/): for enhanced sampling
