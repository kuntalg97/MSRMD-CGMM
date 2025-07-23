# MSRMD-CGMM
This comprises the essential input files and scripts necessary for conducting constrained density functional theory (CDFT)-based fitting, as well as Multiscale Reactive Molecular Dynamics/Coarse-Grained Molecular Mechanics (MS-RMD/CG-MM) simulations.

Key steps:

1) Constrained Density Functional Theory (CDFT) calculations: These were computed using CP2K ([download here](https://www.cp2k.org/)).

2) Diabatic and force matching: Codes provided

3) MS-RMD/AA-MM and MS-RMD/CG-MM simulations: LAMMPS patched with RAPTOR and PLUMED
- LAMMPS ([download here](https://www.lammps.org/#gsc.tab=0)): main classical MD engine
- RAPTOR (download here): required for reactive molecular dynamics (ie, where bonding topology can change)
- PLUMED (download here): for enhanced sampling
