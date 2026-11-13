# MSRMD-CGMM
This comprises the essential input files and scripts necessary for conducting constrained density functional theory (CDFT)-based fitting, as well as Multiscale Reactive Molecular Dynamics/Coarse-Grained Molecular Mechanics (MS-RMD/CG-MM) simulations.

Corresponding paper: [J. Chem. Theory Comput. 2025, 21, 17, 8456–8467](https://doi.org/10.1021/acs.jctc.5c00930)

## Overview

The MS-RMD/CG-MM approach describes a chemical reaction with a reactive model (MS-RMD) while treating the surrounding solvent at a coarse-grained (CG) level. This keeps the chemically important region accurate while greatly reducing the cost of simulating the environment. Results are compared against the corresponding MS-RMD/AA-MM simulations, in which the solvent is represented with all-atom (AA) detail.

The repository covers two S<sub>N</sub>2 reactions in acetone:

| Directory | Reaction |
|---|---|
| `met_chloride/` | Cl⁻ + CH₃Cl → CH₃Cl + Cl⁻ |
| `met_bromide/` | Br⁻ + CH₃Cl → CH₃Br + Cl⁻ |

## Workflow

Key steps:

1) Constrained Density Functional Theory (CDFT) calculations: These were computed using [CP2K](https://www.cp2k.org/)

2) Diabatic and force matching: Codes provided

3) Coarse-graining: Using [OpenMSCG](https://software.rcc.uchicago.edu/mscg/) for multiscale coarse-graining force-matching

4) MS-RMD/AA-MM and MS-RMD/CG-MM simulations: LAMMPS patched with RAPTOR and PLUMED
- [LAMMPS](https://www.lammps.org/#gsc.tab=0): main classical MD engine
- [RAPTOR](https://software.rcc.uchicago.edu/raptor/home.php): required for reactive molecular dynamics (i.e, where bonding topology can change)
- [PLUMED](https://www.plumed.org/): for enhanced sampling

### 1. CDFT reference calculations

CDFT is used to construct two diabatic states for each reaction, one for the reactant (nucleophile free, leaving group bonded to carbon) and one for the product. Each state is defined by a Becke charge constraint on the halogen fragment. A mixed-CDFT calculation then provides, along the reaction coordinate:

- the adiabatic ground-state energy (`Eg`)
- the diabatic energies (`E1`, `E2`)
- the electronic coupling between the states (`V12`)
- the overlap between the states (`S12`)
- the expansion coefficients of the ground state in the diabatic basis (`c1`, `c2`, and their normalized forms)

The reaction coordinate is the difference between the two carbon–halogen distances, `d(C–X1) − d(C–X2)`.

### 2. Diabatic matching and force matching

The MS-RMD model is parameterized against the CDFT data in two stages:

1. **Diabatic matching** (`CDFT/optimization/`): The parameters of the diabatic states (Morse bond terms and a diagonal energy shift) and of the off-diagonal coupling are optimized so that the MS-RMD ground-state eigenvector reproduces the CDFT diabatic coefficients. The overlap matrix from CDFT is included explicitly when the eigenvectors are computed, by Löwdin-orthogonalizing the MS-RMD Hamiltonian.
2. **Force-matching refinement** (`CDFT/optimization/fm_refine_diag/`): With the diabatic parameters fixed, the off-diagonal coupling is refined so that the MS-RMD ground-state forces along the reaction coordinate reproduce the CDFT forces.

At each optimization step the script writes a new LAMMPS/RAPTOR input from a template, evaluates the MS-RMD energies on the CDFT reference geometries, and updates the parameters with the Nelder–Mead algorithm (SciPy).

### 3. Coarse-graining the solvent

The CG acetone model is obtained with multiscale coarse-graining (MS-CG) force matching in OpenMSCG. Each acetone molecule is mapped to CG sites of type `CH3` and `CO`. The resulting tabulated pair interactions are provided as `Pair_*.table` files and are read directly by LAMMPS.

### 4. Reactive simulations and free energy profiles

Umbrella sampling is performed along the same reaction coordinate used in the CDFT calculations, with a harmonic restraint applied through PLUMED. The reactive force field is handled by RAPTOR (`fix evb`), and the simulations are run in the NVT ensemble at 300 K. Both all-atom solvent (`RMD_AA`) and coarse-grained solvent (`RMD_CG`) setups are provided, along with the resulting potentials of mean force.

## Repository structure

Both `met_chloride/` and `met_bromide/` follow the same layout:

```
<system>/
├── CDFT/
│   ├── cdft_files/
│   │   ├── cdft_skeleton/        # CP2K input templates
│   │   ├── total_file_ts.xyz     # reference geometries along the reaction coordinate
│   │   ├── ci_all.dat            # CDFT energies, couplings, overlaps and coefficients
│   │   └── ci_best.dat
│   └── optimization/
│       ├── fit_forces_diag.py    # diabatic matching
│       ├── extract_ci.py         # extract MS-RMD coefficients from RAPTOR output
│       ├── extract_energies.py   # extract MS-RMD energies from RAPTOR output
│       ├── check_fit.py          # compare MS-RMD and CDFT coefficients
│       ├── run_fit.sb            # runs the full fitting sequence
│       ├── sn2_template.inp      # LAMMPS template used during fitting
│       ├── evb.*                 # RAPTOR configuration, topology and parameters
│       └── fm_refine_diag/       # force-matching refinement of the coupling
└── reactive_runs/
    ├── RMD_AA/                   # MS-RMD with all-atom acetone
    │   ├── skeleton/             # input files for one umbrella-sampling window
    │   └── pmf_aa_*.dat          # potential of mean force
    └── RMD_CG/                   # MS-RMD with coarse-grained acetone
        ├── skeleton/             # input files for one umbrella-sampling window
        │                         # (including CG tabulated potentials)
        └── pmf_cg_*.dat          # potential of mean force
```

### Key files

**CP2K (CDFT)**

| File | Description |
|---|---|
| `control.inp~` | Main CP2K input template for single and mixed CDFT calculations |
| `becke_state1.inc`, `becke_state2.inc` | Becke charge constraints defining the two diabatic states |
| `dft-common-params.inc` | Shared DFT settings |
| `subsys.inc` | System and basis set definitions |

**RAPTOR / LAMMPS**

| File | Description |
|---|---|
| `evb.cfg` | Main RAPTOR configuration (output settings, molecule and reaction definitions) |
| `evb.top`, `evb_cg.top` | Reactive topology for the AA and CG systems |
| `evb.type`, `evb_cg.type` | Mapping between RAPTOR and LAMMPS atom types |
| `evb.par`, `evb_model.par` | MS-RMD parameters (diabatic and off-diagonal terms) |
| `in.coeff`, `in.coeff_cg` | Force field coefficients |
| `in.charges`, `in.charges_cg` | Atomic charges |
| `*.data` | LAMMPS data files |
| `sn2*.inp` | LAMMPS input scripts |
| `us.inp` | PLUMED input for umbrella sampling |

## Usage

### Fitting the MS-RMD model

From `CDFT/optimization/` (and then from `fm_refine_diag/` for the refinement stage), run:

```bash
bash run_fit.sb
```

This performs the following steps:

1. `python fit_forces_diag.py` optimizes the parameters and writes the optimized values to `opt_Nelder-Mead.dat`.
2. `python extract_ci.py evb.out` writes the MS-RMD coefficients to `ci.dat`.
3. `python extract_energies.py evb.out` writes the MS-RMD energies to `energy.dat`.
4. `python check_fit.py ci.dat` plots the MS-RMD coefficients against the CDFT coefficients and prints the mean squared deviation.

The fitting requires Python with NumPy, SciPy and Matplotlib, and a LAMMPS executable (`lmp` or `lmp_mpi`) built with RAPTOR.

### Running umbrella sampling

The `skeleton/` directories contain the input for a single umbrella-sampling window. To set up a window, copy the directory and replace the `XXX` placeholders:

- in `us.inp`, with the restraint center along the reaction coordinate
- in the LAMMPS input script, with the label of the starting configuration for that window

Each window is then run with LAMMPS (built with RAPTOR and PLUMED), for example:

```bash
lmp -in sn2_cg_morse.inp
```

The collective variable is written to `CV_K20.dat`. The data from all windows can then be combined into a potential of mean force with a standard method such as WHAM.

## Citation

If you use these files, please cite:

[J. Chem. Theory Comput. 2025, 21, 17, 8456–8467](https://doi.org/10.1021/acs.jctc.5c00930)

## Additional resources

Additional codes for analyzing RAPTOR output files, and the trajectories themselves can be obtained upon reasonable request.
