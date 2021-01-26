# AMBER EDA

[![GPL license](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat)](https://github.com/CisnerosResearch/AMBER-EDA/blob/main/LICENSE)

This repository contains a Fortran 90 program and associated files for
performing an energy decomposition analysis (EDA) of AMBER molecular dynamics
simulations.

## Compiling `Residue_E_Decomp_openmp.f90`

This program has been tested with the `ifort` compiler.
```
ifort Residue_E_Decomp_openmp.f90 -o Residue_E_Decomp_openmp.x -qopenmp
```

Compilation with `gcc` will likely result in an error.

> **Note**: If you are not using a parm7 (i.e., you have a very old prmtop),
> you will need to comment out the sections starting with the
> comments
> `!  type = 'ATOMIC_NUMBER'`
> `!  type = 'SCEE_SCALE_FACTOR'`, and
> `!  type = 'SCNB_SCALE_FACTOR'`
> prior to compilation.

## :warning: NetCDF Warning :warning:

This program will **not** work on NetCDF files, only the ASCII mdcrd.
Thus, they must be converted with `cpptraj` to `mdcrd`.
While converting, **do not** use `autoimage` or stripping, as
these will affect the results and give absurdly high energies.
Stripping the system, while saving on time, will return radically different
values for both Coulomb and van der Waals energies.
Using autoimage will have an effect on the significant figures, so it should
also not be used.
If your trajectory files are originally written to an ASCII mdcrd, you're good
to go.

Starting with AMBER16, the default format for all generated trajectory and
restart files is NetCDF and not ASCII.
Unless you explicitly set `ioutfm=0` in your `mdin` files,
your trajectory files will be written in NetCDF format.

> **Tip**: You can check if your file is actually an mdcrd by doing
> `file thing.mdcrd` or `head file.mdcrd`.
> `file.mdcrd` will print `file.mdcrd: ASCII text`.
> `head file.mdcrd` will print readable rows of numbers.
>
> If these are secretly a NetCDF in disguise, `file.mdcrd` will print
> `file.mdcrd: data` and `head file.mdcrd` will print absolute garbage
> (non-readable characters).

An example input file for `cpptraj` for this conversion:
```
parm system.prmtop
trajin system_1.nc
trajin system_2.nc

trajout system_full.mdcrd crd
```

## Additional Files

### `EDA_new.inp`

This is the input file for the program that contains information on the
number of atoms and frames.
You can name it however you want, as the name is entered either interactively
or using with a file like `ans.txt`.

The input file contains information on the number of atoms, residues, and files.
While the informational comment specifies `protein` residues and atoms, it is
intended for any atom of interest (including ligands, nucleic
acids, and metals).
The program assumes anything important starts with residue 1 and goes until
the end of the number of protein atoms/protein residues.
The current value for max number of types should be fine without adjustment.
You can also "chain read" multiple trajectory files by having each file on a new
line.

```
473 !number of protein residues
10 !number of files
59663 !total number of atoms
7295 !number of protein atoms
17929 !number of total residues
2000 !max number of types
solvated_complex_md1.mdcrd
solvated_complex_md2.mdcrd
solvated_complex_md3.mdcrd
solvated_complex_md4.mdcrd
solvated_complex_md5.mdcrd
solvated_complex_md6.mdcrd
solvated_complex_md7.mdcrd
solvated_complex_md8.mdcrd
solvated_complex_md9.mdcrd
solvated_complex_md10.mdcrd
```

> **Note**: The chain reading method assumes that either all the trajectory
> files have the same number of samples (frames), **or** that the largest
> trajectory file is read in first.
> If you have trajectory files with differing numbers of frames, read in the
> file with the most frames first.
> The calculation relies on an ensemble average across all frames, so time
> sequence isn't important here.

### `ans.txt`

The program has file prompts, which can be answered as input through the
command-line when running.
This is useful when running the program non-interactively.

The first line of `ans.txt` should contain the filename of the input file
(an example is `EDA_new.inp`, described above).
The second line is the filename of the `prmtop`.

> **Warning**: The program will not work if there are comments in this file.

### `EDA_script.sh`

This is an example PBS script for running the EDA program through
a queuing system.

## Program Output

The program will have 3 output files: `fort_sanity_check.txt`,
`fort_coulomb_interaction.dat`, and `fort_vdw_interaction.dat`.
`fort_coulomb_interaction.dat` contains the Coulomb energies and
`fort_vdw_interaction.dat` contains the van der Waals energies.

## Data Analysis

A number of R-based analysis scripts and their explanations can be found on
[emleddin's GitHub](https://github.com/emleddin/research-scripts/tree/main/amber-analysis/EDA).

## Related Publications

These publications include an explanation of the program/equations and present
results obtained through it.

- [Graham S.E., Syeda FS. Cisneros G.A., “Computational prediction of residues involved in fidelity checking for DNA synthesis in DNA polymerase I”, *Biochemistry*, 51, 2569-2578, 2012.](https://doi.org/10.1021/bi201856m)
- [Elias A., Cisneros G.A., “Computational Study of Putative Residues Involved in, DNA Synthesis Fidelity Checking in Thermus aquaticus DNA Polymerase I” , in *Advances in Protein and Structural Biology: Biomolecular Modelling and Simulations*, Vol. 96, Ed. T. Karabencheva-Christova, pp. 35–79, Elsevier, 2014.](https://doi.org/10.1016/bs.apcsb.2014.06.003)
- [Dewage S.W., Cisneros G.A., “Computational Analysis of Ammonia Transfer along Two Intra-molecular Tunnels in Staphylococcus aureus Glutamine-dependent Amidotransferase (GatCAB)”, *J. Phys. Chem. B*, 119, 3669-3677, 2015.](https://doi.org/10.1021/jp5123568)
- [Leddin E., Cisneros G.A., “Comparison of DNA and RNA Substrate Effects on TET2 Structure”, in “Advances in Protein Chemistry and Structural Biology”, T. Karabencheva-Christova and C. Christov Eds., Elsevier, DOI: 10.1016/bs.apcsb.2019.05.002, 2019.](https://doi.org/10.1016/bs.apcsb.2019.05.002)
