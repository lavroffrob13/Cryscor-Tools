# Cryscor-Tools
Tools for ease of use in aperiodic embedding, etc in Cryscor

1. make_inputs2.py: run "python make_inputs2.py INPUT" where INPUT is a Crystal input file. Will generate and submit jobs on the Rigi cluster for ghost-basis and full-basis HF calcs, Wannierization, and a dummy embedding calc to help you pick atoms of the embedded fragment

2. posthf_inp_den.py: run "python posthf_inp_den.py FCIDUMP" to generate a molpro ccsd(t) input for the FCIDUMP file generated in Cryscor. pySCF version to come. Ccsd(t) can of course be changed to any ground or excited state energy calc using canonical orbitals (local MP2, RPA, CIS, CCSD, DCSD and ACD(2) are available using a different interface). Since the FCIDUMP file does not contain any derivative integrals, forces and any other properties besides energy cannot be calculated. Analytic gradients in Cryscor will be available in the future.

3. vasp2crystal_gen.py: Converts a vasp slab POSCAR (ie with a vacuum) to a Crystal slab input/d12 file (Gaussian basis, not plane wave, so no vacuum needed). Run as "python vasp2crystal_gen.py #" where # is the lattice vector of the long (ie vacuum) axis

4. displace_along_bond.py: Moving atoms in a symmetric way (i.e. if a neighbor to an anionic vacancy is repulsed from the anion but all other bond lengths should remain symmetric) gets annoying to do by hand, especially if the bond doesn't run parallel to any lattice vectors. Use this script instead.

5. More to come, for OSV-LPM2 + periodic embedding release, etc

6. An informal tutorial was given on 3/20/2025 and the recording is available here https://ucla.box.com/s/h6rv2ioudvoo12mm5nk8xa2nz6pnxcjv

7. The general workflow is outlined in the cheat sheet in this repo
