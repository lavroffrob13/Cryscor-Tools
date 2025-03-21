#!/usr/bin/env python
#this is an example of how to read in a FCIDUMP file (i.e. from aperiodic embedding) and use it to do DMRG-based multireference calculations with pySCF and block2. it is adapted from https://block2.readthedocs.io/en/latest/user/dmrg-scf.html#dmrg-sc-nevpt2-multi-state
#another strength that pySCF has over Molpro is allowing states besides singlets in EOM-CCSD, including starting from a UCCSD ground state. The inputs after RHF are much simpler for this (https://pyscf.org/user/cc.html#equation-of-motion-coupled-cluster-theory)
#########READ BEFORE USING##########
#pySCF has some key differences compared to Molpro in terms of what's necessary for reading a FCIDUMP file
#as far as I can tell, you only need mol.atom and mol.basis objects and it doesn't really matter what's in them, whereas Molpro actually cares
#likewise, mol.spin, and mol.charge/mol.nelectron seem to be overwritten by whatever's in the FCIDUMP. If your FCIDUMP file is huge, I recommend using sed to edit the header instead of vim
#YOU NEED TO DELETE THE ORBSYM LINE FROM THE FCIDUMP FILE OR ELSE THE CALCULATION WONT WORK
#I've found DMRG-CASSCF to take a lot of cycles to converge, at least compared to Molpro CASSCF. Set "max_cycle_macro" high

#YOU ALSO NEED TO MAKE SMALL EDITS TO THE PYSCF SOURCE CODE
#1. comment Line 393 to Line 416 in "pyscf/pyscf/mcscf/casci.py" or set a lower verbosity, or set an identity matrix for ovlp_ao, and it will not affect any computation, 
#but if you don't you will get an error here related to " orth_coeff = orth.orth_ao(mc.mol, 'meta_lowdin', s=ovlp_ao)"
#2. replace line 376 "aaaa = ao2mo.kernel(mc.mol, mocas)" in the same file with the following two lines:
#eri_or_mol = mc._scf._eri if mc._scf._eri is not None else mc.mol
#aaaa = ao2mo.kernel(eri_or_mol, mocas)

#some tips from Huanchen
#1. If you want to do CASCI-DMRG, you have to make sure the active space only has < 50 orbitals.
#2. If you need to do CASSCF-DMRG, the typical active space in literature is < 40 orbitals (see https://pubs.acs.org/doi/epdf/10.1021/acs.inorgchem.3c03689, for example).
#3. If you need to do DMRG with dynamic correlation, for example, DMRG-NEVPT2, the active space cannot be larger than 30. https://pubs.acs.org/doi/full/10.1021/acs.jpclett.1c04078 has some MRCI info too

import numpy as np
import pyscf
from pyscf import gto, scf, ao2mo, cc, symm, mcscf, mrpt, dmrgscf, lib
import os
from pyscf.tools import fcidump
#from pyscf.gto.basis import parse_molpro

dmrgscf.settings.BLOCKEXE = os.popen("which block2main").read().strip()
#dmrgscf.settings.BLOCKEXE_COMPRESS_NEVPT = os.popen("which block2main").read().strip() #this can make the calculation faster but was throwing some errors for me on at least one cluster, test
dmrgscf.settings.MPIPREFIX = ''

mol = gto.Mole()
mol.atom = '''
h0 0 0 0
h0 10 0 0
h0 20 0 0
h0 30 0 0
h0 40 0 0
h0 50 0 0
h0 60 0 0
h0 70 0 0
h0 80 0 0
h0 90 0 0
h0 100 0 0
h0 110 0 0
h0 120 0 0
h0 130 0 0
h0 140 0 0
h0 150 0 0
h0 160 0 0
h0 170 0 0
h0 180 0 0
h0 190 0 0
h0 200 0 0
h0 210 0 0
h1 220 0 0
'''
mol.unit = 'bohr'
mol.basis = {'h0': gto.basis.parse('''
C    S
     1              1
C    S
     2              1
C    S
     3              1
C    S
     4              1
C    S
     5              1
C    S
     6              1
C    S
     7              1
C    S
     8              1
C    S
     9              1
'''),'h1': gto.basis.parse('''
C    S
     1              1
C    S
     2              1
C    S
     3              1
C    S
     4              1
C    S
     5              1
''')}
mol.charge = 0
mol.spin = 0
mol.nelectron = 42
mol.verbose = 4
#mol.symmetry = 'C1'
name = 'out'
mol.build()

mf = pyscf.tools.fcidump.to_scf('FCIDUMP', molpro_orbsym=True) #, mf = scf.RHF(mol).newton())
mf.mol.verbose = 5
mf.max_cycle = 300
mf.kernel()

# state average casscf for lowest singlet, triplet, and quintet
lib.param.TMPDIR = os.path.abspath(lib.param.TMPDIR)
solvers = [dmrgscf.DMRGCI(mol, maxM=500, tol=1E-6) for _ in range(3)]
weights = [1.0 / len(solvers)] * len(solvers)

solvers[0].spin = 0
solvers[1].spin = 2
solvers[2].spin = 4
for i, mcf in enumerate(solvers):
    mcf.runtimeDir = lib.param.TMPDIR + "/%d" % i
    mcf.scratchDirectory = lib.param.TMPDIR + "/%d" % i
    mcf.threads = 7
    mcf.memory = int(mol.max_memory / 1000) # mem in GB
mc = mcscf.CASSCF(mf, 15, 16)
mcscf.state_average_mix_(mc, solvers, weights)
mc.max_cycle_macro = 150
mc.canonicalization = True
mc.natorb = True
mc.kernel()
mf.mo_coeff = mc.mo_coeff

# need an extra casci before calling mrpt
mc = mcscf.CASCI(mf, 15, 16)
mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=500, tol=1E-10)
mc.fcisolver.runtimeDir = os.path.abspath(lib.param.TMPDIR)
mc.fcisolver.scratchDirectory = os.path.abspath(lib.param.TMPDIR)
mc.fcisolver.threads = 7
mc.fcisolver.memory = int(mol.max_memory / 1000) # mem in GB
mc.fcisolver.conv_tol = 1e-14
mc.fcisolver.nroots = 3
mc.natorb = True
mc.kernel()

# canonicalization for each state
ms = [None] * mc.fcisolver.nroots
cs = [None] * mc.fcisolver.nroots
es = [None] * mc.fcisolver.nroots
for ir in range(mc.fcisolver.nroots):
    ms[ir], cs[ir], es[ir] = mc.canonicalize(mc.mo_coeff, ci=mc.ci[ir], cas_natorb=False)

# mrpt
for ir in range(mc.fcisolver.nroots):
    mc.mo_coeff, mc.ci, mc.mo_energy = ms[ir], cs, es[ir]
    #mr = mrpt.nevpt2.NEVPT(mc).set(canonicalized=True).compress_approx(maxM=200).run(root=ir)
    mr = mrpt.nevpt2.NEVPT(mc).set(canonicalized=True).run(root=ir)
    print('root =', ir, 'E =', mc.e_tot[ir] + mr.e_corr)
