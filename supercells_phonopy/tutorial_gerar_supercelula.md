# How to run
The procedure of QE-phonopy calculation is shown below using the NaCl example found in example/NaCl-QE directory.

Read a QE-PW input file and create supercells with --qe option:

> % phonopy --qe -d --dim="2 2 2" -c NaCl.in
---
&control
    calculation = 'scf'
    tprnfor = .true.
    tstress = .true.
    pseudo_dir = '/home/togo/espresso/pseudo/'
 /
 &system
    ibrav = 0
    nat = 8
    ntyp = 2
    ecutwfc = 70.0
 /
 &electrons
    diagonalization = 'david'
    conv_thr = 1.0d-9
 /
ATOMIC_SPECIES
 Na  22.98976928 Na.pbe-spn-kjpaw_psl.0.2.UPF
 Cl  35.453      Cl.pbe-n-kjpaw_psl.0.1.UPF
ATOMIC_POSITIONS crystal
 Na   0.0000000000000000  0.0000000000000000  0.0000000000000000
 Na   0.0000000000000000  0.5000000000000000  0.5000000000000000
 Na   0.5000000000000000  0.0000000000000000  0.5000000000000000
 Na   0.5000000000000000  0.5000000000000000  0.0000000000000000
 Cl   0.5000000000000000  0.5000000000000000  0.5000000000000000
 Cl   0.5000000000000000  0.0000000000000000  0.0000000000000000
 Cl   0.0000000000000000  0.5000000000000000  0.0000000000000000
 Cl   0.0000000000000000  0.0000000000000000  0.5000000000000000
CELL_PARAMETERS angstrom
 5.6903014761756712 0 0
 0 5.6903014761756712 0
 0 0 5.6903014761756712
K_POINTS automatic
 8 8 8 1 1 1
---

In this example, 2x2x2 supercells are created. `supercell.in` and  `supercell-xxx.in` (xxx are numbers) give the perfect supercell and supercells with displacements, respectively.
In the case of the NaCl example, two files supercell-001.in and supercell-002.in are created. In these supercell files, lines only relevant to crystal structures are given.
phonopy_disp.yaml is also generated, which contains information about supercell and displacements.

To make QE-PW input files, necessary setting information is added to supercell-xxx.in files, e.g., by:

> % for i in {001,002};do cat header.in supercell-$i.in >| NaCl-$i.in; done
---
 &control
    calculation = 'scf'
    tprnfor = .true.
    tstress = .true.
    pseudo_dir = '/home/togo/espresso/pseudo/'
    disk_io = 'none'
 /
 &system
    ibrav = 0
    nat = 64
    ntyp = 2
    ecutwfc = 70.0
 /
 &electrons
    diagonalization = 'david'
    conv_thr = 1.0d-9
 /
K_POINTS automatic
 2 2 2  1 1 1
---


where header.in is specially made for this NaCl example and this file is found in example/NaCl-QE directory.
This setting is of course dependent on systems and has to be written for each interested system.
Note that supercells with displacements must not be relaxed in the force calculations, because atomic forces induced by a small atomic displacement are what we need for phonon calculation.

Then QE-PW supercell calculations are executed to obtain force on atoms, e.g., as follows:

% mpirun pw.x -i NaCl-001.in |& tee NaCl-001.out
% mpirun pw.x -i NaCl-002.in |& tee NaCl-002.out
