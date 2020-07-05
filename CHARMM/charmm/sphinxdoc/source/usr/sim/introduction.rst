
How to run CHARMM
=================

While you can run CHARMM interactively, one usually tells the program what to
do by means of a script. Under Unix (at least for non-parallel versions of the
program), this means that in order to execute a (short) CHARMM calculation, one
runs from the command line (Unix Shell prompt)

 charmm_executable < charmm_input_script.inp 

exploiting input redirection available under all Unix shells. Since as we shall
see shortly CHARMM output tends to be verbose, one normally also redirects the
output to a file, thus ending up with


 charmm_executable < charmm_script.inp > charmm_output.out 


Of course, instead of **charmm_executable** use the path to the CHARMM
executable you have installed on your computer and replace
**charmm_input_script.inp** and **charmm_output_file.out** by the
names of the actual script which you want to run and the file to which you want
to save your output.

Data Structures
---------------

* *Residue Topology File (RTF)* This file defines groups by including the
  atoms, the properties of the group, and bond and charge information. CHARMM
  has standard Residue Topology Files for nucleic acids, lipids, proteins and
  carbohydrates. An example of a simple RTF, which describes a single residue
  (TIP3P water) is given below.

.. code-block:: rtf

     * Residue topology file for TIP3 water
     *
     31 1

     MASS     4 HT     1.00800 H ! TIPS3P WATER HYDROGEN
     MASS    75 OT    15.99940 O ! TIPS3P WATER OXYGEN

     RESI TIP3         0.000 ! tip3p, generate using noangle nodihedral
     GROUP
     ATOM OH2  OT     -0.834
     ATOM H1   HT      0.417
     ATOM H2   HT      0.417
     BOND OH2 H1 OH2 H2 H1 H2    ! the last bond is needed for shake
     ANGLE H1 OH2 H2             ! required
     ACCEPTOR OH2
     PATCHING FIRS NONE LAST NONE

     END


As you can see, this file containes a title and immediately following it the
rather crypting string "31 1". This is the version number of the topology file,
which is tied to the CHARMM version it was released with. Next comes two MASS
statements, each of which define an atom type. Atom numbers 4 and 75 are
assigned to TIP3P hydrogen and oxygen, respectively. Next comes the actually
definition of the residue, which should be fairly self-explanatory, and then
the file ends with the END keyword.

* *Parameter File (PARA or PARM)* This file determines the energy associated
  with the structure by defining bond, angle and torsion force constants and
  van der Waals parameters. CHARMM has standard parameter files for nucleic
  acids, lipids, proteins  carbohydrates, and water. An example of a parameter
  file with all ofthe parameters needed to simulate a TIP3 water molecule as
  defined above is given here. Note that the atom naming convention in the
  parameter file matches that in the topology file. Failure to uphold the atom
  naming and numbering conventions will yield incorrect results, which is why
  topology and parameter files are released together and it is generally not a
  good idea to mix yopologies and parameters (however, it is possible to append
  one set of topologies and parameters to another).

.. code-block:: rtf

    * parameter file needed to simulate TIP3 water
    *

    BONDS
    !atom type Kb          b0
    HT   HT      0.000     1.5139 ! ALLOW WAT
                    ! FROM TIPS3P GEOMETRY (FOR SHAKE/W PARA
    OT   HT    450.000     0.9572 ! ALLOW   WAT
                    ! FROM TIPS3P GEO

    ANGLES
    !atom types     Ktheta    Theta0   Kub     S0
    HT   OT   HT     55.000   104.5200 ! ALLOW WAT
                    ! TIP3P GEOMETRY, ADM JR.

    DIHEDRALS

    IMPROPER

    NONBONDED nbxmod  5 atom cdiel shift vatom vdistance vswitch -
      cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5

    !atom  ignored    epsilon    Rmin/2 ignored   eps,1-4  Rmin/2,1
    OT     0.000000  -0.152100     1.768200 ! ALLOW   WAT
                  !TIP3P OXYGEN PARAMETERS, adm jr., NBFIX obsolete
    HT     0.000000  -0.046000     0.224500 ! ALLOW WAT
                  !TIP3P HYDROGEN PARAMETERS, adm jr., NBFIX obsolete

    END

Note that there are no dihedral or improper dihedrals parameters necessary for
TIP3 water as there are only 3 atoms in the residue. Some parameter files also
contain CMAP parameters, which are 2-dimensional grid corrections for dihedral
angles (see [MacKerell04]_ for further details).

* *Coordinates (COOR)* These are the standard Cartesian coordinates of the
  atoms in the system. These are typically read in or written out in PDB or
  CHARMM card (CRD -- the default file format used throughout CHARMM) file
  format. The card format keeps track of additional molecule information that
  can be useful for structure manipulation (*i.e.* residue name, segment name,
  segment id, resdiue id, etc.). Below is an example of a .crd file and the
  information in contains::

    title = * WATER
    title = *  DATE:     4/10/07      4:25:51      CREATED BY USER: USER
    title = *
    Number of atoms (NATOM)       = 6
    Atom number (ATOMNO)          = 1 (just an exmaple)
    Residue number (RESNO)        = 1
    Residue name (RESName)        = TIP3
    Atom type (TYPE)              = OH2
    Coordinate (X)                = -1.30910
    Coordinate (Y)                = -0.25601
    Coordinate (Z)                = -0.24045
    Segment ID (SEGID)            = W
    Residue ID (RESID)            = 1
    Atom weight (Weighting)       = 0.00000

now what the CHARMM crd file containing that information looks like...

.. code-block:: crd

    * WATER
    *  DATE:     4/10/07      4:25:51      CREATED BY USER: USER
    *
        6
     1    1 TIP3 OH2   -1.30910  -0.25601  -0.24045 W    1      0.00000
     2    1 TIP3 H1    -1.85344   0.07163   0.52275 W    1      0.00000
     3    1 TIP3 H2    -1.70410   0.16529  -1.04499 W    1      0.00000
     4    2 TIP3 OH2    1.37293   0.05498   0.10603 W    2      0.00000
     5    2 TIP3 H1     1.65858  -0.85643   0.10318 W    2      0.00000
     6    2 TIP3 H2     0.40780  -0.02508  -0.02820 W    2      0.00000

* *Protein Structure File (PSF)* The PSF holds lists of every bond, bond angle,
  torsion angle, and improper torsion angle as well as information needed to
  generate the hydrogen bonds and the non-bonded list. It is essential for the
  calculation of the energy of the system.

* *Internal Coordinates (IC)* This data structure defines the internal
  coordinates for atoms and can be used for analysis. Internal coordinates
  represent the position of atoms relative to one another rather than relative
  to Cartesian axes. In many cases, it is not necessary to deal directly with
  the internal coordinate data structure, however it is possible to manipulate
  it within a CHARMM script.

* *Non-Bonded list (NBONds)* This is a atoms which are not bound to each other.
  It is used in calculating the non=bonded energy terms and electrostatic
  properties. The non-bonded list does contain atoms that are in atom-to-atom
  contact and engaging in van der Waals interactions.

* *Constraints (CONS)* Constraints fix atoms in exactly one position during the
  simulation. This information is stored internally in the IMOVe array.

* *Images Data Structures (IMAGe)* This data structure is used to help create
  symmetrical structures and contains bond information. This is a general image
  support system that allows the simulation of almost any crystal and also
  finite point groups. There is also a facility to introduce bond linkages
  between the primary atoms and image atoms. This allows infinite polymers,
  such as DNA to be studied. For infinite systems, an asymmetric unit may be
  studied because rotations and reflections are allowed transformations.

* *Crystal Data Structures (CRYStal)* The crystal module is an extension of the
  image facility within CHARMM that allows calculations on crystals to be
  performed. It is possible to build a crystal with any space group symmetry,
  to optimize its lattice parameters and molecular coordinates and to carry out
  analysis of the vibration spectrum of the entire crystal similar to normal
  mode analysis. All crystal commands are invoked by the keyword CRYStal.

Basic CHARMM script elements
============================

Titles
------

First, let's do something really silly and start up charmm reading from an
empty file; which can be easily accomplished by executing

 charmm_executable < /dev/null

CHARMM prints a header telling you copyright info, version and some more stuff,
followed by a warning

::

             Chemistry at HARvard Macromolecular Mechanics
               (CHARMM) - Developmental Version 35b3     August 15, 2008
   Copyright(c) 1984-2001  President and Fellows of Harvard College
                          All Rights Reserved
  Current operating system: Linux-2.6.18-128.7.1.el5(x86_64)@n138.lobos.[+ 15]
             Created on 12/ 2/ 9 at  2:23:40 by user: tim

        Maximum number of ATOMS:    360720, and RESidues:      120240
        Current HEAP size:        10240000, and STACK size:  10000000

  RDTITL> No title read.

      ***** LEVEL  1 WARNING FROM <RDTITL> *****
      ***** Title expected.
      ******************************************
      BOMLEV (  0) IS NOT REACHED. WRNLEV IS  5

The job finishes by printing some status info. The interesting part is the
warning from which we learn that CHARMM expected a "title". Indeed, each CHARMM
script should start with a title, and if the main script tells CHARMM to read
from another file, the program also expects to find a title at the beginning of
that file. 

A title should not be confused with comments. E.g., it can only occur at the
beginning of a file (we'll explain the apparent exceptions when we encounter
them). Title lines start with a star or asterisk (*); to indicate the end of
the title give a line containing only a star. (A title can consist of up to 32
consecutive lines) 
Thus,

.. code-block:: chm

 * This would be a short title
 * 

If you start CHARMM with a short file containing the above snippet (=title),
you get the title echoed in uppercase letters ::

 RDTITL> * THIS WOULD BE A SHORT TITLE
 RDTITL> *


instead of the warning when using the empty file.

Comments
--------

Having blabbered so much about titles, what are comments: A comment in a CHARMM
script is everything following an exclamation mark *i.e.*, 

.. code-block:: chm

 ! this is a comment on a line by itself

and this would be a line containing a CHARMM command, followed by a comment

.. code-block:: chm

 ENERgy ! as you might expect, this command calculates an energy

Ending a CHARMM script
----------------------

So far, CHARMM finished when it reached the end of the script file (as the line ::

                    NORMAL TERMINATION BY END OF FILE

in the output informs you. It's OK to end a CHARMM script in this manner, but the preferred way of stopping CHARMM is by issuing the 

.. code-block:: chm

 stop

command. We, thus, can create a first, completely useless CHARMM script, looking like 

.. code-block:: chm

 * A CHARMM script doing nothing
 *

 ! we really should be doing something, but in the meantime 
 ! all we know is

 stop

In addition to the title, the comment is echoed as well. Note that CHARMM now
prints upon finishing ::

  NORMAL TERMINATION BY NORMAL STOP

This indicates that the CHARMM script has finished successfully; an
"abnormal" termination message will print if CHARMM exits with an error.

