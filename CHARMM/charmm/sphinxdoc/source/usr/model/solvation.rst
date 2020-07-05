.. _tut-solvation:

Solvation
=========

Solvation is surrounding the solute, generally a protein, nucleic acid strand,
or other macromolecule with a solvent, typically water. This is very useful in
biomolecular simulations since biochemical reactions generally take place in a
viscous environment. Though CHARMM supports methods of estimating solvent
effects implicitly, in many simulations it is useful to be able to explicitly
place the system being studied into a solvated environment. There is no
specific command for solvation, but it can be done through a series of
commands. The basic procedure for solvation, as implemented by CHARMMing, is:

* Read in the system to be solvated and center it
* Read in a large box of pre-equilibrated waters (CHARMMing uses TIP3, which is
  the residue name for the default TIP3 water model)
* Delete the waters that overlap the solute (any water whose Oxygen is within
  2.5 Å of the solute) and are more than a certain distance from the center
  (the exact distance depends on the size of the solvent box chosen by the
  user), leaving a spherical water structure
* Unless the shape of the water structure is to be a sphere, delete any
  remaining waters that are outside of the final unit cell

This procedure will be discussed below. A full-fledged example showing the
exact CHARMM commands is given in the [[FINAL Full example|full example]].

How much solvent is needed
--------------------------

As a general rule, biochemists are not particularly fascinated in the behavior
of bulk solvent, however, they are very interested with the way that various
macromolecules behave in its presence. This creates something of a conundrum:
there needs to be enough solvent present to allow the system of interest to
interact with it naturally, but not so much that the system becomes too big to
simulate on the computing power available (the solvent-solvent interactions add
to the computing power needed to simulate the system).  Using periodic boundary
conditions (PBC) does not relieve the user of the need to deal with this
problem: if there is not enough solvent to fill the unit cell then unphysical
"vacuum pockets" will be created; if the unit cell is too small then the
macromolecule may interact with images of itself.

There is no hard and fast rule on how much solvent is needed. The decision
should based on the size and shape of structure and how much it is expected to
move during the simulation. If the system of interest is expected to move
significantly and PBC is not being used, then it is necessary to provide a
buffer of solvent sufficient to accommodate its movement in any direction plus
a buffer zone of at least 12 Å. It is important to note that there can be
substantial artifacts in the first couple of layers of solvent, and therefore a
sufficiently large layer is needed to reproduce bulk solvent conditions. If PBC
is being used, then in general the solvent structure and unit cell size should
be at least the length of the longest axis of the molecule plus twice the
cutoff distance from nonbond interactions so that no solvent water molecule
interacts with both the solute and its image. Bear in mind when using "long and
narrow" unit cells that some portions of the macromolecule may rotate "outside"
of the solvation shell during simulation.

Consideration on shapes
-----------------------

As mentioned in the tutorial section on [[FINAL The energy function|CHARMM's
energy function]], CHARMM supports PBC using any valid crystal unit cell (in
practice, however, only a few different shapes of unit cells are actually
used). It is convenient to set up the crystal structure (if one is needed)
while solvating the system.  Bear in mind that regular shapes must be used as
unit cells, *i.e.* you cannot have a spherically shaped unit cells because
spheres cannot be packed to completely fill a given volume of space,

CHARMMing, our Web based interface to CHARMM, contains functionality to
automatically determine an efficient shape. To do so, CHARMMing examines the
longest and shortest axis of the structure. If the longest axis is more than
30% longer or the shortest axis is more than 30% shorter when compared to the
middle axis, then a hexagonal prism is used. Otherwise, a rhombic dodecahedron
is chosen. The rhombic dodecahedron (RHDO) is a good choice for globular
proteins because it most closely approximates a sphere and thus is the most
efficient crystal shape (not a lot of excess water). The hexagonal prism works
well for long, thin structures for a similar reason, however it is necessary to
be careful to make sure that the macromolecule does not rotate outside of the
prism during the simulation. It is generally not a good idea to put restraints
on molecules undergoing in MD, however to keep long and narrow structures from
rotating MMFP cylinder restraints (see `mmfp.doc
<http://www.charmm.org/html/documentation/current/mmfp.html>`_ might be better
than simply restraining the end residues. Ideally, however, the molecule should
be able to rotate freely without moving outside the unit cell.

Performing the solvation
------------------------

Orienting the structure
***********************

Once you have the structure read into CHARMM, it is desirable to orient it with
the "COOR ORIEnt" command. This will rotate the molecule so that the x axis is
the longest axis, the y axis is the second longest, and the z axis is the
shortest. You can then figure out the lengths of each of the three axes via:

.. code-block:: chm

 coor stat
 calc xdist = abs( ?xmax - ?xmin )
 calc ydist = abs( ?ymax - ?ymin )
 calc zdist = abs( ?zmax - ?zmin )

Reading in and cutting down the water
*************************************

Once you have the macromolecule read in and oriented properly, the next thing
that is needed is a water box large enough to accommodate the desired crystal
structure. You can download the one used by CHARMMing as a CHARMM coordinate
file `here <http://www.charmmtutorial.org/static/files/water.crd>`_. Once you
have the file, you can read it in with your existing structure. To do so, you
must do two things: (1) append the waters to the PSF as a new segment and (2)
read their coordinates in. Since there are 46656 waters in the box you can do
this via the following commands:

.. code-block:: chm

 ! append to the PSF -- the generate command
 ! automatically appends the BWAT segment to
 ! your existing structure
 read sequence tip3 46656 
 generate bwat noangle nodihedral

 ! now read the coordinates in using APPEND
 ! to add them to the current set
 read coor card append name water.crd

When this is done, it is necessary to delete the water molecules that overlap
with the solute. In CHARMMing, we accomplish this by deleting all waters whose
oxygen atom is within 2.5 angstroms of the solute. The command to do this
(assuming the segment of bulk waters is name BWAT) is:

.. code-block:: chm

 delete atom sort sele .byres. ( segid BWAT .and. type oh2 .and. - 
   (( .not. (segid BWAT .or. hydrogen)) .around. 2.5 )) end

After deleting the atoms that overlap, we recommend cutting the water structure
down to a sphere just large enough to contain the final crystal structure. A
reasonable way to do this is to set the diameter of this sphere just large
enough to circumscribe a cube with an edge length  of the longest exis of the
macromolecule (@xdist) plus two times the padding distance. All waters outside
this diameter are then deleted. This is done for efficiency reasons; it will be
much faster to remove extraneous water surrounding the unit cell if most of the
extraneous waters have already been removed, The CHARMM commands to do so are:

.. code-block:: chm

 calc caxislen = ?xdim + ( 2 * @padding )
 calc caxislsq = @caxislen * @caxislen
 calc spherer = ( sqrt( 3 * @caxislsq ) ) / 2
 delete atom sort sele .byres. ( .not. ( point 0.0 0.0 0.0 cut @spherer ) .and. ( segid bwat ) ) end

In the above commands, *@caxislen* is set to the longest dimension of the
macromolecule plus two times the desired padding and *@caxislsq* is the square
of this value. The *@padding* is often set between 5 and 15 Å, depending on how
much of a solvation shield is needed around the molecule. The *@spherer*
variable is the radius of a sphere circumscribing the cube of edge length
*@caxislsq* (the requisite diameter is calculated via the Pythagorean Theorem
and divided by 2 to get the radius).  The complex part is the selection inside
the delete command. It select all residues that are not within a radius
*@spherer* of point (0,0,0) (this is another reason why it's important to do a
*COOR ORIEnt* before adding the water, as it will center the molecule at the
origin), and that are part of the BWAT segment (in the file from CHARMMing, the
waters have their segment ID set as BWAT, adjust as necessary for your own
structures). It is important to select by residue (*.BYRES.*) because otherwise
you might wind up deleting only part of a water molecule.

Building the crystal
********************

When we have the water structure cut down to size, we can go ahead and create
the unit cell. This is done in two steps: the crystal structure is defined and
then built. This tutorial assumes that you are using one of the standard shapes
built into CHARMM (cube, hexagon, tetragonal, rhombic dodecahedron, etc). To
define a crystal shape use the CRYStal DEFIne command. The basic syntax is:

.. code-block:: chm

 cryst defi <type> A B C α β γ


where a. b. and c are the edge length and α, β, and γ are the
angles. For example, to define a cube with sides of 20 angstroms you would
write:

.. code-block:: chm

 cryst defi cubic 20. 20. 20. 90. 90. 90.

to define a rhombic dodecahedron with edge length 30 you would write :

.. code-block:: chm

 cryst defi rhdo 30. 30. 30. 60. 90. 60.

The correct alpha, beta, and gamma values for each supported crystal type are
given in `crystl.doc
<http://www.charmm.org/documentation/current/crystl.html>`_. Note that a high
degree of precision for the angles is necessary to construct the truncated
octahedron structure (*i.e.* if you try to round off the angles, you will not get
a proper octahedron).

Once the crystal is defined, you can go ahead and build it with:

.. code-block:: chm

 cryst build noper 0

The *NOPERations* option specifies how many crystal operations need to be
performed. A regular shape centered at the origin with only translational
symmetry does not need any crystal operations (see the discussion on crystal
structure from the [[FINAL The Energy Function|energy page of this tutorial]].

Removing waters outside the crystal structure
*********************************************

It is now necessary to remove the waters that lie outside of the unit cell.
This can be done by setting up images and forcing an image update (via the
UPDAte command). When the image update is run, all atoms that are outside of
the unit cell boundaries will be moved back inside them. By copying the
coordinates to the comparison set before doing the update and finding the
difference between the main set (the coordinates after the update) and
comparison set (the coordinates after the update), we can detect which atoms
moved and are, therefore extraneous. We can then delete these so as to preserve
the correct density of the equilibrated water structure. This procedure is
shown in detail in the [[FINAL Full example|worked-out example]].

The solvation procedure is now complete.

After solvation
---------------

It is desirable to do a quick steepest-descent minimization (only a few tens of
steps) to remove bad van der Waals contacts.

When using the solvated structure with PBC in a different script, you will need
to recreate the crystal structure. To assist with this task, CHARMM can write
out the crystal transform file as follows:

.. code-block:: chm

 open unit 50 write card name crystal.xtl
 cryst write card unit 50
 * crystal structure -- it might be a good idea to put a, b, c, 
 * alpha, beta, gamma in the header
 *

This file can then be read back before the *CRYStal BUILd* command.
Alternatively, you can just re-run *CRYStal DEFIne* with the exact same a, b,
c, α, β, γ, and then re-run *CRYStal BUILd*. It is also desirable to put the
lattice type and dimensions into the PSF and coordinate files that are written
out after solvation.

Remember *CRYStal* by itself sets up the image transformations and everything
needed for the energy calculation. You need to use the *IMAGe* command to
specify which atoms in your system are subject to image centering.

