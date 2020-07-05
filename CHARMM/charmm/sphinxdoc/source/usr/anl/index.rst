.. _usr-anl-index:

.. index:: Analysis

Data Analysis
=============

Trajectories have to be processed before a meaningful analysis is
possible. It is at this step that (intimate) knowledge of statistical mechanics
is most needed; one might argue that your command of statistical mechanics (and
statistics in general) is the limiting factor for your analyses (note: stat.
mech. provides the bridge between microscopic quantities you can derive from
the trajectories and macroscopic quantities you can deduce from experiments).
In as much as there is no limit to the analyses you can carry out, there can be
no guide covering all aspects. There is, however, one type of analysis that is
carried out after (or during) almost all MD simulations: monitoring how much
the structure of your protein changes compared to the starting structure (be
that the original pdb structure, or the structure after some minimization /
equilibration): In the following, we'll focus on this task and use it to
present some of the most important analysis facilities of CHARMM. In
particular, we'll compute the root mean square deviation (RMSD) from the
starting structure and monitor the solvent accessible surface area (SASA) and
the radius of gyration. As a more advanced example, we'll look at the
conservation of secondary structure elements (only recent versions of CHARMM
can do this with built-in tools, so we can use this example to illustrate how
external tools can be used to aid in the analysis of trajectories (if no
suitable built in tool exists!).

Trajectories are not the only useful source of information about a run; many
properties such as average potential or kinetic energy can be computed from the
output file. Rick Venable has written a useful script for doing this, which may
be found `in this CHARMM forum post
<http://www.charmm.org/ubbthreads-7-5-5/ubbthreads.php?ubb=showflat&Number=24462>`_.

On trajectory organization
--------------------------

In this section, we try to summarize things to know that are common to all
(most) analyses in CHARMM. 

A production MD simulation of a solvated protein can require weeks and months
of computer time. Obviously, one doesn't trust the computer not to crash during
this time, so, usually, the calculation is split into shorter pieces by means
of restart files. Thus, at the end of the simulation, one doesn't have a single
trajectory, but ten or fifty or even more trajectories. (Trajectories contain
vast amounts of data, so even if one had an "ultra-stable" computer with a
guaranteed uptime of 6 months, one most likely would have to split the
computation since a full trajectory might easily exceed the capacity of some
widely used file systems!)  The analysis commands in CHARMM can handle such
series of trajectories and actually check whether pieces are processed in
chronological order. For CHARMM version c35 and earlier, the maximum
**guaranteed** number of files (including standard input and output) that can
be open at any one time is 99. For the newer, Fortran 95 versions (c36 and up,
which are as of this writing not available to the general public) the limit is
whatever the maximum number of open files is for given process as imposed by
the operating system. Handling a series of trajectories is done as follows.
First, as with all other I/O, trajectories need to be opened
first

.. code-block:: chm

 OPEN UNIT 11 NAME system.1.dcd READ UNFO

with normal syntax; remember, however, that they are binary files, so you need
keyword UNFOrmatted or FILE (instead of FORMatted / CARD) for most other files.
If you have only one trajectory, the unit number is (mostly) up to you (under
Unix/Linux, you may want to avoid unit numbers 0, 1, 5, 6!). Suppose you have
10 trajectories that you want to process by some command later in one sweep
(i.e., the 10 trajectories contain chronological coordinate information about
your system). You then have to open the ten files with continuous unit numbers.
Let's say that your first unit is 11, then your sequence of OPEN statements
should look like:

.. code-block:: chm

 OPEN UNIT 11 NAME system.1.dcd READ UNFO 
 OPEN UNIT 12 NAME system.2.dcd READ UNFO 
 OPEN UNIT 13 NAME system.3.dcd READ UNFO 
 ...
 OPEN UNIT 19 NAME system.9.dcd READ UNFO 
 OPEN UNIT 20 NAME system.10.dcd READ UNFO 

Let command be a CHARMM command to process trajectories. You then process your
trajectories by a statement like:

.. code-block:: chm

 command FIRStu 11 NUNit 10 <other options to command>

The option FIRStu tells *command* at which unit number to find the first
trajectory file; NUNit provides the information about how many trajectories
follow, *i.e.*, we tell *command* to read 10 trajectories found at unit numbers
11, 12, ... 20. Command now reads all frames from the 10 trajectories in order.
Upon switching to a new trajectory file, CHARMM parses header information in
the trajectories to carry out some sanity checks. In particular, CHARMM expects
coordinate frames in chronological order, since many properties one is
interested in are time-dependent. This traps, *e.g.*, the following mistake:
Suppose you confused the order of trajectories, as in:

.. code-block:: chm

 ! BAD EXAMPLE
 OPEN UNIT 11 NAME system.1.dcd READ UNFO 
 OPEN UNIT 12 NAME system.3.dcd READ UNFO ! ERROR: example of mistake
 OPEN UNIT 13 NAME system.2.dcd READ UNFO 
 ...

 command FIRStu 11 NUNit 10 <other options to command>

CHARMM will process the first trajectory. When switching to unit 12, connected
to system.3.dcd, CHARMM will note that the first coordinate set corresponds not
the one it expects following the end of system.1.dcd, ... and you'll crash!

Unit numbers are actually a fairly scarce resource. Once you are done with the
trajectories, you should, therefore, CLOSe all trajetories! BTW, OPENing and
CLOSing of 10 or 50 trajectories one at a time is not fun. You should quickly
consider writing a loop, which is described in the [[FINAL Basic CHARMM
Scripting|basic CHARMM scripting]] section.

Going back for one step, we should say a few words about the prerequisites for
doing a trajectory analysis. The general rule of the thumb is: Set up your
system ''exactly'' as you set it up during the generation of the trajectories.
In particular, read the same PSF and set the same energy options. If you used
PBC and you're performing an analysis that uses crystal and image properties,
set these up with CRYSTAL/IMAGe as during MD. For good measure you may want to
add SHAKe and other restraints you had during the MD.  There are exceptions and
border line cases where this doesn't apply, but in 90% of the cases this should
work.  Thus, an analysis run has typically a structure like the following:

.. code-block:: chm

 ! read rtf
 ! read params
 ! read PSF
 ! read some starting coordinates (necessary when using PBC (CRYStal),
        ! otherwise optional

 CRYS DEFI <as during MD>
 CRYS BUIL
 IMAG BYSE ...
 IMAG BYRE ...

 ENER <all nonbonded options> 

 ! start analysis

 ! Let there be 10 trajectories 
 OPEN UNIT 11 NAME system.1.dcd READ UNFO 
 OPEN UNIT 12 NAME system.2.dcd READ UNFO 
 OPEN UNIT 13 NAME system.3.dcd READ UNFO 
 ! ... 7 more OPEN statements

 command FIRStu 11 NUNit 10 <other options to command>

 CLOSe UNIT 11
 CLOSe Unit 12
 ! ... 8 more CLOSe statements

Carefully onsulting the documentation of the analysis commands that you want to
use is highly recommended. Pay particular attention to whether CRYSTal and
IMAGEs need to be set up for that particular analysis and if so be sure to set
them up exactly as they were for MD. If these are not needed then you can save
some time by not setting them up, particularly if you are running lots of
analyses back to back.

Some basic analyses
-------------------

CORREL
******

CORREL is a generalized facility for analyzing trajectories. The subsystem is
invoked with the CORREL command, which may take several arguments. It is
usually necessary to set the maximum number of atoms (MAXAtoms), and it may be
necessary to also adjust the maximum number of timesteps (MAXTimesteps) and
series (MAXSeries) if the defaults are insufficient. However, settings these
values too high will result in lots of memory being allocated, possibly leading
to a crash if the system does not have sufficient RAM.

What CORREL can track
*********************

One of the main concepts used by CORREL is the concept of the time series. A
time series describes the value of a property such as a bond length, the total
energy of the system, or the distance between two atoms, over the course of the
simulation. Time series are created via the ENTEr command, which takes a number
of different options depending on what time series data is needed. The basic
syntax for the command is:

.. code-block:: chm

 ENTEr <name> <type of data> -
   <atom selection and subcommands>

Each time series must be given a unique name. Then the user must specify the
type of data that is to be collected. An exhaustive listing of supported data
types is in `correl.doc
<http://www.charmm.org/documentation/current/correl.html>`_, but some commonly
used ones are:

* BOND for the length between two bonded atoms
* ANGLe for the angle between three atoms
* DIHEdral and IMPRoper
* RMS, for the root mean squared deviation of the system
* ENERgy, for total system energy
* SDIP for the dipole moment of the solvent shell

The atom selection and subcommands differ based on the type of data that is
needed, *e.g.* DIHEdral requires a selection of four atoms while ENERgy requires
no further arguments. Consulting the documentation is recommended as a number
of data types have different subcommands (*e.g.* whether or not to mass weight
the RMSD).

Once all desired time series are set up via ENTEr, it is necessary to tell
CORREL which trajectory or trajectories will be used to populate the time
series data. This is done via the TRAJectory command. The basic syntax is:

.. code-block:: chm

 TRAJectory FIRSTu <x> NUNIt <y> SKIP <skip> -
   BEGIn <start step> STOP <end step>

Where <x> is the first unit in the series of trajectories and <y> is the number
of units to read; following the previous example x is 11 and y is 10. This is
why it is important that trajectory unit numbers be consecutive and in order!

Once the TRAJectory command is issued, the time series data will be populated
and can be further used.

Display and Manipulation of Trajectory Data
-------------------------------------------

Once the time series are generated, they can be written out to a file or
manipulated further. To write out trajectory data, simply open a unit and tell
CHARMM to write the data:

.. code-block:: chm

 correl maxt X
 enter psi torsion segid resid atomtype1 segid resid atomtype2 segid resid atomtype3 segid resid atomtype4 geometry
 traj firstu 10 nunit 1 begin 100 stop 21000000 skip 10000
 write psi card unit 21
 * psi of residue 1
 *

Obviously, you must replace segid, resid, and the atom numbers with the correct
segments, residues, and atom types for your system.

The resulting data file may be viewed or plotted in an external program such as
GNUPlot or XMGrace.

CHARMM also containes facilities for manipulating the data, via the MANTIME
command.

Useful properties that can be calculated
----------------------------------------

Based only on the limited work with CORREL that has been presented so far,
several useful properties may be computed from simulation.

Mean energy
***********

The average energy, denoted <E> or :math:`\bar{E}` is just the average energy
value over the course of the simulation.

:math:`\bar{E} = \frac{1}{n} \sum_{i=1}^n E_{i}`

This can be calculated via extracting the energy at each time step from the
trajectory and averaging them. However, a more precise way of calculating this
value is to look in the CHARMM output file for lines beginning with "DYNA
AVER>". Usually, CHARMM prints out the average values (including the energy)
every 1000 steps (unless you tell CHARMM to behave differently -- see the
[[FINAL MD|molecular dynamics]] page for details). These average values are
only for the past 1000 (or however many) steps. Using the "DYNA AVER" values
will give you a much greater sample size, unless you are saving out every frame
into the trajectory file.

The average energy is of primary interest in NVT simulations to determine how
much total energy fluctuates at a given temperature, which is useful for
determining various statistical mechanics properties. However, it is also
interesting in NVE simulations to see the proportion of the total, constant
energy that is tied up as potential energy.

Mean structure
**************

The mean structure is just the average coordinate of each molecule, which may
be found in a similar manner as the mean energy. The only difference is that
the average structure can be mass-weighted.

Root mean square fluctuation
****************************

The root mean squared fluctuation is just the average difference between all
particles and their average position for a given time step.

Radius of gyration
******************

The radius of gyration measures the average distance between an atom and its
center of mass at a given time step. It is used as another measure of how much
atoms move around during a simulation.

Next steps
**********

An example of analyzing a molecular dynamics trajectory is given in the [[FINAL
Full example|full example]]. It will show you how to use *CORREL* to graph
various properties. We will also introduce more advanced features, such as time
correlation functions and manipulating trajectories. However, a detailed
discussion of how to analyze a simulation is beyond the scope of this tutorial.

