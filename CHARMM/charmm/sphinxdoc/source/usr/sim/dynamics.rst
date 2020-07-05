.. index:: DYNA, Dynamics

.. _usr-sim-dynamics:

Molecular Dynamics
==================

Molecular dynamics simulates the motion of a system by numerically integrating
Newton's second law of motion. The objective of a MD simulation is usually to
determine the time correlation between events or to get a statistically-valid
collection of structures that satisfies ergodicity (this is called sampling
from the ensemble). Currently, MD simulations tend to run from several
nanoseconds to around a microsecond (although millisecond scale MD simulations
are possible on advanced hardware). This is often too short to simulate
biological processes such as the folding of complex polypeptides. Biophysicists
often therefore look to use the results of MD simulations to make statistical
arguments about the behavior of a structure instead of directly observing long
scale behavior. This tutorial will focus on the practical issues of setting up
a molecular dynamics simulation in CHARMM. We will begin with a discussion of
the prerequisites for undertaking a successful simulation.

Prerequisites
-------------

In order to begin an MD simulation, you should have a structure that has been
reasonably well minimized under the **exact same** conditions as you are
planning to use for the dynamics. This means that the non-bonded set-up
(including Ewald) and periodic boundary conditions should be exactly the same.
However, in many cases it is undesirable to minimize the structure too much, as
it may deform in an undesirable manner. Please see the [[FINAL
Minimization|Minimization]] page for more details.

If there are problems with the way that the model was built, they will likely
manifest themselves in the dynamics run. In many cases, the potential energy of
the system will explode to unrealistic levels. In other cases, the structure
will twist into physically impossible conformations. You should keep an eye on
the energy and position of your system during dynamics to make sure that this
does not happen. Visualizing the MD trajectory is highly recommended as a check
for possible problems. Also, pay careful attention to the warnings CHARMM
prints out and resist the temptation to disable ECHECK (or set it to an
unreasonably high value). Exceeding the energy change tolerance (ECHECK) and
deviations in SHAKE can indicate model building errors or an incorrect dynamics
set-up. 

It is usually necessary to heat and equilibrate a structure at the desired
temperature after minimization. This topic is discussed further below. For
finicky structures, it might be necessary to begin the equilibration at a
shorter time step (e.g. 0.5 fs instead of 1 fs) to keep within the energy
change tolerance. It may also be desirable in certain circumstances to restrain
the structure during heating in order to prevent it from deviating too far from
the minimized structure. However, it is generally **not** a good idea to fix
atoms in place (using cons fix) during any type of dynamics. CHARMMing fixes
atoms whose parameters were estimated via GENRTF, but this is not a good idea
for real production dynamics (atoms without "real" parameters should either be
deleted or real parameters generated for them using, for example, the CHARMM
General Force Field).

Choosing an ensemble
--------------------

As mentioned above, one of the purposes of running molecular dynamics is to
sample a collection of structures so as to be able to make a valid statistical
argument about the system of interest. CHARMM supports several different types
of ensembles which you can sample from. A comparison of these ensembles is
beyond the scope of this tutorial; consult your favorite statistical mechanics
textbook for more details!

The main types of ensembles used are:

* Canonical (NVT): In this ensemble, number of particles (N), volume (V), and
  temperature (T) are held constant. 
* Microcanonical (NVE): This ensemble holds number of particles, volume, and
  total energy (E) constant.
* Isothermic-isobaric (NPT): For this ensemble number of  particles, pressure
  (P), and temperature are constant.

CHARMM has several temperature and pressure control mechanisms which are
discussed below. After initial thermal equilibration, an NVE ensemble can be
simulated by disabling any further heating (setting *IHTFRQ* and *IEQFRQ* to
0). This is discussed in further detail below.

A note about the DYNAmics command
---------------------------------

Molecular dynamics simulations in CHARMM are run via the *DYNAmics* command.
This command takes a lot of options, and for that reason the first six
characters of each subcommand name are significant (as opposed to the first
four elsewhere in CHARMM). All of the key words listed on this page of the
tutorial are options to the dynamics command. Some of the more important
options that will be used in almost every MD run are:

* NSTEP <integer>: specifies the number of steps to be run
* TIMESTep <real>: The time step in picoseconds, 0.001 is 1 fs.
* NPRINT <integer>: specifies the frequency at which the energy should be
  printed out, e.g. NPRINT 100 prints the energy every one hundred steps.
* NSAVC <integer>: frequency to write the coordinates to a trajectory file
* IUNCRD <integer>: Unit number of the coordinate trajectory file (must by
  opened before DYNAmics are invoked)
* NTRFRQ <integer>: Number of steps to check for and cancel external
  translation and rotation forces

A complete dynamics set-up might look like:

.. code-block:: chm

 ! open the trajectory file
 open unit 50 write unform name trajectory.trj

 ! IHTFRQ & IEQFRQ are, by default 0, so since
 ! we don't set them here, this runs with the
 ! NVE ensemble
 dyna leap start nstep 1000000 timestep 0.001 -
      nprint 1000 nsavc 50 iuncrd 50 ntrfrq 5000

Note: it is not usually correct to run NVE dynamics without prior heating and
equilibration (discussed below); this example merely shows the correct way to
invoke the dynamics command,

This will run 1 ns of molecular dynamics, saving a coordinate trajectory frame
to a file called trajectory.trj every 50 steps (0.05 picoseconds). It is also
possible to save the energy, temperature, and velocities to a trajectory, but
this is less used in practice. We will talk much more about trajectory files
and what to do with them once you have them as the tutorial progresses. The
only option that you have not seen here is *START*, which tells CHARMM to start
a new run (as opposed to continuing an old run). We'll talk about restarting
molecular dynamics simulations below.

DYNAmics, like other CHARMM commands, will use the previously established
non-bond configuration unless these are over-ridden in the *DYNAmics* command
itself. The recommended practice is to set up your non-bond options in advanced
so as not to clutter up the DYNA command (it has enough options already!).

What actually happens during molecular dynamics
-----------------------------------------------

Choosing an integrator
**********************

The integrator is the method used to numerically integrate Newton's second law
of motion during molecular dynamics. CHARMM supports five integrators:

* The leapfrog verlet integrator (keyword *LEAP*): this integrator is similar
  to the standard 3-step Verlet integrator, but provides additional accuracy.
* The original verlet integrator (keyword *ORIG*): this is the standard verlet
  integrator. In most cases, LEAP is preferred for its higher accuracy
* The velocity verlet integrator (keyword *VVER*): implements the Verlet
  algorithm differently. However, it does not print the Hamiltonian during
  dynamics so results validation is more difficult.
* The 4D verlet (*VER4*) integrator: this algorithm is primarily used with 4-D
  molecular dynamics, which is beyond the scope of this tutorial
* The velocity verlet 2 (*VV2*) integrator: this integrator is required for use
  with polarizable force fields using drude particles due to the fact that the
  drude particles must be integrated independently.

Most molecular dynamics simulations use the LEAP integrator (which supports
Langevin dynamics as well).

How the leapfrog Verlet algorithm work
**************************************

Given the atomic positions (:math:`X`) at timestep :math:`t` and the velocities
(:math:`V`) at :math:`t - \frac{1}{2}\Delta t`, the leapfrog verlet integrator
computes the positions at :math:`t + \Delta t` and the velocity at :math:`t +
\frac{1}{2}\Delta t` using the following procedure:


* Calculate the acceleration on each atom :math:`i` (:math:`a_{i}`), using the
  formula :math:`a_i = \frac{F_i}{m_i}` where :math:`F_i` is the force on atom
  :math:`i`, which is the negative gradient (first derivative) of [[FINAL The
  energy function|the energy function]] and :math:`m_{i}` is the mass of atom
  :math:`i`.
* Compute the velocity at :math:`t + \frac{1}{2}\Delta t` via the formula
  :math:`V_{t + \frac{1}{2} \Delta t, i} = V_{t - \frac{1}{2} \Delta t, i} +
  a_{i}`
* Compute the positions at :math:`t + \Delta t` via the equation :math:`X_{t +
  \Delta t, i} = X_{t, i} + \Delta t V_{t + \frac{1}{2}\Delta t, i}`

In this case :math:`\Delta t` represents the time step specified by the users;
in general, the larger the time step, the less accurate the numerical
integration will be. We can see why the algorithm is called a leap-frog since
the velocities are computed at midpoints between the time steps.

Further discussion of the leapfrog Verlet and other molecular dynamics
algorithms with formulas may be found on the `embnet Theory of Molecular
Dynamics tutorial <http://www.ch.embnet.org/MD_tutorial/pages/MD.Part1.html>`_.

Velocity assignment
-------------------

If you study the procedure for leapfrog Verlet outlined above, you will notice
a potential problem: for the first time step in the simulation, we need to have
velocities for all of the atoms. Fortunately, CHARMM has the ability to
generate initial velocities for the system or to read them from the COMP
coordinate set. In most cases, CHARMM generates initial velocities from a
distribution. The distribution is determined by the value of *IASVEL*; if this
is greater than 0, a gaussian distribution is used, if it is less than 0 a
uniform distribution is used, and if it is equal to 0, initial velocities are
read from the COMP coordinate set. Subsequent velocity assignment and rescaling
(during heating and equilibration) are controlled by the *IASORS* and *ISCVEL*
options to *DYNAmics*. When *IASORS* is 0, velocities are rescaled, otherwise
they are reassigned. In the latter case, the method of assignment is determined
by the *IASVEL* option, which works the same way as it does for initial
velocity assignment. When velocities are rescaled, the rescaling mechanism is
controlled by the *ISCVEL* option; when it is 0, a single scale factor is used
for all atoms otherwise the scale factor for each atom is dependent on the
ratio of the average kinetic energy along each degree of freedom of the atom.
Note that velocities will only be assigned or rescaled at the start of a
dynamics run so long as *IHTFRQ*, *IEQFRQ*, and *ICHECW* are all set to 0
(these options are described below).

Temperature and pressure control
--------------------------------

Primitive temperature control
*****************************

CHARMM employs several methods of controlling temperature and pressure, however
not all of these are equally good! A primitive method of temperature control
can be obtained by using the basic temperature control subcommands. These are
*IHTFRQ*, *IEQFRQ*, *FIRSTT*, *FINALT*, *TBATH*, *TEMINC*, and *ICHECW*. The
*FIRSTT* and *FINALT* option set the initial and ending temperatures of the
simulation (if the temperature of the simulation is expected to remain constant
then these should be the same. The *TBATH* sets the temperature of the external
heat bath that is coupled to the simulation. In general it is a good idea to
set this to *FINALT*. *TEMINC* and *IHTFRQ* control heating and should be
omitted in constant temperature simulations. During a heating run, the
temperatures are increased by *TEMINC* degrees every *IHTFRQ* steps, which is
accomplished by either reassigning or rescaling (depending on the value of
*IASORS* -- see above) the velocities accordingly. *IEQFRQ* behaves similarly
when equilibrating a simulation at constant temperature; velocities are
reassigned or rescaled to match the desired temperature every *IEQFRQ* steps.
The primary difference with *IEQFRQ* is that velocity adjustments are not tied
to the heating frequency (*IHTFRQ*) or the temperature increment (*TEMINC*), It
is possible to further control this rescaling using *ICHECW*; if *ICHECW* is 0,
the velocities will always be rescaled, but if if it is 1 they will only be
rescaled if they are outside a range around the desired temperature. The upper
and lower bounds of this range are set by the *TWINDH* and *TWINDL* options.

CPT dynamics and the Hoover Thermostat
**************************************

It is important to note that the primitive temperature control above is
**not** sufficient to ensure a NVT statistical ensemble, however it is often
good enough for thermal heating and equilibration prior to an NVE dynamics run.
For proper NVT dynamics, however, some sort of thermostat must be used. The
best thermostat available in CHARMM is the Hoover thermostat, which is part of
the constant pressure and temperature (CPT) method, which is enabled by the
*CPT* option to the dynamics command. The hoover thermostat often is used with
constant pressure (the *PCONS* keyword), but this is not a requirement. NVT, as
opposed to NPT, can be run even when *PCONS* is used (see the next paragraph).
The reason to do this is to make CHARMM print out pressure statistics (these
will not be displayed if *PCONS* is not used). The weak-coupling Berendsen
thermostat (the *TCONS* option) has a number of known problems and its use is
not recommended '''under any circumstances'''. Both the constant pressure and
temperature methods in CPT use the Langevin piston method.

The options to use when applying the hoover thermostat and constant pressure
are:

.. code-block:: chm

 dyna cpt pcons pref AAA pgamma BBB pmass CCC hoover reft XXX tmass YYY ...

These options can be divided into two parts, those controlling temperature, and
those controlling pressure:

Temperature Options
*******************

* *HOOVER*: tells CHARMM to use the Hoover thermostat
* *REFT*: the temperature at which the thermostat is set, i.e. this is the
  temperature at which the thermostat will keep the system.
* *TMASS*: the mass of the temperature piston. The *TMASS* option is not in
  atomic mass units but :math:`\frac{kcal * ps^2}{mol}`. We have found that a
  good value for the piston mass is generally about 20% of the system's mass
  (which can be found in the *?STOT* variable after performing the *SCALar MASS
  STAT* command), and this is the value that CHARMMing uses.

Pressure Options
****************

* *PCONS*: Turns on pressure reporting. The pressure of the system will be held
  constant unless *PMASS* is set to 0 (in this case, pressure statistics are
  still displayed, but constant pressure is not enforced and an NVT ensemble is
  obtained).
* *PREF*: The reference pressure in atmospheres. It is usually set to 1 to
  mimic biological conditions.
* *PMASS*: The mass of the pressure piston, which **is** measured in atomic
  mass units (amu). This value can be set to 0 which will disable the pressure
  control (but pressure will still be reported in the output). In cases when
  *PMASS* is non-zero, a good value is generally 2% of the system mass. 
* *PGAMMA*: The Langevin collision frequency for the pressure piston. This is
  usually set to 0 except during heating to allow the pressure piston to move
  freely. During heating it is desirable to damp the piston (to prevent the
  volume from changing too rapidly) by setting a positive *PGAMMA*.

Writing out trajectories
------------------------

CHARMM provides a method to save coordinates and velocities of all atoms at
fixed intervals during a simulation. In most cases it is not necessary to store
velocities, so most users only save coordinate trajectories. The energies and
forces for a given time step can be derived from the coordinate trajectory
frame. These so-called trajectory files can be used for later analysis, a topic
which is discussed further in the analysis section of this tutorial. In order
to write out trajectory files, it is necessary to open a unit for writing
unformatted (binary) data. This can be done with:

.. code-block:: chm

 open unit 20 write unform name velocities.trj
 open unit 21 write unform name coordinates.trj

The *IUNCRD* and *IUNVEL* options then specify which units are used for
coordinate and velocities trajectories. The *NSAVC* and *NSAVV* options
determine the frequency at which coordinates and velocities are saved,
respectively. Therefore, to write velocities and coordinates out every 100
steps, you would do:

.. code-block:: chm

 dyna ... iuncrd 21 iunvel 20 nsavc 100 nsavv 100 ...

Some people like to give their trajectory files .dcd extensions because certain
third-party software automatically recognizes these as CHARMM-format
trejectories, however CHARMM itself does not care what you name the files. One
thing to be concerned about is trajectory portability between 32 and 64 bit
machines. This is not a problem with newer versions of CHARMM (generally
defined as c34 and later), but trajectories generated by older versions
compiled with certain versions of gcc might only be readable on the same type
of computer as the one that generated them. Trajectory I/O issues caused by
size mismatches can often be dealt with using the *DYNA FORMAT* command. Please
refer to `dynamc.doc
<http://www.charmm.org/documentation/current/dynamc.html>`_ for further
details.

Re-starting MD runs
-------------------

Because molecular dynamics runs can take a long time (potentially months of
wall clock time), CHARMM provides the functionality to stop and restart runs
via the use of a restart file. The restart file contains the coordinates,
velocities, and other state information. If you have a restart file, you can
restart a dynamics run by using "*DYNA RESTart ...*" rather than "*DYNA START
...*". Please note that dynamics parameters (number of steps, trajectory units,
etc.) are not saved in the restart file, and therefore you need to write out
the entire *DYNAmics* command again (or better yet, paste it from your previous
script and modify as necessary). If the *RESTart* options is used than the
*IUNREAd* option must be given to tell CHARMM which unit to read the restart
file from. To write out  restart file, it is necessary to open the unit (card
format) and use the *IUNWRIte* option to tell CHARMM which unit to use. You can
use the *ISVFRQ* option to determine how often the re-start file will be saved.
However, CHARMM will use the same file name for each restart file so if the
program crashes while writing the restart file, the file will be corrupted (you
can write a script to take snapshots of the restart file, however). Many users
prefer to write one restart file per run and utilize many short runs rather
than a single long one.

Here is an example of starting a run and writing out a restart file:

.. code-block:: chm

 open unit 20 write card name restart.res
 open unit 50 write unform name run1.dcd

 dyna start leap cpt -
   nstep 1000000 timestep 0.001 -       ! 1 ns in 1 fs timesteps
   pcons pmass 0.0 pgamma 0 pref 1.0 -  ! report pressure, but do NVT
   hoover reft 300.0 tmass 500.0 -      ! constant temperature
   iuncrd 50 iunwri 20 nsavc 1000 -     ! trajectory and restart
   iasors 1 iasvel 1 ihtfrq 0 -         ! assign velocities, but
   ieqfrq 0 ichecw 0                    ! do no rescaling (NVT)

This command will produce a restart file after one nanosecond of dynamics. To
continue the run, you would write:

.. code-block:: chm

 open unit 20 read card name restart.res
 open unit 21 write card name restart2.res
 open unit 50 write unform name run2.dcd

 dyna restart leap cpt -
   iunrea 20 -                         ! unit to read restart file
   nstep 1000000 timestep 0.001 -      ! 1 ns in 1 fs timesteps
   pcons pmass 0.0 pgamma 0 pref 1.0 - ! constant pressure
   hoover reft 300.0 tmass 500.0 -     ! constant temperature
   iuncrd 50 iunwri 21 nsavc 1000 -    ! trajectory and restart
   iasors 1 iasvel 1 ihtfrq 0 -        ! assign velocities,
   ieqfrq 0 ichecw 0                   ! but do no rescaling

Note that it is not possible to append to the trajectory file given by the
first dynamics command. Instead, you must create a second trajectory and then
merge it with the first one (alternatively, CHARMM's analysis tools have very
good support for reading a sequence of trajectory files -- this process is
described in the analysis section of the tutorial). Also, although *IASVEL* is
set to 1, no velocity re-assignment will be done because the velocities are
read from the restart file and  *IHTFRQ*, *IEQFRQ*, and *ICHECW* are all 0.

Summary and recommendations
---------------------------

Although the *DYNAmics* command is complex, most of its options can be broken
down into simple groups controlling the number of steps and the size of each
step, temperature and pressure control, setting up of initial velocities, and
and reading/writing various control and trajectory files. We recommend setting
up non-bond parameters before the main *DYNA* commands. CHARMM's `dynamc.doc
<http://www.charmm.org/documentation/current/dynamc.html>`_ provides
recommendations of what options to use for various purposes. These are
reproduced below with some minor fixes to represent currently suggested best
practices (in particular #3 has been changed to use the Hoover thermostat, as
we do not recommend using the weak-coupling Berendsen thermostat):

**1)    For heating and early equilibration:**

.. code-block:: chm

 DYNAMICS LEAP VERLET RESTART(*)  NSTEP 20000 TIMESTEP 0.001 -
     IPRFRQ 1000 IHTFRQ 100 IEQFRQ 5000 NTRFRQ 5000  -
     IUNREA 30 IUNWRI 31 IUNCRD 50 IUNVEL -1 KUNIT 70 -
     NPRINT 100 NSAVC 100 NSAVV 0 -
     FIRSTT 100.0 FINALT 300.0 TEMINC 10.0   -
     IASORS 1 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 10.0 TWINDL -10.0

 (*)   Except for first run, then use STRT in place of RESTART


**2)    For late equilibration and analysis runs:**

.. code-block:: chm

 DYNAMICS LEAP VERLET RESTART  NSTEP 20000 TIMESTEP 0.001 -
     IPRFRQ 1000 IHTFRQ 0 IEQFRQ 0(*) NTRFRQ 5000  -
     IUNREA 30 IUNWRI 31 IUNCRD 50 IUNVEL -1 KUNIT 70 -
     NPRINT 100 NSAVC 100 NSAVV 0 -
     FIRSTT 100.0 FINALT 300.0 -
     IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 1 TWINDH 10.0 TWINDL -10.0

 (*) This should probably be positive for an equilibration run.

**3)  For constant temperature and/or pressure dynamics**

.. code-block:: chm

  DYNA LEAP VERLET STRT(*) CPT NSTEP 20000 TIMESTEP 0.001 -
     IPRFRQ 1000 IHTFRQ 0 IEQFRQ 0 NTRFRQ 0  -
     IUNREA 30 IUNWRI 31 IUNCRD 50 IUNVEL -1 KUNIT 70 -
     NPRINT 100 NSAVC 100 NSAVV 0 IHBFRQ 0 INBFRQ 25  -
     PCONS PMASS(+) 50.0 PGAMMA 0 PREF 1.0 -
     HOOVER REFT 300.0 TMASS(+) 500.0 - 
     IASORS 0 IASVEL 1 ISCVEL 0 ICHECW 0 TWINDH 0.0 TWINDL 0.0

 (*)   For first run, use RESTART otherwise
 (+)   PMASS and TMASS must be adjusted to match your system size. PMASS should be set to 0 if NVT, instead of NPT, is desired. See above for more information about this.

Langevin Dynamics
=================

Standard Langevin Dynamics
--------------------------

Langevin dynamics (LD) is similar to molecular dynamics (MD), except that
instead of numerically integrating over Newton's second law of motion (the
familiar :math:`F = ma`), the Langevin equation is used instead. This equation
may be written as::

:math:`F_{i} = m_{i}a_{i} = - \nabla V_{i} - \gamma_{i} + R(t)`

In the above equation, :math:`- \nabla V` is the negative gradient of the
system (i.e. the force), :math:`\gamma` represents the frictional drag on the
system, and :math:`R(t)` represents jostling forces at time :math:`t`. These
random forces have an expected value of 0 and are uncorrelated between time
steps. In CHARMM, the random force is specified via setting a collision
frequency using the *SCALar* command to set the *FBETa* value,
*e.g.*:

.. code-block:: chm

 scalar fbeta set 60.0

sets the collision frequency of 60 per picosecond. This is approximately the
viscosity of water at room temperature.

The same sorts of integrators that we saw earlier can be used to perform
Langevin Dynamics, although in CHARMM only the original Verlet and leap-frog
Verlet ones support this feature.

In Langevin dynamics, the system is coupled to an external heat bath, providing
fairly constant temperature throughout the simulation. The temperature of the
external bath may be specified via the TBATh option to the DYNAmics command.

Self-guided Langevin dynamics
-----------------------------

In self-guided Langevin dynamics (SGLD), the standard LD equation is used, but
an additional guiding force, which is a momentum-derived force that pushes
atoms, is added. Because of this factor, more conformations will be sampled in
a given interval and thus it should only be used if determining the timescale
in which events occur is not important. However, the current working hypothesis
is that the **order** that events occur is unchanged with SGLD (*i.e.* only the
timescale is distorted). The original SGLD implementation does not preserve the
canonical ensemble, however re-weighting may be used to calculate canonical
ensemble properties from an SGLD trajectory. A new SGLDfp method, which
preserves the ensemble at the expense of some sampling, is also available.

Considerations on using Langevin Dynamics
-----------------------------------------

The DYNAmics command is used to run both Langevin dynamics and molecular
dynamics, so once you have learned how to set up a molecular dynamics
simulation, it should be fairly straightforward to run LD as well. However,
there are a few differences that must be noted.

* The FBETA scalar value must be set (see above). If this is not done, a
  default value of 0 will be used and you will effectively be doing regular MD
  (albeit with the extra "kicking" forces). CHARMM did not used to warn you
  about this, but more recent versions do. 
* An *FBETA* value of 1 or less (significantly lower viscosity than
  water) is often used for explicitly solvated systems since the presence of
  water molecules will provide all of the desired frictional effects. A value
  of 2 is often used for vacuum systems which are not otherwise restrained by
  the presence of solvent. One might reasonably ask why anyone would want to
  run LD on a solvated system. In many cases, LD is used for its temperature
  coupling. Also, although Langevin Dynamics mimics the frictional effects of
  water, other properties of a solvated environment (e.g. charge screening for
  electrostatics) are not present.
* Since Langevin dynamics has its own temperature control, it should not be
  used with the constant pressure and temperature (CPT) keyword. Likewise, all
  CPT specific features (such as PCONS and the HOOVER commands) must be removed
  from the DYNAmics command. Similarly, IHTFRQ and IEQFRQ should be set to 0.
* Usually, Langevin dynamics uses the leap-frog integrator (key word LEAP), but
  ORIG will work as well.

To initiate Langevin dynamics, use the LANGevin keyword. The command usually
looks something like :

.. code-block:: chm

 dyna leap langevin tbath XXX ...

where XXX is the desired bath temperature.

