.. _usr-minimization:

.. index:: Minimization

Minimization
============

Energy minimization, as the name implies, is a procedure that attempts to
minimize the `energy <Energy calculation>` of the system to
the lowest possible point. This can be a challenging problem, as there is only
one global minimum of the energy surface, but many local minima. If steps are
not taken to avoid it, a minimization algorithm might get stuck in a local
minimum without ever finding the global minimum. In most cases, this is what
will actually happen in CHARMM, but as explained below, finding a global
minimum is often not necessary. CHARMM's minimization algorithms examine the
first (and in some cases second) derivatives to determine whether they are at a
minimum. More complex methods of exploring the energy surface (*e.g.*
conformational space annealing, or methods in which the system is repeatedly
heated and cooled) are beyond the scope of this tutorial.

Motivation
----------

One might reasonably ask why you would want to perform energy minimization in
the first place. There are a number of reasons, but they mostly center around
removing nonphysical contacts / interactions. For example, with a structure
that has been solved via X-ray crystallography, contacts with neighbors in the
crystal can cause changes from the *in vitro* structure. Missing coordinates
obtained from the internal coordinate facility or the HBUIld module of CHARMM
may be far from optimal.  Additionally, when two sets of coordinates are merged
(*e.g.*, when a protein is put inside a water box or sphere) it is possible that
there are steric clashes present in the resulting coordinate set. Such *'bad
contacts* may cause large, unphysical changes in potential energy in
subsequent dynamics simulations, or possibly other artifacts. A brief
minimization can remove such potential problems. As a result, minimization is
generally not done as an end unto itself, but to prepare for another type of
calculation and must be done under the same conditions as the later
calculation. For example, when planning a QM/MM calculation (*i.e.* when one part
of your system is modeled with quantum mechanics while the rest uses the
standard CHARMM force field)), it is important to re-minimize the system with
the QM/MM conditions in place.

How much should a structure be minimized?
-----------------------------------------

A minimization is considered converged if the root mean squared deviation of
the gradient (GRMS) is very close to zero (recall that the necessary condition
for a function to have a local minimum  is that all its first derivatives are
zero -- of course, it could also be an extreme, saddle point *etc.*). In
addition, the energy changes from step to step should be very small when the
minimization is close to convergence.

How much minimization is necessary depends upon the intended use of the final
structure. In general, if the output is to be used for solvation or dynamics,
it is not necessary for the minimization to be fully converged. In fact, when
preparing a structure for dynamics, it is not advised to let it move too much
from its original conformation. Exploration of conformational changes should be
done during dynamics itself, not minimization. Over-minimization can lead to
unphysical "freezing" of the structure. It may be desirable to use thermal
B-factors to restrain the structure (as these are hints from the people who
solved the structure about the certainty of the atomic positions). However, if
the structure is to be used for a normal modes calculation, then it should be
as close to a local minimum as possible. This is because a normal mode
calculation treats the potential as a harmonic well. If the structure is not
at, or very close to, a minimum, it will not be at the bottom of that well.
Consequently there will be elements of the matrix made up of second derivatives
of the energy (*i.e.* the Hessian) that will be negative. This this matrix will
not be positive definite and there will be negative eigenvalues and, hence,
imaginary normal modes.

When both minimizing and solvating a structure, it is often best to minimize
the structure in vacuum before performing solvation.It is often difficult to
resolve bad protein-protein contacts and bad protein-solvent contacts within
the same minimization, so it is often desirable to do these minimizations
separately.

What algorithm should be used?
------------------------------

A full description of the minimization options available in CHARMM can be found
in `minimiz.doc <http://www.charmm.org/documentation/current/minimiz.html>`_
which also documents the exact syntax for the MINImize command (described in
the next section). Most people use just two of the available  minimizers,
steepest descent (SD) and the adopted basis Newton Raphson method (ABNR). SD is
the simplest algorithm; it simply moves the coordinates in the negative
direction of the gradient. The only adjustable parameter is the step size,
which determines how much the coordinates are shifted at each step. The step
size is adjusted dynamically to achieve rapid convergence. The main problem
with steepest descent is that it does not generally converge to a local
minimum; however, it will rapidly improve the conformation when the GRMS is
high (*i.e.* when the system is far from a minimum). Therefore, it is often used
briefly on a new structure to quickly remove bad contacts and clashes. Many
practitioners ponly use SD for the initial 25-50 steps of minimization to
remove bad van der Waals contact, however it can be run for longer. Once this
is done, it is advisable to switch to a more precise minimizer, usually ABNR,
to finish the calculation. The ABNR minimizer will be able to detect when it
has converged and will exit automatically, therefore there is no danger in
running it for a long time (*i.e.*, specifying a large number of steps; see
below). However, as stated previously, except for when computing normal nodes
or preparing for QM/MM dynamics, there is generally no need to arrive at an
exact minimum. When the energy change from step to step is less than 0.001
kcal/mol, the structure is sufficiently minimized for most purposes.

In certain special cases, it may be necessary to use another minimizer such as
the *CONJugate-gradient* method or the *POWEll* minimizer.
Consulting the documentation carefully is recommended in these cases. In
particular, note that the *POWEll* minimizer does not support the
*INBFrq* and *IMGFrq* keywords (which determine how often the
non-bond and image atom lists are regenerated). Therefore, it will certainly
fail if the structure moves too much during the calculation. It is suggested to
use the *POWEll* method in a loop; however it is probably best avoided
alltogether.

How to perform minimization
---------------------------

Energy minimization can be performed once a valid PSF and coordinate set have
been loaded into CHARMM. The *MINImize* command is used to initiate the
procedure. Immediately after the command itself the method (e.g. *SD*,
*ABNR*) must be specified, followed by any minimizer specific options.
In general, the only options that need to be given are *NSTEp*,
TOLGradient, and *TOLEnergy*. The *NSTEp* option specifies the
maximum number of steps that the minimizer will run. Likewise, the
*TOLGradient* option tells the minimizer to exit once the specified GRMS
is achieved and the *TOLEnergy* option tells the routine to exit after
the energy changes by less than the given amount. If these are not specified,
minimization will run until its own convergence criteria is met. For most
systems, it is sufficient to do a few hundred steps of minimization with no
*TOLG* or *TOLE* specified to get the system away from sharp
peaks in the energy curve, followed by a sufficient number of ABNR steps to
meet the desired *TOLG*. For example:

.. code-block:: chm

 ! perform a basic energy minimization
 mini sd nstep 50
 mini abnr nstep 1000000 tolg 0.01

will be sufficient in most instances. As mentioned above, there is no harm in
setting the large value of *NSTEp*, even without the *TOLG* or
*TOLE* options; the minimizer will exit once it has converged (note,
however, that the defaults for *TOLE* and *TOLG* are 0.00 so
unless these are increased it the minimization will run until the system is
either right on a local minimum or *NSTEp* is reached.

There is some debate as to whether a structure should ever be minimized in
vacuum, although as mentioned above light vacuum minimization is often very
useful when preparing a structure for solvation. This is unavoidable when a gas
phase normal mode calculation is to be performed. General hints for both gas
phase and solvent minimization are given below.

If minimizing directly before a dynamics run, it is necessary to make sure that
the conditions under which the minimization takes place are the same as those
intended to be used for dynamics. For example, if particle-mesh ewald is to be
used in dynamics, it should be used for minimization as well. Consideration
should be made of the system in question; below, general hints for
minimizations in gas phase and in explicit solvent are given.

Gas phase minimization
----------------------

Gas phase minimization may be complex as, without solvent buffering, the
structure can minimize into non-physical conformations. Some find it helpful to
put harmonic restraints on the backbone atoms of the protein, letting only side
chains move. It is possible to go even farther and fix (or, at least, restrain)
all atoms except for hydrogens, allowing just their positions to be minimized.
As mentioned above, thermal B-factors can serve as guides for how strongly to
restrain individual atoms. Minimizing only hydrogens can be desirable since
CHARMM's hydrogen position builder (HBUIld) often does not put hydrogens in the
most physically advantageous positions. It is generally possible to gradually
heat a minimized structure to obtain a higher temperature configuration that is
reasonably energetically stable.

Minimizing a solvated structure
-------------------------------

As mentioned above, the conditions under which minimization takes place should
be the same under which the dynamics run will be performed. Therefore, it is
necessary to set up periodic boundary conditions and nonbonded parameters as
they will be used later. As alluded to above, it is wise to remember that the
goal of the minimization is to reduce bad contacts and relax the structure.
Using too many steps can result in over-minimization, which can cause issues
with subsequent MD simulations. As in the gas phase, it may be desirable to
constrain backbone atoms, allowing the solvent and sidechains to orient
optimally while the basic shape of the protein stays intact. Remember, the
structure can be heated to and equilibrated at the proper temperature before
starting production dynamics, but over minimization will cause these processes
to take a lot longer or possibly fail entirely. One hundred steps of SD
followed by a few hundred steps of ABNR are often sufficient to prepare a
solvated structure for dynamics.

Conclusion
----------

The general difference between vacuum and gas phase minimization is that one
must be even more careful to prevent radical changes to the structure during
gas phase. Ideally, a structure should only be minimized in the presence of
solvent, but there may be issues that need to be resolved before initiating
solvation where a short gas phase minimization would be desirable. The main
thing that must be considered when setting up the minimization of the solvated
system is making sure that, if periodic boundary conditions are to be used in
dynamics then they should be configured identically.

