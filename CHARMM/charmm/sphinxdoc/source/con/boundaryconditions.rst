.. index:: Periodic Boundary Conditions, Minimum Image Convention

.. _con-pbd:

Periodic Boundary Conditions
============================

The textbook view of PBCs
-------------------------

We trust that you have read some background material on molecular mechanics and
molecular dynamics simulations, and, thus, expect you to have some familiarity
with periodic boundary conditions (PBC). The following is intended as a
refresher with emphasis on pointing out where things may not be so obvious for
macromolecular systems. 

Even with today's computers, the system sizes we can handle are microscopically
small. Thus, without any special precautions we would be simulating
"infinitesimal droplets", the properties of which would be dominated by surface
effects. It is standard practice to work around this by employing periodic
boundary conditions (PBC), *i.e.*, by assuming that the central simulation system
is surrounded in all spatial directions by exact copies (this assumes, of
course, that you have a space-filling geometry; obviously, you can't have PBC
for a spherical system. Anticipating CHARMM terminology, we refer to these
surrounding boxes as images (image boxes). Unless otherwise stated, we'll
assume a cubic simulation box (note, though, that CHARMM can handle any valid
crystallographic space group for PBC!)

The role of the surrounding boxes (images) is to replace particles that leave
the central simulation box. Consider a minimization or MD simulation, when the
positions were just updated. You know the basic picture: as a particle leaves
the central (cubic) box, say, to the right, it is replaced by the corresponding
particle from the image of the central box on the left.  Thus, working with
PBCs entails making sure that your particle coordinates always fall within the
central simulation box. This can be accomplished with the following pseudo-code
(which only handles periodic boundaries in the x direction) computer code along
the lines (the center of the box is supposed to be at the origin, x is the x
coordinate of the particle of interest, and xsize is the size of the periodic
box!) ::

 if (periodicx) then
   if (x <  -xsize/2.0) x=x+xsize
   if (x >=  xsize/2.0) x=x-xsize
 endif

So, for example, if we have a particle with a x position of 4 and xsize is 10,
the particle would be within the periodic box, however, if the particle's x
coordinate was 6, it would be outside of the box and it's x coordinate would be
changed to -4 (effectively, when a particle leaves a box on one edge, it
reappears on the opposite edge). Analogous code snippets apply to the y and z
coordinate, and for a cube *xsize = ysize = zsize=L*

There is, however, a second, more important component to PBC. Consider the
interaction (LJ and/or electrostatic) between two particles i and j. If the
replicas of the central box were real, we would have interactions between i and
j in the central box, but also of i with j in all of the (infinite) number of
images of the central box. This is, of course, not what we want to have;
instead, we choose the interaction between i and j having the shortest
distance. This may be the interaction between i and j in the central box, but
it may also be the interaction between i and a particular image atom of j. In
pseudo-code this so-called *minimum image convention* can be expressed as ::

 if (periodicx) then
   dx = x(j) - x(i)
   if (dx >   xsize/2.0) dx = dx - xsize
   if (dx <= -xsize/2.0) dx = dx + xsize
 endif

where x(i) and x(j) are the x coordinates of i and j, respectively, and
analogous statements for y and z. If you consider the above realization of the
minimum image convention carefully, you'll note that the longest possible
inter-particle distance is half the box-length, L/2. Thus, if you used a
cut-off radius greater than L/2, your effective cut-off radius would still be
L/2. As we shall see, CHARMM handles PBC rather differently, and for the
minimum image convention to be obeyed the cutoff radius must not be larger than
L/2, *i.e.*, it is the user's responsibility to shorten CTOFnb if necessary!

Again, we assume that you have heard all this before. The problem with the
above illustrative examples of PBC/minimum image convention is that the code
assumes a collection of atoms. Consider, *e.g.*, a box of waters or a box
containing a protein surrounded by water. Suppose one water edges out of the
box, *i.e.*, the oxygen and one hydrogen are still within the central box, but
the second hydrogen is out of the box. Then naively applying periodic boundary
conditions will tear the molecule apart, *i.e.*, one has to make sure that a
water is shifted periodically only as a complete molecule. Similar
considerations apply of course to the protein, or any molecule present in the
system.  Also, the PBC/minimum image checks have to be carried out for each
energy calculation. For cubic boxes, this is not too bad, and by optimizing the
above code snippets, the overhead is small. However, for more complex periodic
boxes (truncated octahedron, rhombic dodecahedron), this may quickly become
difficult to calculate, however CHARMM will handle this for you as long as you
choose one of the periodic boxes that it supports (although you can define your
own space group, this is not necessary for setting up basic dynamics
simulations).

.. todo:: This page is to be completed by :ref:`developers-rmv`.

