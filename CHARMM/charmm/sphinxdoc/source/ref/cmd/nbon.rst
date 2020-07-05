.. index:: NBON, nbond, nbonds, non-bond, non-bonded, nonbond, nonb

.. _ref-cmd-nbon:

NBONds
======

Non-bonded interactions refer to van der Waals terms and the electrostatic terms
between all atom pairs that are not specifically excluded from nonbond
calculations. This command generates the non-bonded list and sets the
non-bonded options.

.. note:: The default values given below are commonly overridden in the
   nonbonded section of the :ref:`ref-file-prm` in your input script.

**Outline:**

  * :ref:`ref-cmd-nbon-cite`
  * :ref:`ref-cmd-nbon-syntax`
  * :ref:`ref-cmd-nbon-warnings`
  * :ref:`ref-cmd-nbon-examples`

.. _ref-cmd-nbon-cite:

References to Cite
------------------

If you use one of the following methods in a published work, please be sure to
include the relevant citations.


* PME: [Darden93]_ and [Essmann95]_
* IPS: [Wu05]_
* MPOLE: [Simmonett14]_
* ACE:
* LOOKUP: [Nilsson09]_

.. _ref-cmd-nbon-syntax:

Syntax
------

.. command:: nbonds nonbond-spec

   Takes a :spec:`nonbond` as an argument.

-------

.. todo:: Some of these fine grained control literally have *zero* real
   world usefulness to an end user. Would it make sense to undocument such
   keywords?

.. spec:: nonbond method-spec distance-spec misc-spec

   The :spec:`method`, :spec:`distance` and :spec:`misc` collectively define
   the :spec:`nonbond`. This conveniently serves an argument that can be
   passed to many CHARMM commands such as: :cmd:`update`, :cmd:`energy`,
   :cmd:`minimize`, :cmd:`dynamics`, :cmd:`vibran`, :cmd:`correl` and many
   others.

.. spec:: method

   .. keyword:: elec elec-spec
  
      Turns on electrostatic interactions and takes an :spec:`elec` argument.

   .. keyword:: noelec

      Turns off electrostatic interactions.

   .. keyword:: vdw vdw-spec
  
      Turns on vdw interactions and takes a :spec:`vdw` argument.

   .. keyword:: novdwaals
  
      Turns off vdw interactions.

   .. keyword:: ewald ewald-spec

      Computes long range electrostatics using the Ewald summation technique
      and takes an :spec:`ewald` argument. See the :ref:`conceptual guide
      <con-energy-pme>` for more information.

      Mutually exclusive with :kw:`ips`, :kw:`eips` and other long range
      electrostatic corrections.

   .. keyword:: noewald

      Turns off treatment of long ranged electrostatics.

   .. keyword:: ips ips-spec

      Computes long range electrostatics and vdw interactions using the IPS
      approximation, and takes an :spec:`ips` as an argument.

      Mutally exclusive with :kw:`ewald`, :kw:`eips`, :kw:`vips` and other long
      range non-bonded corrections.
  
   .. keyword:: eips ips-spec
  
      Computes long range electrostatics interactions using the IPS
      approximation, and takes an :spec:`ips` as an argument.

      Mutally exclusive with :kw:`ewald`, :kw:`ips`, :kw:`vips` and other long
      range electrostatic corrections.
  
   .. keyword:: vips ips-spec
  
      Computes long range vdw interactions using the IPS approximation, and
      takes an :spec:`ips` as an argument.

      Mutally exclusive with :kw:`ips`, :kw:`eips` and other long
      range vdw corrections.
  
   .. keyword:: lookup atom-selection lookup-spec

      Takes an :spec:`atom-selection` and a :spec:`lookup` as arguments. Use
      lookup tables for fast calculation of non-bonded energy and force
      calculations. See the :ref:`user guide <con-lookup>` for more details.

   .. keyword:: extended exelec-spec

      Enables extended electrostatics treatment for long ranged electrostatics.

      Requires :kw:`group` based lists.

      Mutually exclusive with :kw:`ewald` and other long ranged electrostatic
      corrections.

   .. keyword:: noextended
  
      Disables extended electrostatics.

   .. keyword:: fma fma-spec
  
      Fast multipole approximation. Takes a :spec:`fma` argument.

      .. todo:: Deprecated?

   .. keyword:: nofma

      Turns off fast mutlipole approximation.

   .. keyword:: lrc
  
   .. keyword:: lrc_ms

   .. keyword:: list
  
   .. keyword:: nolist
  
   .. keyword:: bycubes
  
   .. keyword:: bygroup
  
   .. keyword:: bycbim
  
   .. keyword:: bycc
   
.. spec:: elec

   .. keyword:: atom
      :default: True

      .. warning:: Do not combine atom based loops with :kw:`switch` or other
         switching keywords.

      Mutually exclusive with :kw:`group`.
  
   .. keyword:: group
      :default: False
      
      .. warning:: Do not combine group based loops with :kw:`shift` or other
         shifting keywords.

      Mutually exclusive with :kw:`atom`.

   .. keyword:: cdie
      :default: True
      
      Mutually exclusive with :kw:`rdie`

   .. keyword:: rdie
      :default: False
      
      Coulomb's law is now :math:`\frac{q_i q_j}{r^2}`.

      Mutually exclusive with :kw:`cdie`.

   .. keyword:: shift

      .. warning:: Do not combine group based loops with :kw:`shift` or other
         shifting keywords.

      Mutually exclusive with :kw:`switch` and other switching/shifting
      keywords.

   .. keyword:: switch
  
      .. warning:: Do not combine atom based loops with :kw:`switch` or other
         switching keywords.

      Mutually exclusive with :kw:`shift` and other switching/shifting
      keywords.

   .. keyword:: vshift

      .. warning:: Do not combine group based loops with :kw:`vshift` or other
         shifting keywords.

      Mutually exclusive with :kw:`switch` and other switching/shifting
      keywords.


   .. keyword:: vswitch
  
      .. warning:: Do not combine atom based loops with :kw:`vswitch` or other
         switching keywords.

      Mutually exclusive with :kw:`shift` and other switching/shifting
      keywords.

   .. keyword:: gshift

      .. warning:: Do not combine group based loops with :kw:`gshift` or other
         shifting keywords.

      Mutually exclusive with :kw:`switch` and other switching/shifting
      keywords.

   .. keyword:: mshift

      .. warning:: Do not combine group based loops with :kw:`mshift` or other
         shifting keywords.

      Mutually exclusive with :kw:`switch` and other switching/shifting
      keywords.

   .. keyword:: ace ace-spec
    
      Takes a :spec:`ace` argument and calculates solvation free energy and
      forces based on a continuum description of the solvent, in particular the
      analytical continuum electrostatics (ACE) potential.

   .. keyword:: ace2 ace-spec ace2-spec
  
      Takes :spec:`ace` and :spec:`ace2` arguments and invokes a modified
      treatment of the Born solvation radii which are limited by un upper bound
      -- :kw:`mxbsolv`. This takes account of the overestimation of the
      desolvation of charges by the pairwise de-screening potential in ACE1.

.. spec:: vdw

   .. keyword:: vgroup
  
   .. keyword:: vswitched
  
   .. keyword:: vatom

   .. keyword:: vshifted

   .. keyword:: vswitched

   .. keyword:: vfswitch

   .. keyword:: vgshift

      Use GROMACS style shifting for vdw interactions.

   .. keyword:: vtrunc

      Only works when compiled with ##MMFF flag.
  
      .. todo:: Can this KW be removed?

   .. keyword:: ctvt real

      Only works when compiled with ##MMFF flag.
  
      .. todo:: Can this KW be removed?

.. spec:: ips

   Isotropic periodic sum. See the :ref:`Conceptual guide <con-ips>` for more information.

   .. keyword:: raips real
      :length: 5
  
   .. keyword:: rips real

   .. keyword:: nipsfrq real
  
   .. keyword:: dvbips real
   
   .. keyword:: mipsx int
      :length: 5

   .. keyword:: mipsy int
      :length: 5

   .. keyword:: mipsz int
      :length: 5

   .. keyword:: mipso int
      :length: 5

   .. keyword:: PXYZ

   .. keyword:: PYZ

   .. keyword:: PZX

   .. keyword:: PYX

   .. keyword:: PZY

   .. keyword:: PXZ
      
   .. keyword:: PX

   .. keyword:: PY

   .. keyword:: PZ

.. spec:: ewald

.. spec:: mpole

.. spec:: distance

   Cutoff values for building non-bond lists determining switching function
   regions and ignoring non-bonded terms are all given in Angstrom.

   See the :ref:`conceptual guide <con-energyfunctions-cutoffs>` for more information.

   
   .. keyword:: cutnb real
      :default: 8

      Cutoff value for building the list of non-bonded interactions.

   .. keyword:: ctonnb real
      :default: 6.5

      Defines the start of the switching region for CHARMM style switching function.

   .. keyword:: ctofnb real
      :default: 7.5

      Cutoff for CHARMM switching function.

   .. keyword:: cgonnb real
      :default: 0.0

      Defines the start of the switching region for GROMACS style switching function.

   .. keyword:: cgofnb real
      :default: 10.0

      Cutoff for GROMACS style switching function.

   .. keyword:: wmin real
      :default: todo
      
      ???
      
   .. keyword:: wrnmxd real
      :default: todo
      
      ???

.. spec:: lookup

   Keywords for fast lookup tables.

   .. keyword:: interpolate
      :default: True
      
      Use linear interpolation for lookups.
      
   .. keyword:: nointerpolate
  
   .. keyword:: tabincr int
      :default: 20
      
      Determines the size of the lookup table. Size = :kw:`tabincr`
      :math:`\times` :kw:`ctofnb`:sup:`2`.

   .. keyword:: noenergy
      :default: True
      
      Energies will only be evaluated when non-bond list is updated.
      
   .. keyword:: energy
      :default: False
      
      Energies will always be evaluated.

   .. keyword:: novu
      :default: False
      
      Do not use lookup for the sol **V** ent-sol **U** te interactions.

   .. keyword:: nouu
      :default: False
      
      Do not use lookup for the sol **U** te-sol **U** te interactions.
  
   .. note:: 1-4 interactions are always handled by standard routine

.. spec:: rxnfld

   Reaction field.
  
   .. keyword:: norxn
      :default: True

   .. keyword:: rxnfld
      :default: False
      
   .. keyword:: rxnnb
      :default: False
      
   .. keyword:: epsext real
      :default: todo

   .. keyword:: orderrexnfld int
      :default: wtf

   .. keyword:: shell real
      :default: todo

.. spec:: ace

   .. keyword:: ieps real
      :default: 1.0

      Dielectric constant for the space occupied by the atoms that are treated
      explicitly, *e.g.*, the space occupied by the protein. 

   .. keyword:: seps real
      :default: 80.0

      Dielectric constant for the space occupied by the solvent that is treated
      as a continuum (*i.e.*, the complement of the space occupied by the protein). 

   .. keyword:: alpha real
      :default: 1.3
      
      The volumes occupied by individual (protein) atoms are described by
      Gaussian density distributions. The factor ALPHa controls the width of
      these Gaussians.

   .. keyword:: sigma real
      :default: 0.0

      The ACE solvation potential includes a hydrophobic contribution which is
      roughly proportional to the solvent accessible surface area.  The factor
      SIGMa scales the hydrophobic contribution. 

   .. keyword:: ideal
      :default: True
      
      ACE calculates the nonbonded exclusion list distances from ideal bond
      length and angles where possible; the distances for 1-4 atom pairs in the
      exclusion list are calculated from the current atom positions at the
      first ACE energy call.

      See :kw:`current`.

      Mutually exclusive with :kw:`current`.
      
   .. keyword:: current
      :default: False

      All the distances between atoms in the nonbonded exlusion list are
      calculated from the current coordinates of the atoms.

      See :kw:`ideal`.

      Mutually exclusive with :kw:`ideal`.

   .. keyword:: fvscale real
      :default: 1.0

      Volume scale factor, to reduce systematic errors of ACE1 model.

.. spec:: ace2

   .. keyword:: mxbsolv real
      :default: 14.0

      The Born solvation radii of all atoms (charges) are limited by the upper
      bound parameter.

   .. keyword:: tbsolv real
      :default: 8.4

      In the ACE2 potential, the conventional conversion of the atomic
      solvation to the Born solvation radii is applied until a Born radius of
      TBSOlv is obtained ("turning point") 

   .. keyword:: tbshy real
      :default: 3.85

      This parameter has the same meaning as :kw:`tbsolv`, but applies only to
      hydrogens.

.. spec:: fma
  
   .. todo:: Deprecated?
  
   .. keyword:: level int
  
   .. keyword:: terms int

.. spec:: exelec

   .. todo:: Deprecated?

   .. keyword:: gradient
  
   .. keyword:: nogradient
  
   .. keyword:: quadrupole
  
   .. keyword:: noquadrupole

   .. keyword:: ctexnb real
      :default: 999

      Defines the cutoff distance beyond which interaction pairs are excluded
      from the Extended Electrostatics calculation.
      
.. spec:: misc

   Garbage collector for other non-bond options.

   .. keyword:: inbfrq int
      :default: -1
     
      The frequency in steps for recomputing the non-bonded list. The default
      value of `-1` uses a heuristic updating frequency determined by :kw:`cutnb`
      and :kw:`ctofnb`.

   .. keyword:: eps real
      :default: 1.0
      
   .. keyword:: e14factor real
      :default: 1.0
      
   .. keyword:: nbxm int
      :default: todo
      
   .. keyword:: nbscale real
      :default: todo
      
   .. keyword:: imscale real
      :default: todo
  
   .. keyword:: exoforce
  
   .. keyword:: init

      .. todo:: Deprecated?

   .. keyword:: reset

      .. todo:: Deprecated?


.. _ref-cmd-nbon-warnings:

Warnings
--------

The following combinations of keywords known to produce bad results.

=========== ========== ============ === ================================
:kw:`atom`  :kw:`cdie` :kw:`shift`  no  (obsolete, but used in the past)
:kw:`atom`  :kw:`cdie` :kw:`switch` no  Very bad - do not use
:kw:`group` :kw:`cdie` :kw:`shift`  no  (obsolete)
:kw:`group` :kw:`cdie` :kw:`switch` no  Very bad with non-neutral groups
:kw:`atom`  :kw:`rdie` :kw:`shift`  yes but do you really want RDIE?
:kw:`atom`  :kw:`rdie` :kw:`switch` no  switch is bad here.
=========== ========== ============ === ================================


.. _ref-cmd-nbon-examples:

Examples
--------

GROMACS style shifting functions:

.. code-block:: chm

   NBONds CUTNb 14 CGONnb 0.0 CGOFnb 12.0 ATOM CDIE GSHIft -
       CTONnb 9.0 CTOFnb 12.0 VATOM VGSHift
