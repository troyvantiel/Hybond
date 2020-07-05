.. _usr-script-io:

File I/O
========

CHARMM provides a consistent set of commands for reading and writing files. In
general, a file must be *OPENed*, read from or written to, and then, in
some cases, *CLOSEd*. To provide a basic example, here is how to read an
already-created protein structure file (PSF) into CHARMM:

.. code-block:: chm

 open read card unit 10 name myfile.psf
 read psf card unit 10
 close unit 10

Already, from this simple example, there are several features that must be
noted. First and foremost, it is necessary to associate each open file with a
unit number (FORTRAN programmers will recognize this immediately). Unit numbers
must be in the range from 1-99 for CHARMM versions c35 and older. It is
generally good practice to use numbers greater than 10, as the lower number
units may be in use by the system (e.g. the standard output stream is usually
written to unit 6 on modern systems). Another interesting point is the use of
the *CARD* subcommand. CHARMM can read and write two types of files, those
containing ASCII text are designated as *CARD* or *FORMatted*, and those
containing binary data are designated *FILE* or *UNFOrmatted*. By specifying
*CARD* in the *OPEN* and *READ* commands, we are telling CHARMM that we expect
the data to be ASCII text (as, indeed, PSF files are). Also in the *READ*
command we must tell CHARMM which unit we want to read from, since we can have
a number of files open at a given time.

Once we are done with the file it is good practice to *CLOSe* it. If you
attempt to re-open the unit number before closing the file, CHARMM will close
the file for you, but it's cleaner and less confusing to explicitly close it
yourself. Note that when writing files, they are generally closed automatically
after writing is complete.

In this particular simple example, we can compress the syntax to

.. code-block:: chm

 read psf card name myfile.psf

In this case, CHARMM will assign a unit number to myfile.psf, open it, read it,
and close it all with one command.

Writing works similarly to reading, the file is opened written to, and closed.
For example, one may wish to write out the current coordinate set. This can be
done as follows:

.. code-block:: chm

 open unit 11 write card name new-coordinates.crd
 write coor card unit 11
 * title of the coordinate file
 *

Note that in this case, we did not close the file as CHARMM closes it
automatically once the writing is done.

To give an example of opening an unformatted/binary file, consider the case of
opening a file to write a coordinate trajectory for dynamics (trajectories are
discussed further in the dynamics and analysis sections):

.. code-block:: chm

 open unit 20 write unfo name my-trajectory.trj
 dyna ... iuncrd 20 ...

In this case no title is written (after all the file contains binary data) and
the unit number is passed to the *IUNCRD* option to the
*DYNAmics* command, which specifies the unit number on which to write
out the coordinate trajectory.

A complete treatment of how CHARMM handles file I/O can be found in
`io.doc <http://www.charmm.org/documentation/current/io.html>`_. Full-fledged
examples are given on the [[FINAL Full example|complete example page]], and an
example of opening trajectory files is given on the [[FINAL Analysis|analysis
page]].

Main and Comparison Coordinate Sets
-----------------------------------

During a run, CHARMM stores two sets of coordinates; the coordinates being
acted upon are called the main set (*MAIN*), but there is another set
available called the comparison set (*COMP*). This is extremely useful
for measuring how much atoms have moved around during a particular action, for
example...

.. code-block:: chm

 ! copy the contents of the main set into the comp set
 coor copy comp
 
 ! do something to move the coordinates such as
 ! minimization
 .
 .
 .
 
 ! get the difference between the main and the
 ! comparison set and save it in the
 ! comparison set
 coor diff comp
 
 ! print out the differences
 print coor comp

You can also use

.. code-block:: chm

 coor copy

to copy the *COMP* set into the *MAIN* set and

.. code-block:: chm

 coor swap

to swap the main and comparison sets.

One important caveat is that some calculations such as molecular dynamics may
overwrite the *COMP* coordinate set. If in doubt, you should save the
coordinates to file and then read them back into the *COMP* set, *e.g.*:

.. code-block:: chm

 write coor card name temp.crd ! write set out to a file

 ! do something that might overwrite the comp set, e.g. dynamics
 .
 .
 .

 ! re-fill the COMP set from the file we wrote before
 read coor comp card name temp.crd

This way, you will be sure you have the coordinates that you're expecting in *COMP*!

.. _usr-script-io-select:

Atom Selection
--------------

One of the most useful features of CHARMM is its powerful atom selection
functionality. This can be used to modify the effects of certain commands,
making them apply to a subset of the atoms or the whole system. For example, if
you want to write only some of the coordinates out, there's a command for that.

All atom selections are done with the *SELEct* subcommand.
*SELEct* is not a command itself, but is often used as an integral part
of other commands to determine which atom(s) will be operated upon (there will
be plenty of examples of this as we move through the tutorial). The basic
syntax for *SELEct* is:

.. code-block:: chm

 sele <criteria> end

Where <criteria> is the specification of which atoms you want selected (this
tutorial will give plenty of examples of which criteria you can use). When
experimenting with the SELEct command, it is helpful to put the subcommand SHOW
directly before the END statement. This will tell CHARMM to list which atoms
were actually selected, which is very helpful for debugging SELEct statements!

A full description of atom selection is contained in
`select.doc <http://www.charmm.org/documentation/current/select.html>`_.

Controlling which atoms are selected: BYREsidue and BYGRoup
***********************************************************

Usually CHARMM will only mark those atoms which match the selection criteria.
However, sometimes you want to select all atoms in a group or residue in which
one atom matches the criteria. The *.byres.* and *.bygroup.* key
words (note the leading and trailing dots) allow you to do this. For example:

.. code-block:: chm

  sele .byres. bynum 12 end

will select all atoms in the same residue (*.byres.*) as atom number 12 (the
*bynum* factor is described later on, but its use should be self-evident here).

Basic operators: .AND., .OR., and .NOT.
***************************************

It is possible to use basic `boolean logical operators
<http://en.wikipedia.org/Boolean_logic>`_ (.AND., .OR., and .NOT.: the periods
are optional) with atom selections, and they behave exactly as one would
expect. For example, if you need to select atom numbers 12 and 15, this can be
done by:

.. code-block:: chm

  sele bynum 12 .or. bynum 15 end

A commonly seen mistake is to do:

.. code-block:: chm

  sele bynum 12 .and. bynum 15 end

This selection will return no results because it is impossible for an atom to
be both number 12 and number 15 at the same time!

Likewise, you can select all atoms except for atom number 12 by doing:

.. code-block:: chm

  sele .not. bynum 12 end

Basic factors: atom, segment, and residue types
***********************************************

There are a number of keywords that let you select an atom based on a
particular characteristic. We've already seen one of these, the *bynum*
command which select an atom based on (surprise, surprise) its number. An
advanced feature of the *bynum* token is that you can select a range of
numbers by seperating them with a colon, for example:

.. code-block:: chm

 sele bynum 12 : 20 end

selects atoms 12-20. One of the most important factors for day to day use is
the ATOM token. The syntax is:

.. code-block:: chm

 sele atom <segment ID> <residue ID> <type> end

Any of the segment ID, residue ID, or type can use wildcards. A "*" matches any
string of characters, a "#" matches any string of numbers, and "%" and "+"
match a single character or number, respectively. So for example:

.. code-block:: chm

 sele atom A * C* end

will select all carbon atoms in all residues of segment A (note that C* matches
C, CA, CT, etc.).

It is possible to select by the segment ID, residue ID and type individually
with the *SEGId*, *RESId*, and *TYPE* factors, so for
example:

.. code-block:: chm

 sele segid A .and. type C* end

will provide equivalent results to the previous atom selection.

Another example of what you can do is select a range of residues. Suppose, for
example, residues 22 through 48 of a protein make up an alpha helix and you
want to select all of them to perform some action, you can do:

.. code-block:: chm

 sele resid 22:48 end

The final basic factor that is commonly used is the *RESName* criteria, which
selects atoms based on their residue name. For example:

.. code-block:: chm

  sele resn GLY end

will select all atoms in all Glycine residues.

Advanced operators: AROUnd and BONDed
*************************************

It is also possible to select atoms based on their spatial or bonded
relationship with other atoms. To select atoms within a given distance of a
group of atoms, one can use the .AROUND. followed by a real number immediately
after another factor. For example:

.. code-block:: chm

 sele resid 10 .around. 5.0 end

will select all atoms within 5 angstroms of residue number 10.

Likewise, you can select the atoms that are bonded to a particular atom with the *.BONDED.* modifier:

.. code-block:: chm

 sele atom * * H* .and. .bonded. type C* end

will select all hydrogens bonded to a carbon atom.

Additional factors
******************

There are several other keywords that you can use to select atoms

* INITial: This keyword selects all atoms that do not have coordinates initialized (CHARMM sets uninitialized coordinates to 9999.9999).
* HYDRogen: selects all of the hydrogen atoms in the system
* CARBon: selects all of the carbon atoms in the system

Properties
**********

The PROPerty key word allows for atom selection based on CHARMM's SCALar
properties, a full list of which may be found in `scalar.doc
<http://www.charmm.org/documentation/current/scalar.html>`_. Both the main and
comparison coordinate and weighting arrays are permitted. For example:

.. code-block:: chm

 sele property wcomp .eq. 1 end

would select those atoms having a weight of 1 in the comparison weighting
array. Other scalar values that may be selected include the X, Y, and Z
coordinates (called XCOMP, YCOMP, and ZCOMP for the comparison coordinate set),
the mass of the atom (MASS), and several others. The complete list is in
`select.doc <http://www.charmm.org/documentation/current/select.html>`_.

One interesting note is that you can use this to print atoms that move more
than a particular amount during a simulation, *e.g.*

.. code-block:: chm

 ! copy current coordinates to the comparison set
 coor copy comp

 ! minimize, or do dynamics or whatever

 ! compute the differences between the new coordinates and the ones saved previously, 
 ! storing these differences in the comp set
 coor diff comp 

 ! print out those atoms that are displaced more than 10 angstroms
 ! in any directions
 define bigmove sele prop xcomp .gt. 10 .or. prop ycomp .gt. 10 -
  .or. prop zcomp .gt. 10 show end

There is a new command in the example above -- DEFIne. DEFIne allows you to
associate a name for an atom selection and then use it later, for example:

.. code-block:: chm

 define sel1 <big nasty long atom selection>
 .
 .
 .
 print coor select sel1 end

In the above example the &lt;big nasty long atom selection&gt; can be used
multiple times without having to retype it. To give a concrete example, if you
are going to do operations repeatedly on all atoms withing 5 angstroms of
residues 10 through 20, you can do

.. code-block:: chm

 define critreg select resid 10:20 .around. 5 end

 ... ! do some stuff

 coor stat select critreg end ! get coordinate statistics
                              ! of atoms in critreg only

Note that after defining *critreg* it is necessary to encapsulate it within a
*SELEct* statement, *i.e.* *SELEct critreg END*

Variables in CHARMM
-------------------

User-settable Variables
***********************

So far, the CHARMM scripting language seems to be a concatenation of individual
commands. It does contain, however, most (if not all) elements of a
(imperative) programming language (even though it is not an extremely
comfortable one). One key element of a programming language is the ability to
set and manipulate variables. In the context of running CHARMM, it is extremely
useful to pass certain variables into the program from the outside. Run the
following miniscript (we name it title2.inp)  

.. code-block:: chm

 * Example of passing a variable from the command line  
 * The variable myvalue was initialized to @myvalue  
 *  
 stop 

by executing ::

 charmm_executable myvalue=500 < title2.inp 

and study the output (you could also state *myvalue:500*): Before the
title is echoed, you see that an argument was passed  ::

 Processing passed argument "myvalue:500"  
 Parameter: MYVALUE <- "500" 

and in the echoing of the title, @myvalue is replaced by the value assigned to
the variable "myvalue", *i.e.*, 500 in our case. The little example highlights
something to keep in mind. The first thing that is done when a line in a script
(containing a title or a command) is parsed, is to scan for @variables, which
are then replaced by their value (no variable replacement is done in
comments!).  

Substitution Variables
**********************

For completeness sake, a preliminary comment on a second type of variable is
needed. Many CHARMM commands, aside from producing more or less direct output,
place some key values in "internal" variables. *e.g.*, the energy of a system
calculated with the most recent *ENERgy* command is put into a variable named
*ENER*. To avoid name clashes with variables set by users, the values of these
variables can be accessed by preceding the variable name by a question mark
"?". For example, if CHARMM encounters (in a title or in a command) the string
*?ENER*, it will attept to replace it with the energy value from the most
recent *ENERgy* command. This is quite handy, and we'll see many examples later
on, once we have reached the stage where we can use CHARMM to work with
biomolecular systems.

Loops and flow control
----------------------

The beginning of loops in CHARMM is marked by a *label*

.. code-block:: chm

 label someplace

basically anywhere in a CHARMM script. If (before or after) that *label* CHARMM
encounters a 

.. code-block:: chm

 goto someplace

command, the script continues from "*label someplace*"; *i.e.*, control is
transferred to that line of the script.  The *label/goto* tandem, combined with
*ifs*, make it possible to construct almost arbitrary control structures, in
particular loops (see below). First, however, two simpler examples.

**Example 1:** We just showed how to provide a default value to a variable that
is expected to be passed from the command line. Obviously, this is not possible
or sensible in all cases. The script *simplemini.inp* expects two parameters,
*@nstep* and *@system*. While it makes sense to set a default for the former
(if the user forgets to specify a value), the second variable has to be set to
a sensible value. Thus, we may want to check whether *@system* is initialized,
and if not, write a meaningful error message and exit the script. This is done
by querying *@?system* -- the *@?variable* operator is 1 if *@variable* is set
and 0 if it is not. The following script fragment shows how to use this

.. code-block:: chm

 * The title
 * ...
 *
 if @?nstep eq 0 set nstep 100 ! provide default for nstep
 if @?system eq 0 goto systemerror
 
 ... ! continue with normal script up to normal
 
 stop
 
 label systemerror
 
 echo You need to pass a value for variable system
 echo
 echo Aborting execution
 
 stop

**Example 2:** It was pointed out that many CHARMM scripts do identical things
(read rtf, params, psf, coords, some more stuff) whereas the "genuine" work
part (minimization, molecular dynamics, etc.) consists of just a few lines.
Thus, the first twenty to forty lines of many CHARMM scripts contain
essentially the same commands. Using *label/goto* statements one can reorganize
such scripts, so that the boring stuff is placed after a label at the end of
the scripts. Thus, the more unique parts of the script can be seen earlier in
the file (note that **stream files provide a better way of accomplishing the
same goal**, but the technique may prove useful in other cases). Take a look at
`simple_mini2.inp
<http://www.mdy.univie.ac.at/charmmtutorial/simple_mini2.inp>`_.  After
checking whether *@nstep* and *@system* were passed from the command line,
command is transferred to *label setupsystem*. The corresponding code is at the
end of the script; from there, control is transferred back to the beginning of
the script (*label fromsetupsystem*). The first interesting line (*mini sd*)
moves up about 15 lines; one can understand a bit more quickly what the script
does.

**Example 3:** As mentioned above, *label* statements provide for a simple way
to make a loop. The following example shows a loop that repeats an action ten
times over

.. code-block:: chm

 set i = 0

 label beginloop
 incr i by 1
 ! ... this is the body of the loop ...

 if @i .lt. 10 then goto beginloop

 ! commands below this point are executed
 ! after the loop finishes

Pay attention to how the *INCRement* command is used to add 1 to the value of
*@i* for each iteration of the loop.

Constraints and Restraints
--------------------------

Basic overview
**************

CHARMM has the ability to constrain or restrain atoms using the CONStraint
command. This can be used to ensure that atoms stay close to a given point or
another atom during a simulation. As an example, it is often desirable to
constrain or restrain protein backbone atoms during minimization and only
optimize the positions of the various side chains while the basic trace of the
backbone remains fixed. Constraints and restraints can also be used to reduce
high-frequency motions in the system as in the case of SHAKE, which constrains
bond lengths.

The basic difference between a constraint and a restraint is that a constraint
completely removes the given degree(s) of freedom from the system while a
restraint retains them but applies a penalty function to the energy if the
atom(s) involved move away from the desired configuration. More information can
be found in the `cons.doc
<http://www.charmm.org/documentation/current/cons.html>`_ CHARMM documentation
file. Below we describe three of the most commonly seen restraints.

Harmonic restraints
*******************

Harmonic restraints are used to hold atoms near a given location by applying a
harmonic penalty function as they move away from the desired location. There
are three types of harmonic restraints: absolute, relative, and bestfit.

* Absolute restraints require the atom(s) to stay near a given cartesian
  position. A set of reference coordinates must be given. By default, the
  current position in the main coordinate set is used, but the comparison set
  may also be used. It is not possible to specify the reference coordinates
  directly in the command; they must be in either the main or comparison
  coordinate set.

* Bestfit restraints are similar to absolute restraints, but rotation and
  translation are allowed to minimize the restraints energy.This makes it
  useful for holding inter-atom geometry intact while allowing for block
  motion. Second derivatives are not supported for bestfit restraints.

* Relative positional restraints are similar to bestfit restraints, but there
  are no reference coordinates. It is used to force two sets of atoms to have
  the same shape

One important caveat is that no atom may participate in more than one restraint
set. It is possible in the case of *ABSOLUTE* restraints to specify the
exponent for the harmonic restraint (*i.e.* where the penalty function is linear,
quadratic, cubic, quartic, *etc.*) and to scale the restraint separately in the
x, y, and z directions.

Some examples are:

.. code-block:: chm

 define backbone sele type N .or. type CA .or. type C end
 cons harm absolute sele backbone end

restrains the protein backbone based on its current coordinates. This can be
handy in minimization where you don't want to change the basic shape of the
atom from that in the crystal structure.

.. code-block:: chm

 cons harm absolute expo 4 xscale 1.0 yscale 0.5 zscale 0.0 sele backbone end

restrains the backbone using a quartic penalty function in the x direction. The
restraint energy in the y direction is halved, and the restraint is not applied
in the z direction. Such an unusual restraint would probably not be used much
in normal circumstances. 

.. code-block:: chm

 cons harm bestfit sele resid 20 end

restrains residue #20 to its current geometry. For example, if you have a
ligand whose internal conformation should remain rigid but that is allowed to
rotate and translate, a restraint like this could be used.

.. code-block:: chm

 cons harm relative sele segid 1 end sele segid 2 end

restrains segment 1 to maintain the same geometry as segment 2.

Note that the *MASS* keyword may be used to mass weight the restraints.

All harmonic constraints can be cleared with the command

.. code-block:: chm

 cons harm clear

Distance Restraints
*******************

Distance restraints require two atoms to remain within a given proximity of
each other. The desired distance, force constant, and exponent must be given to
the RESDistance command. For example:

.. code-block:: chm

 resd kval 10.0 rval 5.0 ival 2 sele atom 1 1 CA end sele atom 1 2 CA end

restrains the alpha carbons of residues 1 and 2 of segment 1 to be 5 angstoms
apart, with a quadratic penalty function. The functional form of this energy
term is:

 :math:`E_{resd} = Kval(r-rval)^{ival}`

This energy term is very similar to the bonded energy term!

For individual pairs of atoms, NOE restraints can also be used (this is very
good for distances which depend on one another). Discussion of this
functionality is beyond the scope of this tutorial, but you can find more
information in the NOE sections of `const.doc
<http://www.charmm.org/documentation/current/cons.html>`_

SHAKE
*****

SHAKe is a constraint that fixes bond lengths. It can be used to fix angles as
well, but this is not recommended. It is most widely used to fix the lengths of
bonds involving hydrogens, since these tend to vibrate at very high
frequencies, which can lead to numerical imprecision during simulations. The
use of SHAKe on bonds involving hydrogens allows for numerical integration at a
1-2 fs step size during dynamics.

To apply shake to all bonds involving hydrogens do

.. code-block:: chm

 shake bonh param sele all end

Note that the *sele all end* is the default atom selection when none is given
(*i.e.* it's a bit redundant here). The *PARAm* keyword uses the values in the
parameter file as the optimal bond lengths. If it is not used, the current
distance from the main coordinate set is used (or the current distance from the
comparison set if the *COMP* keyword is given in place of *PARAm*).

Debugging CHARMM scripts
------------------------

Debugging via prnlev
********************

The *PRNLev* command controls how verbose CHARMM is in its output. The default
value is 5, the minimum level is 0 (which will cause CHARMM to be virtually
silent), and the maximum level is 9. Higher values of *PRNLev* produce more
output. This command also controls which node prints output in a parallel job.
The default is that all nodes produce output, which is usually not desirable.
The *NODE* subcommand to *PRNLev* restricts print out to a single node, *e.g.*:

.. code-block:: chm

 prnlev 6 node 0

Sets the *PRNLev* to 6 and states that only the first node (index 0) should produce output.

At higher values of PRNLev, CHARMM will detail exactly which energy terms are
being computed and give details about which subroutines are being called. This
can be very helpful for debugging incorrect results.

The *BOMBlev* and *WARNlev* commands are similar, except that they control the
severity of errors which will cause CHARMM to abort execution or issue a
warning. Errors and warnings range in severity from +5 (least severe) to -5
(most severe). One common mistake made by novices is to set *BOMBlev* too low.
Level 0 errors and below are serious problems that may affect the validity of
numerical results from the runs. Level 0 and -1 errors should only be accepted
when the user understands what they mean and has determined that they will not
affect the results. Errors at the -2 level or below are very severe and
generally mean that the results of the run are invalid. In no cases should the
BOMBlev ever be set below -2. In general, it is always safe and recommended to
set *BOMBlev* to 0 (although it may need to be temporarily lowered for some
operations).

Printing data structures
************************

CHARMM provides the functionality to print out data structures. For example,
the commands:

.. code-block:: chm

 print psf
 print coor

would print out the current protein structure file and coordinate set.
Parameters can be printed in the same way, however in many cases most of the
parameters that are read in do not apply to the structure being studied. The
command

.. code-block:: chm

 print param used

only prints out those parameters that are currently referenced by the PSF.

It is also possible to display  the current forces acting on each atom, for
example:

.. code-block:: chm

 coor force comp
 print coor comp

Puts the forces into the COMParison coordinate set and then prints them out. It
is possible to combine these with a *PROPerty* based atom selection (described
above) to print all forces over a certain magnitude, *e.g.*:

.. code-block:: chm

 coor force comp
 print coor comp sele prop abs x .gt. 500.0 .or. prop abs y .gt. 500.0 .or. prop abs z .gt. 500.0 end

will print all atoms who have at least one of their force components greater
than 500 in absolute value.

If a more detailed analysis is needed (for example, if the force on one of the
atoms is blowing up for an unknown reason), it is possible to print out all of
the individual terms of the first field, this is done by the *ANAL TERM*
command. By default, *ANAL TERM* only prints out the bonded interaction terms,
*ANAL TERM NONBonded* will print all of the terms. It is possible to see in
which order the energy functions are called by setting the *PRNLev* to 9.

