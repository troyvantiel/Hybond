About
=====

This document describes the Sphinx documentation system in CHARMM.  There are a
number of standards defined below, so please read this before adding new
documentation.



Required software
=================

To get the basic HTML version of the manual to work, you only need Sphinx
itself (and Python).  If you have python installed, you should be able to run

    [sudo] easy_install sphinx

to get this.  To test if you have sphinx installed, try running

    sphinx-build

from the command line.  To be able to generate PDFs, you'll need pdfLaTeX
install, too.



The manual syntax
=================

The manual is written in a meta-language called reStructuredText (ReST), which
is quite common and there are many tutorials and cheat sheets available, e.g.:

    http://www.siafoo.net/help/reST
    http://docutils.sourceforge.net/docs/user/rst/quickref.html

The reason for choosing this language is that it gives us a very easy way to
generate many different formats from a single source, using the Sphinx engine.
Everything can be hyperlinked, which allows us to connect the many different
sections of the manual in a very user-friendly way.  Although the language is
generally insensitive to spacing, careful attention must be paid to get the
indentation correct; this is how paragraphs and section are separated and
parsed.  For example, the following usage of the command block:

/------------------------------------------------------------------------------\
|                                                                              |
|.. command:: some_charmm_command                                              |
|                                                                              |
|   This sentence is aligned with the "c" in command, so it is an argument to  |
|   the "command" function.                                                    |
|   This sentence has the same alignment, so it is also an argument.           |
|                                                                              |
|This sentence is left-aligned, and is a regular paragraph.                    |
|                                                                              |
\------------------------------------------------------------------------------/

Most of the functions in ReST expect a newline between the command and the
arguments to it, as in the example above.




The layout of the manual
========================

The manual is categorized into sections, per the discussion at the CHARMM
meeting (NIH, Sept. 2014).  The directory structure is designed to mirror the
layout of the document(s).  At the top level there are four content folders:-

 * con - The "concepts" section, describing the theory behind methods.  This
         should not reference CHARMM-specific advice at all, but give an
         insightful overview of the underlying theories that we support,
         including references.
 * ref - The "reference" section, which describes command syntax (like the
         existing .doc files do).  This should be a very brief, technical
         description of the syntax expected in CHARMM input files, and each
         keyword should be described using the markers described below.  Take a
         look at the source/ref/cmd/nbon.rst file for a demonstration.
 * usr - The "user manual" section, which describes some recommended practices
         and examples for users.  This should be tutorial-like and should link
         in any germane concepts or reference sections, while providing the
         reader with a set of best practices and recommendations, as well as a
         clear description of any potential pitfalls.  This section can be used
         to demonstrate the advanced features of CHARMM.

All literature citations should go in bib.rst, and should be linked using the
syntax [AuthorYY]_ in the main text.



Namespaces and links
====================

Everything in ReST can be linked.  To make this happen, we should tag our
sections, subsections, equations, images, etc.  For example, the following text
may appear in the "user manual" section of the document.

/--------------------------------\
|                                |
|.. _usr-runningnptdynamics:     |
|                                |
|Running NPT Dynamics            |
|====================            |
|                                |
|Blah blah blah blah, stuff..... |
|                                |
\--------------------------------/

The tag ".. _usr-runningnptdynamics:" allows us to link to this section of text
from absolutely anywhere in the text, by simply inserting the command

  :ref:`usr-runningnptdynamics`

anywhere else in the document.  In order to make the entire document
consistent, we therefore enforce the following conventions.

 * The label comprises the following elements (in order, as demonstrated above):-

    * namespace "usr, ref or tut" then a hyphen (*NOT* an underscore).
    * the name of the subdirectory (if any) that the current file resides in.
    * the title that the tag labels, but all lowercase and with no spaces.

 * The tag must immediately precede the (sub)section title (whitespace / newlines are fine).

Highlighted blocks with notes/warnings can be placed in the text as follows:

.. note:: This is a note!!

.. warning:: This is a warning!!

Similarly a 

.. todo:: Complete this section.

should be used to mark something as incomplete; a warning will appear in the
manual, and a to-do list will be generated in the table of contents, allowing
us to go back in and complete the sections that need attention.


We have provided some similar syntax for marking up input specs (similar to
what's in the current *.doc files) and keyword options.  See the
source/ref/cmd/nbon.rst file for a demonstration of their usage but, briefly,
we can can specify

/--------------------------------------------------------------\
|.. command:: name-of-charmm-command  options-expected         |
|                                                              |
|                                                              |
|  Some text describing the options expected by this command.  |
|                                                              |
\--------------------------------------------------------------/

This will automatically create an entry in the index for this command, and can
be refered to anywhere else in the text using the syntax
:cmd:`name-of-charmm_command`.  If multiple arguments are possible, they may be
combined into a "spec", as in the current docs.  There's a convenient helper
for these:

/-----------------------------------------------------------------------------------\
|.. spec: name-of-spec                                                              |
|                                                                                   |
|   .. keyword:: first-keyword any-args-or-specs-expected-by-this-keyword           |
|                                                                                   |
|      Brief description of what first_keyword will do, and any options it takes.   |
|                                                                                   |
|   .. keyword:: another-keyword any-args-or-specs-expected-by-this-keyword         |
|                                                                                   |
|      Brief description of what another_keyword will do, and any options it takes. |
|                                                                                   |
\-----------------------------------------------------------------------------------/

The spec function will also automatically create an entry in the index, and can
be linked to using the syntax :spec:`name-of-spec`.  The keyword function
formats keywords in uniform manner, creates index entries, and allows links to
be created using the syntax :kw:`name-of-keyword`.  More options are available
when specifying keywords:

/-------------------------------------------------------\
|.. keyword: keyword-name  keyword-arguments-go-here    |
|   :default: the default arguments are described here  |
|   :length: 5                                          |
|                                                       |
|   The keyword behavior is described here, as normal.  |
|                                                       |
\-------------------------------------------------------/

The default argument is entirely optional; if provided it will print a
formatted statement describing the keyword's default behavior.  The length
keyword is also optional, and is provided to account for the variability in the
number of "significant" characters in CHARMM keywords.  The default is 4.  For
example, the snippet above specifies a length of 5, so our keyword will appear
in the manual as

KEYWOrd-name

and may be referenced by in the manual using :kw:`keyword-name`, :kw:`keywo`,
or anything in between.  The label used to reference the keyword is case
insensitive, and the resulting tag will also be formatted according to the
number significant characters.



Sections and subsections
========================

ReST is very flexible about the types of characters used to underline
(sub)section titles.  However, the order that they are used determines whether
they define a section, subsection, sub-subsection, etc.  Therefore, we enforce
the following conventions.

  Section
  =======

  Subsection
  ----------

  Subsubsection
  *************

  Subsubsubsection
  ................

  Subsubsubsubsection
  ~~~~~~~~~~~~~~~~~~~

If you need even more level than this, you should probably re-think your
document structure.



Syntax Highlighting
===================

Script snippets may be included by using the code-block directive, with the
appropriate language flag:

/--------------------------------\
|                                |
|.. code-block:: chm             |
|                                |
| ! AN EXAMPLE                   |
| OPEN UNIT 11 NAME system.1.dcd |
|                                |
|Blah blah blah blah, stuff..... |
|                                |
\--------------------------------/

The "languages" defined so far are:-

  * chm - Charmm input scripts
  * crd - Charmm coordinate files
  * rtf - RTF file
  * psf - PSF file
  * prm - Parameter file

The same machinery can be used to include code snippets from an external file:

/------------------------------------------\
|                                          |
|.. literal-include:: scripts/example.inp  |
|   :language: chm                         |
|   :linenos:                              |
|                                          |
|Blah blah blah blah, stuff.....           |
|                                          |
\------------------------------------------/

The optional linenos argument will also number the lines.



Building the document
=====================

If you have the appropriate software setup, you can make the web page version
of the docs by running (from this directory)

  make html

and the resulting documentation will appear in build/html/index.html.  To get
the PDF version, you can run

  make latexpdf

and the result will be in build/latex/CHARMMDocumentation.pdf.  We've noticed
that the dependencies don't seem to be completely bulletproof, so we recommend
you remove the entire build directory before making the documentation, if
you've added a significant amount of content.


