.. index:: SELE, selec, select

.. _ref-select:

Atom Selection
==============

Atom selection is used for many commands within CHARMM.  Its existance is one
of the main factors in the versatility of CHARMM.

Syntax
------

Atom selection can be performed by the following CHARMM command

.. code-block:: chm

   SELEct ... END

The ``SELEct ... END`` construct taken together defines the :spec:`atom-selection`.

Atom selection is ubiquitous in CHARMM and can be used with commands such as:
:cmd:`vibr`, :cmd:`ener`, :cmd:`test` and many others.
   
.. spec:: atom-selection

   For further helpful discussion please see the :ref:`user guide
   <usr-script-io-select>` for more information.
