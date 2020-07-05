.. _usr-script-msi:

Hooks to internally stored CHARMM Data
--------------------------------------

These are informally known as "question mark variables" because of the syntax
for retrieving them.

The following snippet will capture the last computed value for the total
energy, and store it in the user defined variable 'myener'.

.. code-block:: bash

   set myener = ?ENER

Listing
*******

Dump `grep setmsi` here


Examples
********

What are some nifty usages of ?VARs that people in the community already use?
Let's put those here.
