.. _con_analysis:

.. index:: Analysis

Analysis
========

In a MD simulation, we follow the atomic motions of a system (such as a
solvated protein) over time. Normally, you save this information (the
coordinates at a specific timestep) in regular intervals, generating
*trajectories*. These trajectories are the key for later analysis. In
passing, we note that in addition to coordinates one sometimes also saves the
velocities for later analysis.

One obvious thing to do with trajectories is to watch a "molecular movie",
using programs such as VMD. Occasionally, watching such a movie may provide
interesting insight; most of the time, however, one is left with "information
overload"; the human mind is not capable of analyzing in raw form the
coordinates of several thousands atoms (or more).

Correlation Analysis
--------------------


Although the theory should be described here, there should be quite a detailed hands-on tutorial in the tut section, too.

.. todo:: This is to be completed by :ref:`developers-ln` and/or :ref:`developers-mc`
