.. index:: fix mbx

    fix mbx command
    ===============

    Syntax
    """"""

    .. code-block:: LAMMPS
        fix ID group-ID mbx num_mon_types ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mbx = style name of this fix command
* num_mon_types = number of monomer types

  .. parsed-literal::

    *json*
    *print/dipoles*
    *print/settings*


Examples
""""""""

.. code-block:: LAMMPS

    # For a system involving ch4 (atom types C=1, H=2) and
    # water (atom types O=3, H=4)
    fix 1 all mbx 2 ch4 1 2 5 1 2 2 2 2 h2o 3 4 3 3 4 4 json mbx.json

Description
"""""""""""

This fix instructs LAMMPS to call the `MBX library <_mbxwebsite>`
in order to simulate MB-nrg models such as MB-pol. This fix 

The MBX library code development is available at
`https://github.com/paesanilab/MBX <https://github.com/paesanilab/MBX>`.

If you have questions not answered by this documentation, please reach
out to us at `https://groups.google.com/g/mbx-users <https://groups.google.com/g/mbx-users>`



Restrictions
""""""""""""

This fix is part of the MBX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info. This fix also relies on the
presence of `pair mbx <pair_mbx>` command.

There can only be one fix mbx command active at a time.

Due to the usage of Partridge and Schwenke charges for MB-pol,
all electrostatic interactions are calculated internally in MBX.
Therefore one should avoid calculating coulombic interactions in
LAMMPS such as using `coul/cut` or `coul/long`.



Related commands
""""""""""""""""

:doc:`pair mbx <pair_mbx>`



-----------

.. _mbxwebsite: https://mbxsimulations.com