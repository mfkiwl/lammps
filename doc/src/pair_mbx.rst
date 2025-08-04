.. index:: pair_style mbx

pair_style mbx command
======================

Syntax
""""""

.. code-block:: LAMMPS

    pair_style mbx cutoff

* cutoff = real-space cutoff for MBX. 9.0 is usually a safe value.


Examples
""""""""

.. code-block:: LAMMPS
    pair_style mbx 9.0
    compute         mbx all pair mbx


    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    compute         mbx all pair mbx
