.. index:: fix graphics/labels

fix graphics/labels command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/label Nevery mode keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/label = style name of this fix command
* Nevery = update graphics information every this many time steps
* zero or more keyword/args groups may be appended
* keyword = *image* or *text*

  .. parsed-literal::

     *image* filename x y z keyword args = display image in visualization
        filename = name of the image file
        x, y, z  = position where the center of the image is located in the visualization
        any of x, y, or z can be a variable (see below)
        one or more keyword/arg pairs may be appended
        keyword = *scale* or *transcolor*
          *scale* value = the image is scaled by this value (default 1.0), can be a variable (see below)
          *transcolor* arg = select color for transparency: *auto* or *none* or *r/g/b*
             *auto* = uses the color in the lower left corner of the image for transparency
             *none* = disables transparency
             *r/g/b* = provide three integers in the range 0 to 255 to select transparancy color in RGB color space
        
     *text* labeltext x y z keyword args = display text in visualization
        labeltext = text for the label, must be quoted if it contains whitespace
        x, y, z  = position where the center of the text is located in the visualization
        any of x, y, or z can be a variable (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix pix all graphics/label 100 image teapot.ppm 5.0 -1.0 -2.0 transcolor auto scale 0.75
   fix lbl all graphics/label 1000 text "LAMMPS Graphics Demo" 5.0 -1.0 -2.0

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add either images or text to :doc:`dump image
<dump_image>` images using the *fix* keyword.  This can be useful to
augment images directly with additional graphics or text directly and
without having to post-process the images.  Since the positioning uses
the coordinate system of the simulation and because the graphics
object use the depth buffer of the image rasterizer, atoms can be
located before or behind any text or image.

The *group-ID* is ignored by this fix.

The *Nevery* keyword determines how often the arrows graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

The *image* keyword...

The filename suffix determines whether LAMMPS will try to read a file in
JPEG, PNG, or PPM format.  If the suffix is ".jpg" or ".jpeg", then
LAMMPS attempts to read the image in `JPEG format <jpeg_format_>`_, if
the suffix is ".png", then LAMMPS attempts to read the image in `PNG
format <png_format_>`_.  Otherwise LAMMPS will try to read the image in
`PPM (aka NETPBM) format <ppm_format_>`_.

The *text* keyword...

There may be multiple *image* or *text* keywords with their arguments
in a single fix *graphics/labels* command.

The arguments for the position or scale of an *image* or *text* can be
specified as an equal-style :doc:`variable <variable>`, namely *x*, *y*,
*z*, or the *scale* value for an image.  If any of these values is a
variable, it should be specified as `v_name`, where `name` is the
variable name.  In this case, the variable will be evaluated each
*Nevery* timestep, and its value used to define the indenter geometry.
Please see the documentation of the :doc:`fix graphics <fix_graphics>`
command for more details on using variables with graphics objects.

.. _jpeg_format: https://jpeg.org/jpeg/
.. _png_format: https://en.wikipedia.org/wiki/Portable_Network_Graphics
.. _ppm_format: https://en.wikipedia.org/wiki/Netpbm

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/label is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix adds images or text to the
visualization.

The color of the text is either that of the first atom type when using
color styles "type" or "element" with the fix command.  With fix color
style "const" the default value of "white" can be changed using
:doc:`dump_modify fcolor <dump_image>`.  The transparency is by default
fully opaque and can be changed with *dump\_modify ftrans*\ .

The *fflag1* and *fflag2* settings of *dump image fix* are ignored.

Restart, fix_modify, output, run start/stop, minimize info
==========================================================

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

To read JPEG or PNG format images, support for the corresponding
graphics libraries must have been compiled and linked into LAMMPS.
Please see the :ref:`instructions for building LAMMPS <graphics>` for
more information on how to do that.

Related commands
""""""""""""""""

:doc:`fix graphics <fix_graphics>`, :doc:`fix graphics/arrows <fix_graphics_arrows>`

Default
"""""""

radius = auto, atoms = yes, bonds = yes, if supported by atom style otherwise no,
no label graphics
