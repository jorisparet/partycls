Radial descriptor
=================

Definition
----------

Let :math:`\mathbf{r}_i` be the position of particle :math:`i` and define :math:`r_{ij} = |\mathbf{r}_j - \mathbf{r}_i|` as the distance between particle :math:`i` and its neighbors :math:`j`. We define 
:math:`n_i(r_m)` as the number of neighbors :math:`j` of particle :math:`i` for which :math:`r_{ij}` is between :math:`r_m = m \times \Delta r` and :math:`r_{m+1} = (m+1) \times \Delta r`. Here, :math:`\Delta r` has the interpration of a bin width in a histogram.

We then consider :math:`n_i(r_m)` for a set of distances separated by :math:`\Delta r`, :math:`d_n = \{ r_\mathrm{min}, r_\mathrm{min} + r_1, r_\mathrm{min} + r_2, \dots, r_\mathrm{max} \}`.

Constructor
-----------

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.gr.RadialDescriptor.__init__

Examples
--------
