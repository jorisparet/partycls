Smoothed bond-orientational descriptor
======================================

.. Important::
	See :doc:`bo` first.

Definition
----------

This is a smooth version of the :doc:`bo`, in which the coefficients :math:`q_{lm}(i)` are multiplied by a weighting function :math:`f(r)` that depends on the radial distance :math:`r` between the central particle :math:`i` and its neighbors,

.. math::
	q_{lm}^{S}(i) = \frac{1}{Z(i)} \sum_{j=1}^{N} f({r}_{ij}) Y_{lm}(\hat{\mathbf{r}}_{ij}) ,


where :math:`Z(i)=\sum_{j=1}^{N} f({r}_{ij})` is a normalization constant and the superscript :math:`S` indicates the smooth nature of the descriptor. We use

.. math::
	f(r_{ij}) = \exp \left[- (r_{ij} / r_{\alpha\beta}^c)^\gamma \right] H(R_{\alpha\beta}^c - r_{ij}) ,


where :math:`r_{\alpha\beta}^c` is the first minimum of the corresponding partial radial distribution function for the pair :math:`(i,j)` and :math:`\gamma` is an integer.
Also, :math:`H` is the `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which ensures, for efficiency reasons, that the descriptor only has contributions from particles within a distance :math:`R_{\alpha\beta}^c = \xi \times r_{\alpha\beta}^c` from the central one, where :math:`\xi > 1` is a scaling factor.

The rotational invariants are defined similarly to the :doc:`bo`.

Constructor
-----------

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.smoothed_bo.SmoothedBondOrientationalDescriptor.__init__

Examples
--------
