Radial bond-orientational descriptor
====================================

.. Important::
	See :doc:`bo` first.

Definition
----------

The radial bond-order descriptor captures the radial dependence of bond order in the most straightforward way :cite:`boattini_2021`. It uses a radial basis of Gaussian functions of width :math:`\delta` centered on a grid of points :math:`\{d_n\}_{n=1 \dots n_\mathrm{max}}`,

.. math::
	G_n(r) = \exp{\left(-\frac{(d_n - r)^2}{2\delta^2}\right)} .

The complex radial bond-order coefficients are defined as

.. math::
  q_{l m n}^{R}(i) = \frac{1}{Z(i)} \sum_{j=1}^{N}
  G_n(r_{ij}) Y_{l m}(\hat{\mathbf{r}}_{ij}) ,

where :math:`Z(i) = \sum_{j=1}^N G_n(r_{ij})` is a normalization constant and the superscript :math:`R` indicates the radial dependence of the descriptor.
In the following, we actually use

.. math::
	G_n(r_{ij}) = \exp{\left(-\frac{(d_n - r_{ij})^\gamma}{2\delta^2}\right)} H(R_\mathrm{max} - r_{ij}) ,

where :math:`\gamma` is an integer, :math:`H` is the `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, which allows to neglect the contributions of particles further than a distance :math:`R_\mathrm{max} = d_{n_\mathrm{max}} + \sigma \times \delta` from the central particle, where :math:`d_{n_\mathrm{max}}` is the largest distance in the grid of points :math:`\{ d_n \}` and :math:`\sigma` is a skin distance.
Then, only the diagonal coefficients of the power spectrum, namely 

.. math::
	Q_{ln}^R(i) = \left( \frac{4\pi}{2l + 1} \sum_{m=-l}^l |q_{l m n}^R(i)|^2 \right)^{1/2} ,

are retained to form the descriptor of particle :math:`i` as :math:`(\dots, Q_{ln}^R(i), \dots)`, which is flattened as a vector composed of :math:`l_\textrm{max}\times n_\textrm{max}` structural features.

Constructor
-----------

.. automethod:: partycls.descriptor.radial_bo.RadialBondOrientationalDescriptor.__init__

Examples
--------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
