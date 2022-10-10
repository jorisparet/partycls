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
	Q_{l,n}^R(i) = \sqrt{ \frac{4\pi}{2l + 1} \sum_{m=-l}^l |q_{l m n}^R(i)|^2 } ,

are retained to form the descriptor of particle :math:`i`.

We then consider :math:`Q^R_{l,n}(i)` for a sequence of orders :math:`\{ l_m \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}` and for a grid of distances :math:`\{ d_n \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{RBO}(i) = (\: \dots \;\; Q^R_{l,n}(i) \;\; Q^R_{l, n+1}(i) \;\; \dots \;\; Q^R_{l_\mathrm{max}, n_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import RadialBondOrientationalDescriptor

	traj = Trajectory("trajectory.xyz")
	D = RadialBondOrientationalDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.radial_bo.RadialBondOrientationalDescriptor.__init__

Demonstration
-------------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
