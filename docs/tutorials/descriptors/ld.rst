Lechner-Dellago descriptor
==========================

.. Important::
	See :doc:`bo` first.

Definition
----------

The complex coefficients :math:`q_{lm}(i)` of particle :math:`i` can be averaged over nearest neighbors, as suggested by Lechner and Dellago :cite:`lechner_2008`,

.. math::
	\bar{q}_{lm}(i) \equiv \langle q_{l m}(i) \rangle ,

and then made invariant,

.. math::
	\bar{Q}_{l}(i) = \sqrt{ \frac{4\pi}{2l + 1}\sum_{m=-l}^l |\bar{q}_{lm(i)}|^2 } ,

to provide an improved descriptor for crystal structure detection.

We then consider :math:`\bar{Q}_l(i)` for a sequence of orders :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{LD}(i) = (\: \bar{Q}_{l_\mathrm{min}}(i) \;\; \dots \;\; \bar{Q}_{l_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import LechnerDellagoDescriptor

	traj = Trajectory("trajectory.xyz")
	D = LechnerDellagoDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.bo.LechnerDellagoDescriptor.__init__

Demonstration
-------------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
