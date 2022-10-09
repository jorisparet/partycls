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
	\bar{Q}_{l}(i) \equiv \left( \frac{4\pi}{2l + 1}\sum_{m=-l}^l |\bar{q}_{lm(i)}|^2 \right)^{1/2} ,

to provide an improved descriptor for crystal structure detection.

Constructor
-----------

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.bo.LechnerDellagoDescriptor.__init__

Examples
--------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
