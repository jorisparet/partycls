Bond-orientational descriptor
-----------------------------

Bond-order parameters :cite:`steinhardt_1983` are standard measures of structure in the first coordination shell.

.. math::
	q_{lm}(i) = \frac{1}{N_b(i)} \sum_{j=1}^{N_b(i)} Y_{l m}(\hat{\mathbf{r}}_{ij})

Rotational invariants:

.. math::
	Q_{l}(i) = \left( \frac{4\pi}{2l + 1}\sum_{m=-l}^l |q_{lm}(i)|^2 \right)^{1/2}

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.bo.BondOrientationalDescriptor.__init__

References
~~~~~~~~~~

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
