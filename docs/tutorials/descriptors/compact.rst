Compactness descriptor
======================

Definition
----------

Reference :cite:`tong_2018`.

.. math::
	\omega_{\langle ijkm \rangle} = \frac{ \sum_{\langle ab \rangle} | r_{ab} - \sigma_{ab} |}{\sum_{\langle ab \rangle} \sigma_{ab}}

.. math::
	\Omega(i) = \frac{1}{N_\mathrm{tetra}(i)} \sum_{\langle ijkm \rangle} \omega_{\langle ijkm \rangle}

Constructor
-----------

.. warning::
	This descriptor uses particles' radii. They must thus be set beforehand or read from the input trajectory file.

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.compactness.CompactnessDescriptor.__init__

Examples
--------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
