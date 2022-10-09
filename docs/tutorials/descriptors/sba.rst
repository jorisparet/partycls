Smoothed bond-angle descriptor
==============================

.. Important::
	See :doc:`ba` first.

Definition
----------

.. math::
	q_\theta^S(i) = \sum_{j=1}^N \sum_{k=1}^N f(r_{ij}, r_{ik}) \delta(\theta - \theta_{jik})

.. math::
	f(r_{ij}, r_{ik}) = \exp \left[ - \left[ ( r_{ij} / r_{\alpha\beta}^c )^\gamma + ( r_{ik} / r_{\alpha\beta}^c )^\gamma \right] \right] H( R_{\alpha\beta}^c - r_{ij} ) H(R_{\alpha\beta}^c - r_{ik})

Constructor
-----------

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.smoothed_ba.SmoothedBondAngleDescriptor.__init__

Examples
--------
