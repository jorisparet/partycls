Tetrahedral descriptor
======================

Definition
----------

.. math::
	T(i) = \frac{1}{N_\mathrm{ba}} \sum_{i \in N_b(i)} \sum_{j \in N_b(i)} | \cos(\theta_{jik}) - \cos(109.5) |

Constructor
-----------

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.tetrahedrality.TetrahedralDescriptor.__init__

Examples
--------
