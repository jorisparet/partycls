Bond-orientational descriptor
=============================

Definition
----------

Bond-order parameters :cite:`steinhardt_1983` are standard measures of structure in the first coordination shell. Let :math:`\mathbf{r}_i` be the position of particle :math:`i` and define :math:`\mathbf{r}_{ij} = \mathbf{r}_j - \mathbf{r}_i` and :math:`r_{ij} = |\mathbf{r}_{ij}|`. Then consider the weighted microscopic density around particle :math:`i`:

.. math::
	\rho(\mathbf{r}; i) = \sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \mathbf{r}_{ij})

where :math:`w_j` is a particle-dependent weight and the sum involves a set of :math:`N_b(i)` particles, which defines the coordination shell of interest for particle :math:`i`.

We project the microscopic density on a unit-radius sphere, that is, :math:`\hat{\rho}(\hat{\mathbf{r}}; i) = \sum_{j=1}^{N_b(i)} w_j \delta(\mathbf{r} - \hat{\mathbf{r}}_{ij})`,
where :math:`\hat{\mathbf{r}} = \mathbf{r} / |\mathbf{r}|` and similarly :math:`\hat{\mathbf{r}}_{ij} = \mathbf{r}_{ij}/|\mathbf{r}_{ij}|`. Expanding in spherical harmonics yields

.. math::
	\hat{\rho}(\hat{\mathbf{r}}; i) = \sum_{l=0}^\infty \sum_{m=-l}^l c_{l m}(i) Y_{l m}(\hat{\mathbf{r}}) ,

with coefficients

.. math::
	c_{l m}(i) =  \int d\mathbf{r} \rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) .

In the conventional bond-order analysis, one sets the weights :math:`w_j` to unity and considers the normalized complex coefficients,

.. math::
	\begin{align}
	q_{lm}(i) & = \frac{1}{N_b(i)} \int d\mathbf{r} \rho(\mathbf{r}; i) Y_{l m}(\hat{\mathbf{r}}) 
	\nonumber \\ & = \frac{1}{N_b(i)} \sum_{j=1}^{N_b(i)} Y_{l m}(\hat{\mathbf{r}}_{ij}) .
	\end{align}

The rotational invariants,

.. math::
	Q_{l}(i) = \sqrt{ \frac{4\pi}{2l + 1}\sum_{m=-l}^l |q_{lm}(i)|^2 },

provide a detailed structural description of the local environment around particle :math:`i`.


We then consider :math:`Q_l(i)` for a sequence of orders :math:`\{ l_n \} = \{ l_\mathrm{min}, \dots, l_\mathrm{max} \}`. The resulting feature vector for particle :math:`i` is given by

.. math::
	X^\mathrm{BO}(i) = (\: Q_{l_\mathrm{min}}(i) \;\; \dots \;\; Q_{l_\mathrm{max}}(i) \:) .

Setup
-----

Instantiating this descriptor on a ``Trajectory`` can be done as follows:

.. code-block:: python

	from partycls import Trajectory
	from partycls.descriptor import BondOrientationalDescriptor

	traj = Trajectory("trajectory.xyz")
	D = BondOrientationalDescriptor(traj)

The constructor takes the following parameters:

.. automethod:: partycls.descriptor.bo.BondOrientationalDescriptor.__init__

Demonstration
-------------

References
----------

.. bibliography:: ../../references.bib
	:style: unsrt
	:filter: docname in docnames
