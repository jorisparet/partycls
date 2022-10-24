Installation
============

The package is compatible with Python >= 3.6 (tested on 3.6, 3.7, 3.8 and 3.9). We currently only support Unix-based systems.

From PyPI
---------

The latest stable release is available on `PyPI <https://pypi.org/project/partycls/>`_. Install it with ``pip``:

.. code-block:: sh

	pip install partycls

From source
-----------

To install the latest development version from source, clone the source code from the official `GitHub repository <https://github.com/jorisparet/partycls>`_ and install it with:

.. code-block:: sh

	git clone https://github.com/jorisparet/partycls.git
	cd partycls
	pip install -r requirements.txt
	make install

Run the tests using:

.. code-block:: sh
	
	make test

or manually compile the Fortran sources and run the tests:

.. code-block:: sh

	cd partycls/
	f2py -c -m neighbors_wrap neighbors.f90
	cd descriptor/
	f2py -c -m realspace_wrap realspace.f90
	cd ../../
	pytest tests/

Dependencies
------------

partycls relies on several external packages, most of which only provide additional features and are not necessarily required.

Required
~~~~~~~~

- Fortran compiler (*e.g.* `gfortran <https://gcc.gnu.org/wiki/GFortran>`_)
- `NumPy <https://pypi.org/project/numpy/>`_
- `scikit-learn <https://scikit-learn.org>`_

Optional
~~~~~~~~

- `MDTraj <https://www.mdtraj.org>`_ (additional trajectory formats)
- `atooms <https://framagit.org/atooms/atooms>`_ (additional trajectory formats)
- `DScribe <https://singroup.github.io/dscribe>`_ (additional descriptors)
- `Matplotlib <https://matplotlib.org/>`_ (visualization)
- `OVITO <https://ovito.org/>`_ < 3.7.0 (visualization)
- `Py3DMol <https://github.com/avirshup/py3dmol>`_ (interactive 3D visualization)
- `pyvoro <https://github.com/joe-jordan/pyvoro>`_ or its `memory-optimized fork <https://framagit.org/coslo/pyvoro>`_ for large systems (Voronoi neighbors and tessellation)
- `tqdm <https://tqdm.github.io/>`_ (progress bars)
