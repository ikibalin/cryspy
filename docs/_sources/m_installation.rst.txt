====================================
Installation
====================================

.. _CrysPy GitHub repository : https://github.com/ikibalin/cryspy
.. _numpy:  http://numpy.org/
.. _scipy:  https://scipy.org/scipylib/index.html
.. _pycifstar:  https://pypi.org/project/pycifstar/
.. _pyqt5:  https://www.riverbankcomputing.com/software/pyqt/download5
.. _matplotlib:  https://matplotlib.org/
.. _sphinx: https://www.sphinx-doc.org
.. _python:  https://python.org
.. _pytest: https://pytest.org/


Prerequisites
~~~~~~~~~~~~~~~

CrysPy works with `Python`_ versions 3.6 and higher.

CrysPy requires the following Python packages, with versions given:
   * `NumPy`_
   * `SciPy`_
   * `PyCifStar`_
   * `PyQt5`_ (only for **CrysPy Editor**)
   * `Matplotlib`_ (only for **CrysPy Editor**)

All of these are readily available on PyPI, and can be installed
with ``pip install <package name>``.

Installation
~~~~~~~~~~~~~~~~~

The latest stable version of cryspy is |release| and is available from `PyPI
<https://pypi.python.org/pypi/cryspy/>`_. 

The easiest way to install cryspy is with::

    pip install cryspy

To get the latest development version from the CrySpy `GitHub repository`_, use::

   git clone https://github.com/ikibalin/cryspy.git

and install it using::

   python setup.py develop

We welcome all contributions to **CrysPy** library. 

Checking 
~~~~~~~~~~~~~~~~~

A battery of tests scripts that can be run with the `pytest`_ testing framework
is distributed with cryspy in the ``tests`` folders. These are automatically run
as part of the development process.
For any release or any master branch from the git repository, running ``pytest``
should run all of these tests to completion without errors or failures.

To run the tests use command::

    python -m pytest 

in the directory where the **CrysPy** library was installed (../python37/Lib/cryspy).


Many of the examples in this documentation are distributed with cryspy in the
``examples`` folder, and should also run for you. 
