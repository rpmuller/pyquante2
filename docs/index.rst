.. pyquante2 documentation master file, created by
   sphinx-quickstart on Tue Jul 23 10:19:02 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to pyquante2's documentation!
=====================================

Contents:

.. toctree::
   :maxdepth: 2

PyQuante_ (`Sourceforge Project Page`_) is an open-source suite of
programs for developing quantum chemistry methods. The program is
written in the Python programming language, but has many
"rate-determining" modules also written in C for speed. The resulting
code, though not as fast as Jaguar_, NWChem_, Gaussian_, or MPQC_, is
much easier to understand and modify. The goal of this software is not
necessarily to provide a working quantum chemistry program (although
it will hopefully do that), but rather to provide a well-engineered
set of tools so that scientists can construct their own quantum
chemistry programs without going through the tedium of having to write
every low-level routine.

.. _PyQuante: http://pyquante.sourceforge.net
.. _Jaguar: http://www.schrodinger.org
.. _NWChem: http://www.nwchem.org
.. _Gaussian: http://www.gaussian.com
.. _MPQC: http://www.mpqc.org
.. _`Sourceforge Project Page`: http://sourceforge.net/projects/pyquante

PyQuante2 is a major rewrite of the code to make the code cleaner and more 
flexible. In PyQuante2, all distribution will be done through the 
`pyquante2 github site`_

.. _`pyquante2 github site`: http://github.com/rpmuller/pyquante2

What works so far: 

* Huzinaga and HGP integral methods, both in Python and with Cython wrappers to C
* RHF, UHF wave functions
* MP2 perturbation theory
* A limited number of basis sets (STO-3G, 6-31G, 6-31G**)
* Line and contour plotting
* Basic IPython notebook support for some objects

Feel free to fork this if it interests you. The PyQuante code is still
around, so I'm not rushing through the process, I'm just taking as
much time as I feel I need to do this properly.

Requirements
============
* Python 2.7
* Numpy >= 1.7
* The test suites run without matplotlib or pyglet, but there is extra
  functionality using both.


License
=======
The software is released under the modified BSD license, which means
that everyone is free to download, use, and modify the code without
charge.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

