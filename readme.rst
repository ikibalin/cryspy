CrysPy
====================

CrysPy_ is a crystallographic library for neutron data analysis. Based on the library *CrysPy* a script *RhoChi* allows to refine polarized neutron diffraction experiments performed with single crystals as well as with powder magnetic compounds. A console version is given together with a simple GUI viewer.

.. image:: cryspy/scripts/rhochi/f_icon/smm.PNG

+--------+-----------------+
| Folder | Description     |
+========+=================+
| cryspy | cryspy library  |
+--------+-----------------+
| rhochi | RhoChi script   |
+--------+-----------------+
| example| examples        |
+--------+-----------------+
| docs   | documentation   |
+--------+-----------------+


Main Features
------------------------
- Analysis of the polarized neutron scattering on crystals by the library *CryPy*;
- Diffraction data refinement for single crystals or powder by *RhoChi*;

Installation and Requirements of *CryPy*
------------------------------------------

CrysPy_ is developed and tested using Python 3.7 and depends on:

- *numpy*
- *scipy*
- *matplotlib*
- *pyqt5*.

It can be installed by

>>> python -m pip install cryspy  # as root (in Windows OS)


Or if you have downloaded a source tarball you can install CrysPy_ by doing the following:

>>> python setup.py develop # as root



Run *RhoChi*
------------------------------------------


The *CrysPy* library should be installed. 

Use the command-line to run the refinement in the folder which contents the data:

>>> python -m cryspy.run_rhochi

Or use a  simple viewer (GUI):

>>> python -m cryspy.gui_rhochi

Collaboration
---------------------------

Any third-party scripts based on the library *CrysPy* can be added.

If you have any suggestions, bug reports or annoyances please report them to our issue tracker at CrysPy_.

Copyright and License
-------------------------------

MIT License

Copyright (c) 2018-2019 Iurii Kibalin
https://github.com/ikibalin/cryspy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

.. _CrysPy: https://github.com/ikibalin/cryspy"GitHub link on CrysPy"