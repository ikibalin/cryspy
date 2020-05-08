CrysPy
====================

CrysPy_ is a crystallographic library for neutron data analysis. 
Based on the library *CrysPy* a script *RhoChi* allows to refine polarized neutron diffraction experiments 
performed with single crystals as well as with powder magnetic compounds. 
A console version is given together with a simple GUI_editor.

.. image:: cryspy/logo.png

+----------+------------------------+
| Folder   | Description            |
+==========+========================+
| cryspy   | cryspy library         |
+----------+------------------------+
| examples | examples               |
+----------+------------------------+
| docs     | html documentation     |
+----------+------------------------+


Main Features
------------------------
- Analysis of the polarized neutron scattering on crystals by the library *CrysPy*;
- Diffraction data refinement for single crystals or powders by *RhoChi*;

Installation and Requirements of *CrysPy*
------------------------------------------

CrysPy_ is developed and tested using Python 3.7 and depends on:

- *numpy*
- *scipy*
- *matplotlib* (for cryspy_editor)
- *pyqt5* (for cryspy_editor).

It can be installed by

>>> python -m pip install cryspy  # cryspy library
>>> python -m pip install cryspy_editor  # simple gui editor of cryspy library


Run *RhoChi*
------------------------------------------


The *CrysPy* library should be installed. 

Use the command-line to run the refinement given in file "main.rcif":

>>> python -m cryspy 

or 

>>> python -m cryspy input.rcif output.rcif

If the cryspy_editor is installed use the command-line to run simple GUI of *CrysPy* library:

>>> python -m cryspy_editor

or 

>>> python -m cryspy_editor main.rcif

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