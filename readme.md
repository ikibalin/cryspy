# RhoChi 
***

[RhoChi] is a software for polarized neutron diffraction analysis. It allows to work with single crystals or powder magnetic compounds. A console version is given together with a simple GUI viewer. 

![logo](rhochi/f_icon/smm.png "logo")



 Folder    |  Description
 -----     |  :--------
 neupy     |  NeuPy library 
 rhochi    |  RhoChi script 
 mem       |  MEM script
 example   |  examples
 docs      |  documentation

## Main Features


 - Crystallographical calculations by library *NeuPy*;
 - Diffraction data analysis for single crystals or powder by script *RhoChi*;
 - MEM reconstruction by script *MEM*.

## Installation and Requirements of *NeuPy*


NeuPy is developed and tested using Python 3.7 and depends on:

- *numpy* 
- *scipy*
- *matplotlib*
- *pyqt5*.

They can be installed by
```
python -m pip install numpy scipy matplotlib pyqt5  # as root 
```

If you have downloaded a source tarball you can install NeuPy by doing the following:
```
python setup.py develop # as root 
```

Or you can install NeuPy using pip like this (not realized):
```
pip install neupy
```


And then just import a library:
```
import neupy
```


## Run RhoChi

### As a script (CLI) 
Use the command-line to run the refinement. E.g.:
```
python ./rhochi.py ./full.rcif ./full.rcif
```

### Simple Viewer (GUI) 
Run the file `rhochi_viewer.pyw`

## Collaboration

Any third-party scripts based on the library *NeuPy* can be added.

If you have any suggestions, bug reports or annoyances please report them to our issue tracker at [site][RhoChi].

## Copyright and License

MIT License

Copyright (c) 2018-2019 Iurii Kibalin   
https://github.com/ikibalin/rhochi

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

[RhoChi]: https://github.com/ikibalin/rhochi "GitHub link on RhoChi"