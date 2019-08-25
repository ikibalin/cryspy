# NeuPy 
***

[NeuPy] is a crystallographic library for neutron data analysis. Based on the library *NeuPy* a script *RhoChi* allows to refine polarized neutron diffraction experiments performed with single crystals as well as with powder magnetic compounds. A console version is given together with a simple GUI viewer. 

<img src="rhochi/f_icon/smm.png">



 Folder    |  Description
 -----     |  :--------
 neupy     |  NeuPy library 
 rhochi    |  RhoChi script 
 example   |  examples
 docs      |  documentation

## Main Features


 - Analysis of the polarized neutron scattering on crystals by the library *NeuPy*;
 - Diffraction data refinement for single crystals or powder by *RhoChi*;
 
## Installation and Requirements of *NeuPy*

NeuPy is developed and tested using Python 3.7 and depends on:

- *numpy* 
- *scipy*
- *matplotlib*
- *pyqt5*.

They can be installed by
```
python -m pip install numpy scipy matplotlib pyqt5  # as root (in Windows OS)
```

If you have downloaded a source tarball you can install NeuPy by doing the following:
```
python setup.py develop # as root 
```

And then just import a library in a python code:
```
import neupy
```


## Run *RhoChi*

The *NeuPy* library should be installed.

### As a script (CLI) 
Use the command-line to run the refinement. E.g.:
```
python ./rhochi.py ./full.rcif ./full.rcif
```

### Simple Viewer (GUI) 
Run the file `rhochi_viewer.pyw`.

## Collaboration

Any third-party scripts based on the library *NeuPy* can be added.

If you have any suggestions, bug reports or annoyances please report them to our issue tracker at [site][NeuPy].

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

[NeuPy]: https://github.com/ikibalin/neupy "GitHub link on NeuPy"