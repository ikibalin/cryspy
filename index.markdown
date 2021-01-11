---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
---

Welcome  to the documentation page

If [python 3][python-downloads] is installed, type in console (or terminal):

- **to install CrysPy:**

{% highlight bash %}
python -m pip install cryspy cryspy_editor

# additional libraries
python -m pip install pyqt5 scipy matplotlib numpy
{% endhighlight %}

- **to run it:**

{% highlight bash %}
python -m cryspy_editor
{% endhighlight %}


- **to upgrade:**

{% highlight bash %}
python -m pip install --upgrade cryspy cryspy_editor pycifstar --user
{% endhighlight %}

Current versions of the CrysPy library is **0.5.4**, and the CrysPy Editor is **1.5.5**. See examples on [GitHub][cryspy-examples].

CrysPy is a crystallographic library for neutron data analysis.
It allows to refine polarized
neutron diffraction experiments performed with single crystals as well
as with powder magnetic compounds. The program has been mainly developed for Rietveld analysis
of polarized neutron powder diffraction data. The program can be also used as a MEM tool
to reconstruct magnetization (or spin) density.
Time-of-flight (TOF) neutron data analysis is also available. 

# Features 

- Analysis of the polarized neutron scattering on crystals by the library CrysPy;
- Diffraction data refinement for single crystals or powder by RhoChi;
- MEM for reconstruction of magnetization density;
- Powder refinement of magnetic space groups (under development).
- TOF neutron data analysis (under development)

# Technical support

The author may provide technical support to the users of the program. If you have any
questions regarding the use of CrysPy library or "CrysPy Editor", troubles with 
its installation or running the program try the following steps:

- Read the relevant manual sections, paying particular attention to examples files (under development).
- Send an e-mail to *iurii.kibalin@cea.fr* or call by Zoom (id: *309 327 7961*).

A reference to the software use the following DOI: [10.5281/zenodo.4271777][cryspy-zenodo]

# Collaboration

Any third-party scripts based on the library CrysPy can be added.

If you have any suggestions, bug reports or annoyances please report them to our issue [tracker on GitHub][cryspy-issue].

[python-downloads]: https://www.python.org/downloads
[cryspy-zenodo]: https://zenodo.org/badge/latestdoi/178411703
[cryspy-issue]: https://github.com/ikibalin/cryspy/issues

[cryspy-examples]: https://github.com/ikibalin/cryspy/tree/examples
