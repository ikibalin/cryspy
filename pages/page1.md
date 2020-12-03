---
layout: page
title: Installation/Upgrade
permalink: /installation/
---


CrysPy is developed and tested using Python 3.7 and depends on:

- numpy
- scipy
- matplotlib
- pyqt5
- pyqtgraph
- pycifstar

# Installation

The CrysPy library can be install through pip:

{% highlight bash %}
python -m pip install cryspy cryspy_editor  # as root (in Windows OS)
{% endhighlight %}

If it's needed install missing libraries

{% highlight bash %}
python -m pip install numpy scipy matplotlib pyqt5 pyqtgraph pycifstar
{% endhighlight %}


# Upgrade

Print in the conosole/ terminal to upgrade the packages
{% highlight bash %}
python -m pip install --upgrade cryspy cryspy_editor --user
{% endhighlight %}
