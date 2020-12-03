---
layout: page
permalink: /usage/
---

The CrysPy library can be used in two ways: 
(i) as a python script (or jupyter notebook) or
(ii) the graphical user interface: "Cryspy Editor". 

To run the program in a console (or terminal) the user has to invoke the name of the executable file or an
appropriate alias, for instance:


# To run GUI "CrysPy Editor"

Use console/terminal to run GUI "CrysPy Editor":

{% highlight bash %}
python -m cryspy_editor
{% endhighlight %}

or put it into bash file (in Windows os) and run it after.  

# To run refinement in console mode

To run data refinement given in the file `file_name.rcif` print in console mode :

{% highlight bash %}
python -m cryspy file_name.rcif
{% endhighlight %}

# Example of python script or Jupyter notebook 

{% highlight python %}
import cryspy

# the name of the input file (in .rcif format)
f_name = "main.rcif" 

rhochi_obj = cryspy.globaln_to_file(f_name)
rhochi_obj.run_refinement()
rhochi_obj.save_to_files()

print(rhochi_obj)
{% endhighlight %}

