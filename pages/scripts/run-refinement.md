---
layout: page
permalink: /scripts/run-refinement/
---
[Back to scripts][scripts]


{% highlight python %}
import cryspy

# the name of the input file (in .rcif format)
f_name = "main.rcif" 

rhochi_obj = cryspy.globaln_to_file(f_name)
rhochi_obj.run_refinement()
rhochi_obj.save_to_files()

print(rhochi_obj)
{% endhighlight %}

File "main.rcif" should content information about [crystal and experiment][rho-chi]

[scripts]: /cryspy/scripts
[rho-chi]: /cryspy/rcif-format/rho-chi