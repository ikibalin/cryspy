---
layout: page
permalink: /content/add-packages/
---
[Back to Content][content]

# Add packages to the CrysPy library

[The own packages][own-package] can be added into the CrysPy library by the script (it's enough to do one time):

{% highlight python %}
import cryspy

folder = "C:\my_package"
cryspy.add_package(folder)
{% endhighlight %}

To load new objects into the library use the method `.load_packages()`:

{% highlight python %}
import cryspy

cryspy.load_packages()

# access to object defined in the package "my_package"
cryspy.MyObject() 
{% endhighlight %}

[content]: /cryspy/content
[own-package]:  /cryspy/content/own-package