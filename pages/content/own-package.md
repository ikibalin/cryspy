---
layout: page
permalink: /content/own-package/
---
[Back to Content][content]

The folder with template files can be downloaded by the [link][template-package].

## 4 steps:

1. In the empty folder create the file `my_items_loops.py` in which [new items and loops][own-item-loop] should be defined (not necessary);

2. Create the file `my_data.py` in which [new data objects][own-data] should be defined (not necessary);

3. Create the file `my_global.py` in which [new global objects][own-global] should be defined (not necessary);

4. Create the file `__init__.py`:

{% highlight python %}
# import new items and loops
from .my_items_loops import MyObject, MyObjectL 
# import new data containers
from .my_data import MyData
# import new global containers
from .my_global import MyGlobal

# add all imported classes to the list `CRYSPY_CLASSES`
CRYSPY_CLASSES = [MyObject, MyObjectL, MyData, MyGlobal, ]
{% endhighlight %}

[content]: /cryspy/content
[template-package]: https://1drv.ms/u/s!AszKwJJFkEkTjMFgKigLj3raX5Vfug?e=wSJ2Ie

[own-item-loop]:  /cryspy/content/own-item-loop
[own-data]:  /cryspy/content/own-data
[own-global]:  /cryspy/content/own-global