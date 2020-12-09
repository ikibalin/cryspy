---
layout: page
permalink: /content/own-global/
---
[Back to Content][content]

**The template of python script for new global container**

{% highlight python %}
from cryspy import GlobalN
from cryspy import Crystal, PdMeasL, RefineLs

class MyGlobal(GlobalN):
    """MyData class.
    """

    CLASSES_MANDATORY = (Crystal, )
    CLASSES_OPTIONAL = (PdMeasL, RefineLs)

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "my_global"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs):
        super(MyGlobal, self).__init__()

        self.__dict__["items"] = []
        self.__dict__["data_name"] = data_name

        for key, attr in self.D_DEFAULT.items():
            setattr(self, key, attr)
        for key, attr in kwargs.items():
            setattr(self, key, attr)

    def my_method():
        pass
{% endhighlight %}

[content]: /cryspy/content
