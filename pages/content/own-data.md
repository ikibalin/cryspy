---
layout: page
permalink: /content/own-data/
---
[Back to Content][content]

**The template of python script for new data container**

{% highlight python %}
from cryspy import DataN
from cryspy import Cell, SpaceGroup, PdPeakL, RefineLs

class MyData(DataN):
    """MyData class.
    """

    CLASSES_MANDATORY = (Cell, SpaceGroup, PdPeakL)
    CLASSES_OPTIONAL = (RefineLs, )

    CLASSES = CLASSES_MANDATORY + CLASSES_OPTIONAL

    PREFIX = "my_data"

    # default values for the parameters
    D_DEFAULT = {}

    def __init__(self, data_name=None, **kwargs):
        super(MyData, self).__init__()

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
