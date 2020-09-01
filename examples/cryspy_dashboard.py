from bokeh.plotting import figure, output_file, show
from bokeh.io import output_notebook
from bokeh.layouts import column, row
from bokeh.models import RangeTool
from bokeh.models.widgets.inputs import FileInput 
from bokeh.models import ColumnDataSource, CustomJS, Slider

from bokeh.models.widgets import TextInput, Button, Paragraph
import numpy

output_file("cryspy_dashboard.html")

np_x = numpy.linspace(0,3,100)
np_y = numpy.square(np_x)


# create some widgets
button = Button(label="Say HI")
input = TextInput(value="Bokeh")
output_p = Paragraph(text="test text")

# add a callback to a widget
    
fi = FileInput(accept=".rcif, .cif")

p1 = figure()

def update():
    output_p.text = "Hello, " + input.value + "filename: " + fi.filename

button.on_click(update)


# create a layout for everything

p1.line(np_x, np_y)



x = [x*0.005 for x in range(0, 200)]
y = x

source = ColumnDataSource(data=dict(x=x, y=y))

plot = figure(plot_width=400, plot_height=400)
plot.line('x', 'y', source=source, line_width=3, line_alpha=0.6)

callback = CustomJS(args=dict(source=source), code="""
    var data = source.data;
    var f = cb_obj.value
    var x = data['x']
    var y = data['y']
    for (var i = 0; i < x.length; i++) {
        y[i] = Math.pow(x[i], f)
    }
    source.change.emit();
""")

slider = Slider(start=0.1, end=4, value=1, step=.1, title="power")
slider.js_on_change('value', callback)

layout = column(slider, plot)

show(column(fi, p1, button, input, output_p, layout))

