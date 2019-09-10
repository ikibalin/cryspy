__author__ = 'ikibalin'
__version__ = "2019_09_10"

import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
 
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot
 


class cwind_central(QtWidgets.QMainWindow):
    def __init__(self):
        super(cwind_central, self).__init__()
        self.title = "program 'Graph'"
        self.setWindowTitle(self.title)
        widg_central = cwidg_central()
        self.setCentralWidget(widg_central)
        self.show()



class cwidg_central(QtWidgets.QWidget):
    def __init__(self):
        super(cwidg_central, self).__init__()
        self.init_layout_central()
        self.setLayout(self.layout_central)
 
    def init_layout_central(self):
        
        lay_main = QtWidgets.QVBoxLayout()
        self.graph = Graph(self, width=5, height=4)
        

        lay_1 = QtWidgets.QVBoxLayout()
        lay_1.addWidget(self.graph)

        lay_main.addLayout(lay_1)

        
        self.layout_central = lay_main
        
    def plot_file(self, x, l_y_mod, l_y_exp, l_y_sig):
        self.graph.ax_pri.cla()
        
        self.graph.data_x = x
        self.graph.data_l_y_mod = l_y_mod # len(l_y_mod) == len(l_y_exp) == len(l_y_sig)
        self.graph.data_l_y_exp = l_y_exp
        self.graph.data_l_y_sig = l_y_sig

        self.graph.set_data_to_graph()


class Graph(matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg):
    
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi = dpi)
        fig.subplots_adjust(left = 0.07,
                            right = 0.97,
                            top = 0.97,
                            bottom = 0.07,
                            wspace = 0.0,
                            hspace = 0.0)
        
        super(Graph, self).__init__(fig)
        
        self.info_press = (False, False)

        self.data_x = []
        self.data_l_y_mod = []
        self.data_l_y_exp = []
        self.data_l_y_sig = []
        
        self.control = parent
        self.figure = fig
        
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.ax_pri = fig.add_subplot(111)
        
        fig.canvas.mpl_connect("button_press_event", self.onclick)
        
    def set_data_to_graph(self):
        col_1 = "#FF0000"
        col_2 = "#0000FF"
        col_3 = "#006400"

        x = self.data_x
        l_y_mod = self.data_l_y_mod
        l_y_exp = self.data_l_y_exp
        l_y_sig = self.data_l_y_sig
        l_colors = [col_1, col_2]
        l_n = max([l_y_mod[0].max(), (l_y_exp[0]+3*l_y_sig[0]).max(), (l_y_exp[0]-l_y_mod[0]+3*l_y_sig[0]).max()])
        for y_mod, y_exp, y_sig, col in zip(l_y_mod, l_y_exp, l_y_sig, l_colors):
            y_max =  max([y_mod.max(), (y_exp+3*y_sig).max(), (y_exp-y_mod+3*y_sig).max()])
            y_min =  min([y_mod.min(), (y_exp-3*y_sig).min(), (y_exp-y_mod-3*y_sig).min()])
            y_diff = y_max - y_min
            self.ax_pri.errorbar(x, y_exp+l_n-y_max, yerr = y_sig, ecolor = col_1,  fmt='.', color=col, linewidth = 0.5)
            self.ax_pri.plot(x, y_mod+l_n-y_max, "k-", linewidth=1.0)
            self.ax_pri.plot(x, y_exp-y_mod+l_n-y_max, "k-",linewidth=0.5)
            l_n -=  y_diff



        #ly1_before = self.data_4_y_exp
        #y_diff_max = max(ly1_before)
        #y_diff_min = min(ly1_before) - y_diff_max + y_down_min
        x_min = min(x)
        x_max = max(x)
        x_diff = x_max - x_min
        
        y_max = max([l_y_mod[0].max(), (l_y_exp[0]+3*l_y_sig[0]).max(), (l_y_exp[0]-l_y_mod[0]+3*l_y_sig[0]).max()])
        self.xlim_orig = (x_min-0.05*x_diff, x_max+0.05*x_diff)
        self.ylim_orig = (l_n, y_max)
        """
        self.ax_pri.plot(self.data_1_x, self.data_1_y_mod, "k-", linewidth=1.0)    
         
        self.ax_pri.errorbar(self.data_1_x, self.data_1_y_exp, yerr = self.data_1_sy_exp, ecolor = col_1, fmt='.', color=col_1, linewidth = 0.5)
        """ 
        self.draw()
        
    def onclick(self, event):
        if event.button == 3:
            
            if self.info_press == (False, False):
                self.info_press = (True, False)
                self.xlim = [event.xdata]
                self.ylim = [event.ydata]
            elif self.info_press == (True, False):
                self.info_press = (True, True)
                self.xlim.append(event.xdata)
                self.ylim.append(event.ydata)
            if self.info_press == (True, True):
                self.info_press = (False, False)
                xlim = (min(self.xlim), max(self.xlim))
                ylim = (min(self.ylim), max(self.ylim))
                self.ax_pri.set_xlim(xlim)
                self.ax_pri.set_ylim(ylim)
                self.xlim = []
                self.ylim = []
                self.draw()
        elif event.button == 2:
                self.info_press == (False, False)
                self.ax_pri.set_xlim(self.xlim_orig)
                self.ax_pri.set_ylim(self.ylim_orig)
                self.xlim = []
                self.ylim = []
                self.draw()
        else:
            self.info_press == (False, False)
            self.xlim = []
            self.ylim = []
            
 
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    fname = "NaCaAlF_exp.out"
    ex = cwind_central()
    
    sys.exit(app.exec_())
