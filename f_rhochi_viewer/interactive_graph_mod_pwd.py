import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
 
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot
 

def read_pwd_ddata(ffile):
    fid = open(ffile, 'r')
    lcontent = fid.readlines()
    fid.close()
    lcontent = [hh1 for hh1 in lcontent if hh1[0] != "#"]
    lname=lcontent[0].strip().split()
    ddata = {}
    for name in lname:
        ddata[name]=[]
    for line in lcontent[1:]:
        if (line.strip()==''): break
        for name,sval in zip(lname, line.strip().split()):
            ddata[name].append(float(sval))
    return ddata



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
        
    def plot_file(self, ffig_full):
        self.graph.ax_pri.cla()
        ddata = read_pwd_ddata(ffig_full)

        self.graph.data_1_x = ddata["ttheta"]
        self.graph.data_1_y_exp = ddata["exp_up"]
        self.graph.data_1_sy_exp = ddata["s_exp_up"]
        self.graph.data_1_y_mod = ddata["mod_up"]
        self.graph.data_1_y_diff = [hh1 - hh2 for hh1, hh2 in zip(ddata["mod_up"], ddata["exp_up"])]
        
        
        self.graph.data_2_y_exp = ddata["exp_down"]
        self.graph.data_2_sy_exp = ddata["s_exp_down"]
        self.graph.data_2_y_mod = ddata["mod_down"]
        self.graph.data_2_y_diff = [hh1 - hh2 for hh1, hh2 in zip(ddata["mod_down"], ddata["exp_down"])]


        self.graph.data_3_x = ddata["ttheta"]
        self.graph.data_3_y_exp = [hh1 + hh2 for hh1, hh2 in zip(ddata["exp_up"], ddata["exp_down"])]
        self.graph.data_3_sy_exp = [(hh1**2 + hh2**2)**0.5 for hh1, hh2 in zip(ddata["s_exp_up"], ddata["s_exp_down"])]
        self.graph.data_3_y_mod = [hh1 + hh2 for hh1, hh2 in zip(ddata["mod_up"], ddata["mod_down"])]
        self.graph.data_3_y_diff = [hh1 + hh2 for hh1, hh2 in zip(self.graph.data_3_y_exp, self.graph.data_3_y_mod)]
        

        self.graph.data_4_x = ddata["ttheta"]
        self.graph.data_4_y_exp = [hh1 - hh2 for hh1, hh2 in zip(ddata["exp_up"], ddata["exp_down"])]
        self.graph.data_4_sy_exp = [(hh1**2 + hh2**2)**0.5 for hh1, hh2 in zip(ddata["s_exp_up"], ddata["s_exp_down"])]
        self.graph.data_4_y_mod = [hh1 - hh2 for hh1, hh2 in zip(ddata["mod_up"], ddata["mod_down"])]
        self.graph.data_4_y_diff = [hh1 - hh2 for hh1, hh2 in zip(self.graph.data_4_y_exp, self.graph.data_4_y_mod)]

        
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

        self.data_1_x = []
        self.data_1_y_exp = []
        self.data_1_sy_exp = []
        self.data_1_y_mod = []
        self.data_1_y_diff = []

        self.data_2_y_exp = []
        self.data_2_sy_exp = []
        self.data_2_y_mod = []
        self.data_2_y_diff = []


        self.data_3_x = []
        self.data_3_y_exp = []
        self.data_3_sy_exp = []
        self.data_3_y_mod = []
        self.data_3_y_diff = []
        
        self.data_4_x = []
        self.data_4_y_exp = []
        self.data_4_sy_exp = []
        self.data_4_y_mod = []
        self.data_4_y_diff = []

        self.data_5_x =  []
        
        self.control = parent
        self.figure = fig
        
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.ax_pri = fig.add_subplot(111)
        
        fig.canvas.mpl_connect("button_press_event", self.onclick)
        
    def set_data_to_graph(self):
        col_1 = "#FF0000"
        col_2 = "#0000FF"
        col_3 = "#006400"
        lx = self.data_3_x
        ly1 = self.data_3_y_exp
        lsy1 = self.data_3_sy_exp
        y_up_min = min(ly1)
        ly1_mod = self.data_3_y_mod
        ly1_exp_mod = [-hh1+hh2+y_up_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax_pri.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax_pri.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)
        self.ax_pri.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1,  fmt='.', color=col_1, linewidth = 0.5)



        ly1_before = self.data_4_y_exp
        lsy1 = self.data_4_sy_exp
        y_down_max = max(ly1_before)
        y_down_min = min(ly1_before) - y_down_max + y_up_min
        ly1 = [hh - y_down_max + y_up_min  for hh in ly1_before]
        ly1_mod = [hh - y_down_max + y_up_min  for hh in self.data_4_y_mod]
        ly1_exp_mod = [-hh1+hh2 + y_down_min for hh1, hh2 in zip(ly1, ly1_mod)]
        
        #linestyle='-'
        self.ax_pri.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax_pri.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)
        self.ax_pri.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, fmt='.', color=col_2, linewidth = 0.5)


        ly1_before = self.data_4_y_exp
        y_diff_max = max(ly1_before)
        y_diff_min = min(ly1_before) - y_diff_max + y_down_min

        x_diff = max(lx)-min(lx)
        y_diff = max(self.data_4_y_exp) - y_diff_min
        self.xlim_orig = (min(lx)-0.05*x_diff, max(lx)+0.05*x_diff)
        self.ylim_orig = (y_diff_min-0.05*y_diff, max(self.data_3_y_exp)+0.05*y_diff)
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
