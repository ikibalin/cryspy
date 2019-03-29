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
    def __init__(self, fname):
        super(cwind_central, self).__init__()
        self.title = "program 'Graph'"
        self.setWindowTitle(self.title)
        widg_central = cwidg_central(fname)
        self.setCentralWidget(widg_central)
        self.show()
        
class cwidg_central(QtWidgets.QWidget):
    def __init__(self, ffig_full = None):
        super(cwidg_central, self).__init__()
        self.init_layout_central(ffig_full)
        self.setLayout(self.layout_central)
        self.ffig_full = ffig_full
 
    def init_layout_central(self, ffig_full):
        
        lay_main = QtWidgets.QVBoxLayout()
        self.graph = Graph(self, width=5, height=4)
        

        lay_1 = QtWidgets.QVBoxLayout()
        lay_1.addWidget(self.graph)

        lay_main.addLayout(lay_1)

        
        self.layout_central = lay_main
        
        
        ddata = read_pwd_ddata(ffig_full)

        self.graph.data_1_x = ddata["ttheta"]
        self.graph.data_1_y_exp = ddata["IntUPexp"]
        self.graph.data_1_sy_exp = ddata["sIntUPexp"]
        self.graph.data_1_y_mod = ddata["IntUPmod"]
        self.graph.data_1_y_diff = [hh1 - hh2 for hh1, hh2 in zip(ddata["IntUPmod"], ddata["IntUPexp"])]
        

        
        self.graph.data_2_y_exp = ddata["IntDOWNexp"]
        self.graph.data_2_sy_exp = ddata["sIntDOWNexp"]
        self.graph.data_2_y_mod = ddata["IntDOWNmod"]
        self.graph.data_2_y_diff = [hh1 - hh2 for hh1, hh2 in zip(ddata["IntDOWNmod"], ddata["IntDOWNexp"])]


        self.graph.data_3_x = ddata["ttheta"]
        self.graph.data_3_y_exp = [hh1 - hh2 for hh1, hh2 in zip(ddata["IntUPexp"], ddata["IntDOWNexp"])]
        self.graph.data_3_sy_exp = [(hh1**2 + hh2**2)**0.5 for hh1, hh2 in zip(ddata["sIntUPexp"], ddata["sIntDOWNexp"])]
        self.graph.data_3_y_mod = [hh1 - hh2 for hh1, hh2 in zip(ddata["IntUPmod"], ddata["IntDOWNmod"])]
        self.graph.data_3_y_diff = [hh1 - hh2 for hh1, hh2 in zip(self.graph.data_3_y_mod, self.graph.data_3_y_exp)]
        
        #self.data_4_x =  ddata[""]    
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
        
        self.data_4_x =  []
        
        self.control = parent
        self.figure = fig
        
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.ax_pri = fig.add_subplot(111)
        
        fig.canvas.mpl_connect("button_press_event", self.onclick)
        
    def set_data_to_graph(self):

        
        col_1 = "#FF0000"
        col_2 = "#0000FF"
        col_3 = "#00FF00"
        lx = self.data_1_x
        ly1 = self.data_1_y_exp
        lsy1 = self.data_1_sy_exp
        y_up_min = min(ly1)
        ly1_mod = self.data_1_y_mod
        ly1_exp_mod = [-hh1+hh2+y_up_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax_pri.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_1, linewidth = 0.5)
        self.ax_pri.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax_pri.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)



        ly1_before = self.data_2_y_exp
        lsy1 = self.data_2_sy_exp
        y_down_max = max(ly1_before)
        y_down_min = min(ly1_before) - y_down_max + y_up_min
        ly1 = [hh - y_down_max + y_up_min  for hh in ly1_before]
        ly1_mod = [hh - y_down_max + y_up_min  for hh in self.data_2_y_mod]
        ly1_exp_mod = [-hh1+hh2 + y_down_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax_pri.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_2, linewidth = 0.5)
        self.ax_pri.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax_pri.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)


        ly1_before = self.data_3_y_exp
        y_diff_max = max(ly1_before)
        y_diff_min = min(ly1_before) - y_diff_max + y_down_min
        lsy1 = self.data_3_sy_exp
        ly1 = [hh - y_diff_max + y_down_min for hh in ly1_before]
        ly1_mod = [hh - y_diff_max + y_down_min  for hh in self.data_3_y_mod]
        ly1_exp_mod = [hh1-hh2 + y_diff_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax_pri.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_3, linewidth = 0.5)
        self.ax_pri.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax_pri.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)
        
        x_diff = max(lx)-min(lx)
        y_diff = max(self.data_1_y_exp) + y_diff_max - y_down_min
        self.xlim_orig = (min(lx)-0.05*x_diff, max(lx)+0.05*x_diff)
        self.ylim_orig = (- y_diff_max + y_down_min-0.05*y_diff, max(self.data_1_y_exp)+0.05*y_diff)
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
    ex = cwind_central(fname)
    
    sys.exit(app.exec_())
    
    
    
"""
    def clean_fig(self):
        self.ldata = []
        self.ldata_simple = []
        self.pcanvas.figure.clf()
        self.pcanvas.make_plot()
        self.pcanvas.draw()
        
    def replot_fig(self):
        lhelp=['#005555','#BB5555','#550055','#55BB55','#555500','#5555BB']
        for iddata, ddata in enumerate(self.ldata_simple):
            lflags = [ (isinstance(hh1,(float, int))&isinstance(hh2,(float, int))&isinstance(hh3,(float, int))&isinstance(hh4,(float, int))) for hh1,hh2,hh3,hh4 in zip(ddata["IntUP"],ddata["sIntUP"],ddata["IntDOWN"],ddata["sIntDOWN"])]
            lx = [hh for hh, flag in zip(ddata["ttheta"], lflags) if flag]
            li_u_exp = [hh for hh, flag in zip(ddata["IntUP"], lflags) if flag]
            lsi_u_exp = [hh for hh, flag in zip(ddata["sIntUP"], lflags) if flag]
            li_d_exp = [hh for hh, flag in zip(ddata["IntDOWN"], lflags) if flag]
            lsi_d_exp = [hh for hh, flag in zip(ddata["sIntDOWN"], lflags) if flag]
            
            #li_u_exp = ddata["IntUP"]
            #lsi_u_exp = ddata["sIntUP"]
            #li_d_exp = ddata["IntDOWN"]
            #lsi_d_exp = ddata["sIntDOWN"]
            ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
            lsdiff_exp = [(hh1**2+hh2**2)**0.5  for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
            lcolors=[lhelp[iddata] for hh in range(3)]
            self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,None],[li_d_exp,lsi_d_exp,None],[ldiff_exp,lsdiff_exp,None],None,lcolors)
        xmin = float(self.le_xmin.text())
        xmax = float(self.le_xmax.text())
        ymin = float(self.le_ymin.text())
        ymax = float(self.le_ymax.text())
        self.pcanvas.set_limits(xmin,xmax,ymin,ymax)          
        self.pcanvas.draw()
"""    