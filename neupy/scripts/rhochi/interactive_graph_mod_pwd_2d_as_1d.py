import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
import numpy
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot
 

def read_data(file_name):

    fid = open(file_name,'r')
    l_cont = fid.readlines()
    fid.close()

    lmat, mat = [], []
    for hh in l_cont:
        if hh.strip() == "":
            if mat != []:
                lmat.append(mat)
                mat = []
        else:
            mat.append(hh)
    if mat != []:
        lmat.append(mat)
        mat = []
    
    ll_int, l_ang1, l_ang2 = [], [], []
    for mat in lmat:
        ang1 = [float(hh) for hh in mat[0].split()[1:]]
        l_int, ang2 = [], []
        for hh in mat[1:]:
            lhelp = hh.split()
            ang2.append(float(lhelp[0]))
            l_int.append([float(hh2) if hh2 != "None" else None for hh2 in lhelp[1:]])
        ll_int.append(l_int)
        l_ang1.append(ang1)
        l_ang2.append(ang2)

    tth = numpy.array(l_ang1[0], dtype=float)
    phi = numpy.array(l_ang2[0], dtype=float)
    int_s_e = numpy.array(ll_int[0], dtype=float).transpose()
    int_s_m = numpy.array(ll_int[1], dtype=float).transpose()
    sint_s = numpy.array(ll_int[2], dtype=float).transpose()
    int_d_e = numpy.array(ll_int[3], dtype=float).transpose()
    int_d_m = numpy.array(ll_int[4], dtype=float).transpose()
    sint_d = numpy.array(ll_int[5], dtype=float).transpose()
    return tth, phi, int_s_e, int_s_m, sint_s, int_d_e, int_d_m, sint_d


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
        tth, phi, int_s_e, int_s_m, sint_s, int_d_e, int_d_m, sint_d = read_data(ffig_full)

        i_s_e_1d = numpy.where(numpy.isnan(int_s_e), 0., int_s_e).sum(axis=1)
        i_s_m_1d = numpy.where(numpy.isnan(int_s_e), 0., int_s_m).sum(axis=1)
        si_s_e_1d = ((numpy.where(numpy.isnan(int_s_e), 0., sint_s)**2).sum(axis=1))**0.5
        i_d_e_1d = numpy.where(numpy.isnan(int_d_e), 0., int_d_e).sum(axis=1)
        i_d_m_1d = numpy.where(numpy.isnan(int_d_e), 0., int_d_m).sum(axis=1)
        si_d_e_1d = ((numpy.where(numpy.isnan(int_d_e), 0., sint_d)**2).sum(axis=1))**0.5

        self.graph.data_1_x = []
        self.graph.data_1_y_exp = []
        self.graph.data_1_sy_exp = []
        self.graph.data_1_y_mod = []
        self.graph.data_1_y_diff = []
        
        
        self.graph.data_2_y_exp = []
        self.graph.data_2_sy_exp = []
        self.graph.data_2_y_mod = []
        self.graph.data_2_y_diff = []


        self.graph.data_3_x = tth
        self.graph.data_3_y_exp = i_s_e_1d
        self.graph.data_3_sy_exp = si_s_e_1d
        self.graph.data_3_y_mod = i_s_m_1d
        self.graph.data_3_y_diff = i_s_e_1d-i_s_m_1d
        

        self.graph.data_4_x = tth
        self.graph.data_4_y_exp = i_d_e_1d
        self.graph.data_4_sy_exp = si_d_e_1d
        self.graph.data_4_y_mod = i_d_m_1d
        self.graph.data_4_y_diff = i_d_e_1d-i_d_m_1d

        
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


        #ly1_before = self.data_4_y_exp
        #y_diff_max = max(ly1_before)
        #y_diff_min = min(ly1_before) - y_diff_max + y_down_min

        x_diff = max(lx)-min(lx)
        y_diff = max(self.data_4_y_exp) - y_down_min
        self.xlim_orig = (min(lx)-0.05*x_diff, max(lx)+0.05*x_diff)
        self.ylim_orig = (y_down_min-0.05*y_diff, max(self.data_3_y_exp)+0.05*y_diff)
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
