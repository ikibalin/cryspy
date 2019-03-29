import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
 
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot
 



def read_sc_data(finp):
    ddata={}
    fid = open(finp,'r')
    lcontentH = fid.readlines()
    fid.close()
    lname = [hh.strip() for hh in lcontentH[0].split()]
    for hh1 in lname:
        ddata[hh1] = []

    lcontent = [hh.strip() for hh in lcontentH[1:] if hh.strip() != ""]
    
    for line in lcontent:
        lhelp = [float(hh) for hh in line.split()]
        for hh1, hh2 in zip(lname, lhelp):
            ddata[hh1].append(hh2)
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
    def __init__(self, ffig_full = None):
        super(cwidg_central, self).__init__()
        self.init_layout_central(ffig_full)
        self.setLayout(self.layout_central)
        self.ffig_full = ffig_full
 
    def init_layout_central(self, ffig_full):
        
        lay_main = QtWidgets.QVBoxLayout()
        self.pcanvas  = PlotCanvas(self, width=5, height=4)
        if ffig_full != None:
            self.plot_file(ffig_full)
        
        lay_1 = QtWidgets.QVBoxLayout()
        lay_1.addWidget(self.pcanvas)
        
        """
        lay_vert_1 = QtWidgets.QHBoxLayout()
        b_open = QtWidgets.QPushButton('Open', self)
        b_open.clicked.connect(self.ask_open_file)
        lay_vert_1.addWidget(b_open)
        lay_vert_1.addStretch(1)
        """
        lay_main.addLayout(lay_1)
        """
        lay_main.addLayout(lay_vert_1)
        """
        
        self.layout_central = lay_main
    """
    def ask_open_file(self):
        if self.ffig_full == None:
            ffig_full, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd(),
                                                                  "Extension file (*.dat)")
            if (not ok):
                return
            self.ffig_full = ffig_full
        self.plot_file(self.ffig_full)
    """
    def plot_file(self, ffig_full):
        ddata = read_sc_data(ffig_full)
        

        lx = ddata["FRmod_ext"]
        ly = ddata["FRexp"]
        lsy = ddata["sFRexp"]


        
        self.pcanvas.plot(lx, ly, lsy)
        #self.pcanvas.set_limits(20,50,0,3000)            
        self.pcanvas.draw()
        
class PlotCanvas(matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg):
 
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi = dpi)
        #self.axes = fig.add_subplot(111)
 
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.__init__(self, fig)
        self.setParent(parent)
 
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding, 
                                   QtWidgets.QSizePolicy.Expanding)
        
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.updateGeometry(self)
        
        self.make_plot()
 
    def make_plot(self):
        self.ax1 = self.figure.add_subplot(111)
    """    
    def set_limits(self,xmin,xmax,ymin,ymax):
        self.ax1.set_xlim((xmin,xmax))
        self.ax1.set_ylim((ymin,ymax))
        self.ax2.set_xlim((xmin,xmax))
        self.ax2.set_ylim((ymin,ymax))
        self.ax3.set_xlim((xmin,xmax))
        self.ax3.set_ylim((-0.5*(ymax-ymin),0.5*(ymax-ymin)))
    """
        
    def plot(self, lx, ly, lsy):
        col_1 = "#000000"
        self.ax1.plot(lx, lx, "k-", linewidth=1.0)    
        self.ax1.errorbar(lx, ly, yerr = lsy, ecolor = col_1, fmt='o', color=col_1, linewidth = 0.5)
        
 
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = cwind_central()
    
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