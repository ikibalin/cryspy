import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
 
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot
 

 
import random

def splitlinewminuses(line):
    """Split line even if between numbers there is no spaces and only minus.\n
    If it is possible transform numbers to float or integer it will be done.\n
    """
    lbefore1,lbefore2='',''
    n1=''
    splitedline=[]
    for letter in line:
        if letter.isdigit():
            n1+=letter
        elif ((letter=='.') and (lbefore1.isdigit())):
            n1+=letter
        elif ((letter=='-') and (lbefore1=='E') and (lbefore2.isdigit())):
            n1+=letter
        elif ((letter=='-') and (lbefore1!='E')):
            if ((n1!=' ') and (n1!='')):
                splitedline.append(n1)
            n1=letter
        elif ((letter==' ') and (not((n1==' ') or (n1=='')))):
            splitedline.append(n1)
            n1=''
        else:
            n1+=letter
        lbefore2=lbefore1
        lbefore1=letter
        if (n1==' '):
            n1=''
    if ((n1!=' ') and (n1!='')):
        splitedline.append(n1)
        n1=''
    fsplitedline=[]
    for elsplitedline in splitedline:
        try:
            if elsplitedline[0]=='-':
                hhh=elsplitedline[1:]
            else:
                hhh=elsplitedline
            if hhh.isdigit():
                hh=int(elsplitedline)
            else:
                hh=float(elsplitedline)
            fsplitedline.append(hh)
        except:
            fsplitedline.append(elsplitedline)
    if (len(fsplitedline)==1):
        fsplitedline=fsplitedline[0]
    return fsplitedline


def strings_to_ddata(lstring):
    lcontent = [hh1 for hh1 in lstring if hh1[0] != "#"]
    lname=lcontent[0].strip().split()
    ddata = {}
    for name in lname:
        ddata[name]=[]
    for line in lcontent[1:]:
        if (line=='\n'): break
        for name,sval in zip(lname, line.strip().split()):
            ddata[name].append(float(sval))
    return ddata

def readexpdata(finp):
    ddata={}
    fid=open(finp,'r')
    lcontentH=fid.readlines()
    fid.close()
    lparam=[line[1:] for line in lcontentH if line[0]=='#']
    if (len(lparam)>1):
        for line in lparam[:-1]:
            lhelp=splitlinewminuses(line)
            if (len(lhelp)>2):
                ddata[lhelp[0]]=lhelp[1:]
            elif (len(lhelp)==2):
                ddata[lhelp[0]]=lhelp[1]
            else:
                print "Mistake in experimental file '{}' in line:\n{}".format(finp,line)
                print "The program is stopped."
                quit()
    lnames=lparam[-1].split()
    for name in lnames:
        ddata[name]=[]
    lcontent=[line for line in lcontentH if line[0]!='#']
    for line in lcontent:
        for name,val in zip(lnames,splitlinewminuses(line)):
            ddata[name].append(val)
    return ddata

class powder_data(object):
    def __init__(self, fname = None):
        self.phases = []
        self.data= {}
        if fname != None:
            self.load_data(fname)
    def load_data(self, fname):
        fid = open(fname,'r')
        lcontentH = fid.readlines()
        fid.close()
        
        lnumb = [ihh1 for ihh1, hh1 in enumerate(lcontentH) if hh1[0] == "\n"]
        lcontent1=lcontentH[0:lnumb[0]]
        lnumb_2 = [ hh1 for ihh1, hh1 in enumerate(lnumb[:-1]) if (hh1+1) != lnumb[ihh1+1]]
        if (len(lcontentH)-1) != lnumb_2[-1]:
            lnumb_2.append(len(lcontentH)-1)


        ddata = strings_to_ddata(lcontent1)
        self.data = ddata
        
        for inum1,num1 in enumerate(lnumb_2[:-1]):
            lcontent1 = lcontentH[num1+1:lnumb_2[inum1+1]]
            dphase = strings_to_ddata(lcontent1)
            self.phases.append(dphase)

 
class App(QtWidgets.QMainWindow):
 
    def __init__(self, fdir_file = None):
        super(App, self).__init__()
        self.title = "program 'Graph'"
        self.initUI()
        self.ldata = []
        self.ldata_simple = []
        if fdir_file == None:
            fdir_file = os.getcwd()
        self.fdir_file = fdir_file
 
    def initUI(self):
        self.setWindowTitle(self.title)
        #self.setGeometry(self.left, self.top, self.width, self.height)
        
        lay_main = QtWidgets.QVBoxLayout()
        canv_up = PlotCanvas(self, width=5, height=4)
        canv_down = PlotCanvas(self, width=5, height=4)
        canv_diff = PlotCanvas(self, width=5, height=4)
        canv_all = PlotCanvas(self, width=5, height=4)
        self.canv_up = canv_up
        self.canv_down = canv_down
        self.canv_diff = canv_diff
        self.canv_all = canv_all
 
        lay_2 = QtWidgets.QVBoxLayout()
        lay_2.addWidget(canv_up)
        lay_2.addWidget(canv_down)
        lay_2.addWidget(canv_diff)
    
        lay_1 = QtWidgets.QHBoxLayout()
        lay_1.addWidget(canv_all)
        lay_1.addLayout(lay_2)
        
        
    
        lay_vert_1 = QtWidgets.QHBoxLayout()
        b_open = QtWidgets.QPushButton('Open', self)
        b_open.clicked.connect(self.open_file)
        lay_vert_1.addWidget(b_open)

        b_open_s = QtWidgets.QPushButton('Open (just data)', self)
        b_open_s.clicked.connect(self.open_file_simple)
        lay_vert_1.addWidget(b_open_s)

        b_clean = QtWidgets.QPushButton('Clean', self)
        b_clean.clicked.connect(self.clean_fig)
        lay_vert_1.addWidget(b_clean)
        
        b_replot = QtWidgets.QPushButton('Replot', self)
        b_replot.clicked.connect(self.replot_fig)
        lay_vert_1.addWidget(b_replot)
        lay_vert_1.addStretch(1)

        lay_vert_2 = QtWidgets.QHBoxLayout()
        l_xmin = QtWidgets.QLabel("x_min:")
        self.le_xmin = QtWidgets.QLineEdit("")
        l_xmax = QtWidgets.QLabel("x_max:")
        self.le_xmax = QtWidgets.QLineEdit("")
        l_ymin = QtWidgets.QLabel("y_min:")
        self.le_ymin = QtWidgets.QLineEdit("")
        l_ymax = QtWidgets.QLabel("y_max:")
        self.le_ymax = QtWidgets.QLineEdit("")
        lay_vert_2.addWidget(l_xmin)
        lay_vert_2.addWidget(self.le_xmin)
        lay_vert_2.addWidget(l_xmax)
        lay_vert_2.addWidget(self.le_xmax)
        lay_vert_2.addWidget(l_ymin)
        lay_vert_2.addWidget(self.le_ymin)
        lay_vert_2.addWidget(l_ymax)
        lay_vert_2.addWidget(self.le_ymax)
        lay_vert_2.addStretch(1)
        
        lay_main.addLayout(lay_1)
        lay_main.addLayout(lay_vert_2)
        lay_main.addLayout(lay_vert_1)
        
        widg_main = QtWidgets.QWidget()
        widg_main.setLayout(lay_main)
        
        
        self.setCentralWidget(widg_main)
 
        self.show()

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
        for icdata, cdata in enumerate(self.ldata):
            lx = cdata.data["ttheta"]
            li_u_exp = cdata.data["IntUPexp"]
            lsi_u_exp = cdata.data["sIntUPexp"]
            li_d_exp = cdata.data["IntDOWNexp"]
            lsi_d_exp = cdata.data["sIntDOWNexp"]
            li_u_mod = cdata.data["IntUPmod"]
            li_d_mod = cdata.data["IntDOWNmod"]
            li_b_mod = cdata.data["IntBKGR"]
            ltth_hkl = cdata.phases[0]['ttheta']
            ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
            lsdiff_exp = [(hh1**2+hh2**2)**0.5 for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
            ldiff_mod =  [hh1-hh2 for hh1, hh2 in zip(li_u_mod, li_d_mod)]
            lcolors=[lhelp[icdata] for hh in range(3)]
            self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,li_u_mod],[li_d_exp,lsi_d_exp,li_d_mod],[ldiff_exp,lsdiff_exp,ldiff_mod],ltth_hkl,lcolors)
        xmin = float(self.le_xmin.text())
        xmax = float(self.le_xmax.text())
        ymin = float(self.le_ymin.text())
        ymax = float(self.le_ymax.text())
        self.pcanvas.set_limits(xmin,xmax,ymin,ymax)
        #self.pcanvas.set_limits(20,50,0,3000)            
        self.pcanvas.draw()
        
    def open_file_simple(self):
        fdir_file = self.fdir_file
        lhelp = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdir_file,
                                                      "Extension file (*.dat; *.out)")
        fdir_file = os.path.dirname("{}".format(lhelp[0]))
        fdir_file = fdir_file.replace("/","\\")
        fname_file = os.path.basename("{}".format(lhelp[0]))
        self.fdir_file = fdir_file 
        ddata = readexpdata(os.path.join(fdir_file, fname_file))
        self.ldata_simple.append(ddata)
        lhelp=['#005555','#BB5555','#550055','#55BB55','#555500','#5555BB']
        for icdata, ddata in enumerate(self.ldata_simple):
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
            lcolors=[lhelp[icdata] for hh in range(3)]
            self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,None],[li_d_exp,lsi_d_exp,None],[ldiff_exp,lsdiff_exp,None],None,lcolors)
        #self.pcanvas.set_limits(20,50,0,3000)            
        self.pcanvas.draw()
        
    def open_file(self):
        fdir_file = self.fdir_file
        lhelp = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdir_file,
                                                      "Extension file (*.out)")
        fdir_file = os.path.dirname("{}".format(lhelp[0]))
        fdir_file = fdir_file.replace("/","\\")
        fname_file = os.path.basename("{}".format(lhelp[0]))
        self.fdir_file = fdir_file 
        #fdir_file = os.getcwd()
        #fname_file = "oFe3O4p1_10K_5T.out"
        cdata = powder_data(os.path.join(fdir_file, fname_file))
        self.ldata.append(cdata)
        lhelp=['#005555','#BB5555','#550055','#55BB55','#555500','#5555BB']
        for icdata, cdata in enumerate(self.ldata):
            lx = cdata.data["ttheta"]
            li_u_exp = cdata.data["IntUPexp"]
            lsi_u_exp = cdata.data["sIntUPexp"]
            li_d_exp = cdata.data["IntDOWNexp"]
            lsi_d_exp = cdata.data["sIntDOWNexp"]
            li_u_mod = cdata.data["IntUPmod"]
            li_d_mod = cdata.data["IntDOWNmod"]
            li_b_mod = cdata.data["IntBKGR"]
            ltth_hkl = cdata.phases[0]['ttheta']
            ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
            lsdiff_exp = [(hh1**2+hh2**2)**0.5 for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
            ldiff_mod =  [hh1-hh2 for hh1, hh2 in zip(li_u_mod, li_d_mod)]
            lcolors=[lhelp[icdata] for hh in range(3)]
            self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,li_u_mod],[li_d_exp,lsi_d_exp,li_d_mod],[ldiff_exp,lsdiff_exp,ldiff_mod],ltth_hkl,lcolors)
        self.pcanvas.set_limits(20,50,0,3000)            
        self.pcanvas.draw()
        
class PlotCanvas(matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg):
 
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = matplotlib.figure.Figure(figsize=(width, height), dpi=dpi)
        #self.axes = fig.add_subplot(111)
 
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.__init__(self, fig)
        self.setParent(parent)
 
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding, 
                                   QtWidgets.QSizePolicy.Expanding)
        
        matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg.updateGeometry(self)
        
        self.make_plot()
        #self.plot(data)
 
    def make_plot(self):
        self.ax1 = self.figure.add_subplot(111)
        
    def set_limits(self,xmin,xmax,ymin,ymax):
        self.ax1.set_xlim((xmin,xmax))
        self.ax1.set_ylim((ymin,ymax))
        self.ax2.set_xlim((xmin,xmax))
        self.ax2.set_ylim((ymin,ymax))
        self.ax3.set_xlim((xmin,xmax))
        self.ax3.set_ylim((-0.5*(ymax-ymin),0.5*(ymax-ymin)))

        
    def plot(self, lx, lly_u, lly_d, lly_diff, ltth_hkl = None,lcolors = None):
        ly1 = lly_u[0]
        lsy1 = lly_u[1]
        ly2 = lly_u[2]
        if lcolors!=None:
            col_1 = lcolors[0]
            col_2 = lcolors[1]
            col_3 = lcolors[2]
        else:
            col_1 = "#FF0000"
            col_2 = "#0000FF"
            col_3 = "#00FF00"
        #self.ax1.plot(lx, ly1, marker = ".", markersize=2.0, color = col_1, linestyle='')
        if (ly2 != None):
            self.ax1.plot(lx, ly2, "k-")
            self.ax1.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='')
            ly3 = [hh1-hh2 for hh1,hh2 in zip(ly1,ly2)]
            self.ax1.plot(lx, ly3, "k-",linewidth=0.5)
        else:
            self.ax1.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_1, linewidth = 0.5)
        if (ltth_hkl != None):
            self.ax1.plot(ltth_hkl, [0. for hh in ltth_hkl], "k|")
        
        ly1 = lly_d[0]
        lsy1 = lly_d[1]
        ly2 = lly_d[2]
        #self.ax2.plot(lx, ly1, marker = ",", color = col_2)
        if (ly2 != None):
            self.ax2.plot(lx, ly2, "k-")
            self.ax2.errorbar(lx, ly1, yerr = lsy1, ecolor = col_2, linestyle='')
            ly3 = [hh1-hh2 for hh1,hh2 in zip(ly1,ly2)]
            self.ax2.plot(lx, ly3, "k-",linewidth=0.5)
        else:
            self.ax2.errorbar(lx, ly1, yerr = lsy1, ecolor = col_2, linestyle='-', color=col_2, linewidth = 0.5)
        if (ltth_hkl != None):
            self.ax2.plot(ltth_hkl, [0. for hh in ltth_hkl], "k|")
        
        ly1 = lly_diff[0]
        lsy1 = lly_diff[1]
        ly2 = lly_diff[2]
        #self.ax3.plot(lx, ly1, marker = ",", color = col_3)
        if (ly2 != None):
            self.ax3.plot(lx, ly2, "k-")
            self.ax3.errorbar(lx, ly1, yerr = lsy1, ecolor = col_3, linestyle='')
            ly3 = [hh1-hh2 for hh1,hh2 in zip(ly1,ly2)]
            self.ax3.plot(lx, ly3, "k-",linewidth=0.5)
        else:
            self.ax3.errorbar(lx, ly1, yerr = lsy1, ecolor = col_3, linestyle='-', color=col_3, linewidth = 0.5)
            
        if (ltth_hkl != None):
            self.ax3.plot(ltth_hkl, [0. for hh in ltth_hkl], "k|")
        #ax.set_title('')
 
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    ex = App()
    
    sys.exit(app.exec_())
    