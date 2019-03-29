import os
import sys

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore

import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure
import matplotlib.pyplot



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
                                                                  "Extension file (*.out)")
            if (not ok):
                return
            self.ffig_full = ffig_full
        self.plot_file(self.ffig_full)
    """
    def plot_file(self, ffig_full):
        cdata = powder_data(ffig_full)


        lhelp=['#005555','#BB5555','#550055','#55BB55','#555500','#5555BB']

        lflag = [(hh == 1) for hh in cdata.data["flag"]]
        lx = [hh for hh, flag in zip(cdata.data["ttheta"], lflag) if flag]
        li_u_exp = [hh for hh, flag in zip(cdata.data["IntUPexp"], lflag) if flag]
        lsi_u_exp = [hh for hh, flag in zip(cdata.data["sIntUPexp"], lflag) if flag]
        li_d_exp = [hh for hh, flag in zip(cdata.data["IntDOWNexp"], lflag) if flag]
        lsi_d_exp = [hh for hh, flag in zip(cdata.data["sIntDOWNexp"], lflag) if flag]
        li_u_mod = [hh for hh, flag in zip(cdata.data["IntUPmod"], lflag) if flag]
        li_d_mod = [hh for hh, flag in zip(cdata.data["IntDOWNmod"], lflag) if flag]
        li_b_mod = [hh for hh, flag in zip(cdata.data["IntBKGR"], lflag) if flag]
        ltth_hkl = cdata.phases[0]['ttheta']
        ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
        lsdiff_exp = [(hh1**2+hh2**2)**0.5 for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
        ldiff_mod =  [hh1-hh2 for hh1, hh2 in zip(li_u_mod, li_d_mod)]
        lcolors = lhelp[:3]
        self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,li_u_mod],[li_d_exp,lsi_d_exp,li_d_mod],[ldiff_exp,lsdiff_exp,ldiff_mod],ltth_hkl,lcolors)




        #li_u_exp = ddata["IntUP"]
        #lsi_u_exp = ddata["sIntUP"]
        #li_d_exp = ddata["IntDOWN"]
        #lsi_d_exp = ddata["sIntDOWN"]
        ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
        lsdiff_exp = [(hh1**2+hh2**2)**0.5  for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
        lcolors=lhelp[:3]

        self.pcanvas.plot(lx,[li_u_exp, lsi_u_exp, li_u_mod],[li_d_exp, lsi_d_exp, li_d_mod],[ldiff_exp,lsdiff_exp, ldiff_mod], ltth_hkl, lcolors)
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

    def plot(self, lx, lly_u, lly_d, lly_diff, ltth_hkl = None,lcolors = None):
        if lcolors!=None:
            col_1 = lcolors[0]
            col_2 = lcolors[1]
            col_3 = lcolors[2]
        else:
            col_1 = "#FF0000"
            col_2 = "#0000FF"
            col_3 = "#00FF00"

        ly1 = lly_u[0]
        lsy1 = lly_u[1]
        y_up_min = min(ly1)
        ly1_mod = lly_u[2]
        ly1_exp_mod = [hh1-hh2+y_up_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax1.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_1, linewidth = 0.5)
        self.ax1.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax1.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)



        ly1_before = lly_d[0]
        lsy1 = lly_d[1]
        y_down_max = max(ly1_before)
        y_down_min = min(ly1_before) - y_down_max + y_up_min
        ly1 = [hh - y_down_max + y_up_min  for hh in ly1_before]
        ly1_mod = [hh - y_down_max + y_up_min  for hh in lly_d[2]]
        ly1_exp_mod = [hh1-hh2 + y_down_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax1.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_2, linewidth = 0.5)
        self.ax1.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax1.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)


        ly1_before = lly_diff[0]
        y_diff_max = max(ly1_before)
        y_diff_min = min(ly1_before) - y_diff_max + y_down_min
        lsy1 = lly_diff[1]
        ly1 = [hh - y_diff_max + y_down_min for hh in ly1_before]
        ly1_mod = [hh - y_diff_max + y_down_min  for hh in lly_diff[2]]
        ly1_exp_mod = [hh1-hh2 + y_diff_min for hh1, hh2 in zip(ly1, ly1_mod)]

        self.ax1.errorbar(lx, ly1, yerr = lsy1, ecolor = col_1, linestyle='-', color=col_3, linewidth = 0.5)
        self.ax1.plot(lx, ly1_mod, "k-",linewidth=1.0)
        self.ax1.plot(lx, ly1_exp_mod, "k-",linewidth=0.5)

        if (ltth_hkl != None):
            self.ax1.plot(ltth_hkl, [0. for hh in ltth_hkl], "k|")


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