# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 09:52:24 2019

Создает главное окно, определяет компоновщик, вызывает первую модель, давая ее представлениям

@author: ikibalin
"""
import os
import sys

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore

import f_rhochi_viewer.interactive_graph_mod_pwd

# =============================================================================
# #rewrite it
# =============================================================================

def take_f_name_from_rcif(f_rcif, s_search="_pd_file_name_output"):
    fid = open(f_rcif, 'r')
    l_cont = fid.readlines()
    fid.close
    l_answ = []
    for line in l_cont:
        hh= line.strip()
        if hh.startswith(s_search):
            l_help = hh.split()
            if len(l_help)>=2:
                l_answ.append(l_help[1])
    f_dir = os.path.dirname(f_rcif)
    l_res = []
    for answ in l_answ:
        f_name = os.path.join(f_dir, answ)
        if os.path.isfile(f_name):
            l_res.append(f_name)
    return l_res


    


def del_layout(layout):
    """
    delete all elements from layouts
    """
    for i in reversed(range(layout.count())):
        if layout.itemAt(i).widget() != None:
           layout.itemAt(i).widget().setParent(None)
        elif layout.itemAt(i).layout() != None:
            del_layout(layout.itemAt(i).layout())
            layout.itemAt(i).layout().setParent(None)
        else:
            layout.removeItem(layout.itemAt(i))
    return



class mythread(QtCore.QThread):
    def __init__(self, parent=None):
        QtCore.QThread.__init__(self,parent)
        self.signal = None
        self.core = None

        #self.startpatam = startparam
    def run(self):
        ccore = self.core
        try:
            lmessage = ccore.run_refinement()
            model = cmodel.cmodel()
            model.take_from_core(ccore)
            self.signal.model = model
        except:
            lmessage = ["Undefined mistake during the refinement procedure.\nSomething wrong."]
        print("Calculations are stopped")
        self.signal.message = lmessage
        self.signal.rename_signal.emit()



class cthread_signal(QtCore.QObject):
    rename_signal = QtCore.pyqtSignal()



class cbuilder(QtWidgets.QMainWindow):
    def __init__(self, f_dir_prog = os.path.dirname(__file__)):
        super(cbuilder, self).__init__()
        self._p_f_dir_prog = f_dir_prog 
        self._p_d_setup = {}

        self.read_setup()
        self.def_actions()
        self.init_widget()

        #self.setStyleSheet("background-color:lightgray;")
        self.show()

    def read_setup(self):
        f_dir_prog = self._p_f_dir_prog
        f_file_prog_setup = os.path.join(f_dir_prog, 'f_rhochi_viewer', 'setup.dat')

        if not(os.path.isfile(f_file_prog_setup)):
            self.write_setup()

        fid = open(f_file_prog_setup, "r")
        l_cont = fid.readlines()
        fid.close()
        d_setup = {}
        for line in l_cont:
            l_help = [hh.strip() for hh in line.strip().split("::")]
            if len(l_help) == 2:
                d_setup.update({l_help[0]: l_help[1]})
        
        self._p_d_setup = d_setup
        
        if "f_dir_data" in d_setup.keys():
            self.setWindowTitle('RhoChi Viewer - {:}'.format(d_setup["f_dir_data"]))

    def write_setup(self):
        l_key = ["f_dir_data"]
        l_key_val = [self._p_f_dir_prog]

        d_setup = self._p_d_setup

        if (("f_name_data" in d_setup.keys()) & ("f_dir_data" in d_setup.keys())):
            self.setWindowTitle('RhoChi Viewer - {:}'.format(d_setup["f_dir_data"]))


        l_key_setup = d_setup.keys()
        for key, key_val in zip(l_key, l_key_val):
            if key not in l_key_setup:
                d_setup.update({key: key_val})
        ls_out = ["{:}::{:}".format(key, val) for key, val in d_setup.items()]

        f_dir_prog = self._p_f_dir_prog
        f_file_prog_setup = os.path.join(f_dir_prog, 'f_rhochi_viewer', 'setup.dat')
        fid = open(f_file_prog_setup, "w")
        fid.write("\n".join(ls_out))
        fid.close()


    def def_actions(self):

        f_dir_prog = self._p_f_dir_prog
        f_dir_prog_icon = os.path.join(f_dir_prog, 'f_rhochi_viewer', 'f_icon')

        open_rcif_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'open.png')),'&Open data', self)
        open_rcif_action.setShortcut('Ctrl+O')
        open_rcif_action.setStatusTip('Open data (1D)')
        open_rcif_action.triggered.connect(self.open_rcif)


        run_rhochi = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'calc.png')),'&Run RhoChi', self)
        run_rhochi.setShortcut('Ctrl+R')
        run_rhochi.setStatusTip('Model optimization of marked parameters')
        run_rhochi.triggered.connect(self.run_rhochi)

        read_rcif_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'read_rcif.png')),'&Read .rcif', self)
        read_rcif_action.setStatusTip('Read .rcif')
        read_rcif_action.triggered.connect(self.read_rcif)

        read_bkgr_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'read_bkgr.png')),'&Read .bkgr', self)
        read_bkgr_action.setStatusTip('Read background')
        read_bkgr_action.triggered.connect(self.read_bkgr)


        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(open_rcif_action)

        calcMenu = menubar.addMenu('&Calculations')
        calcMenu.addAction(run_rhochi)
        calcMenu.addAction(read_rcif_action)


        self.toolbar = self.addToolBar("Open")
        self.toolbar.addAction(open_rcif_action)
        self.toolbar.addAction(run_rhochi)
        self.toolbar.addAction(read_rcif_action)
        self.toolbar.addAction(read_bkgr_action)

    def init_widget(self): 
        """
        module is reloaded for cwidget
        """
        self.location_on_the_screen()
        self._p_widg_central = cwidget(self.info_width)

        self.setCentralWidget(self._p_widg_central)
        self.try_to_plot_from_rcif()

    def try_to_plot_from_rcif(self):
        try:    
            f_name_data = self._p_d_setup["f_name_data"]
            f_dir_data = self._p_d_setup["f_dir_data"]
            f_name_full = os.path.join(f_dir_data, f_name_data)
            l_f_name_data = take_f_name_from_rcif(f_name_full, "_pd_file_name_output")
            if len(l_f_name_data) >= 1:
                self.plot_data(l_f_name_data[0])
        except:
            pass

    def location_on_the_screen(self):
        """
        position and size of main window
        """
        screen = QtWidgets.QDesktopWidget().screenGeometry()
        self.setMinimumSize(screen.width()*1/4, screen.height()*1/4)
        #self.setMaximumSize(screen.width()*3/4, screen.height()*3/4)
        self.info_width = screen.width()*3/4
        self.move(screen.width()/8 , screen.height()/16)

    def open_rcif(self):
        self.open_file(s_constr="Rcif files (*.rcif)")
        self.ask_file(ext=".rcif")
        self.try_to_plot_from_rcif()


    def open_file(self, s_constr="All files (*.*)"):
        #dir_ = QtWidgets.QFileDialog.getOpenFileName(self, 'Select a folder:', 'C:\\', "XML files (*.xml);;HTML files (*.html);;"
        #        "SVG files (*.svg);;User Interface files (*.ui)")
        if "f_dir_data" in self._p_d_setup.keys():
            f_dir_data = self._p_d_setup["f_dir_data"]
            if not(os.path.isdir(f_dir_data)):
                f_dir_data = self._p_f_dir_prog
        else:
            f_dir_data = self._p_f_dir_prog
        f_file_data_new, res = QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file:', f_dir_data, s_constr)

        if f_file_data_new != "":
            f_dir_data = os.path.dirname(f_file_data_new)
            f_name_data = os.path.basename(f_file_data_new)
            self._p_d_setup["f_dir_data"] = f_dir_data
            self._p_d_setup["f_name_data"] = f_name_data
            self.write_setup()

    def ask_file(self, ext=".rcif"):
        cond_ask = False
        cond_cycle = True
        while cond_cycle:
            l_key = self._p_d_setup.keys()
            if (("f_name_data" in l_key) & ("f_dir_data" in l_key)):
                f_name_data = self._p_d_setup["f_name_data"]
                f_dir_data = self._p_d_setup["f_dir_data"]
                f_name_full = os.path.join(f_dir_data, f_name_data)
                cond_1 = f_name_data.endswith(ext)
                cond_2 = os.path.isfile(f_name_full)
                cond_ask = not(cond_1 & cond_2)
            else:
                cond_ask = True

            if cond_ask:
                self.open_file()
            else:
                cond_cycle = False

    def run_rhochi(self):
        self.ask_file(".rcif")
        f_name_data = self._p_d_setup["f_name_data"]
        f_dir_data = self._p_d_setup["f_dir_data"]
        f_name_full = os.path.join(f_dir_data, f_name_data)

        f_dir_prog = self._p_f_dir_prog
        f_prog_full = os.path.join(f_dir_prog, "rhochi.py")
        line = "python {:} {:} {:}".format(f_prog_full, f_name_full, f_name_full)
        os.system(line)
        self.try_to_plot_from_rcif()

            

    def read_rcif(self):
        self.ask_file(".rcif")
        f_name_data = self._p_d_setup["f_name_data"]
        f_dir_data = self._p_d_setup["f_dir_data"]
        f_name_full = os.path.join(f_dir_data, f_name_data)
        os.startfile(f_name_full)

    def read_bkgr(self):
        self.ask_file(".rcif")
        f_name_data = self._p_d_setup["f_name_data"]
        f_dir_data = self._p_d_setup["f_dir_data"]
        f_rcif_full = os.path.join(f_dir_data, f_name_data)
        l_f_bkgr = take_f_name_from_rcif(f_rcif_full, "_pd_file_name_bkgr")
        if len(l_f_bkgr) >= 1:
            os.startfile(l_f_bkgr[0])

    def plot_data(self, f_name_full):
        widg_central = self._p_widg_central
        widg_graph = widg_central._p_widg_graph
        widg_graph.plot_file(f_name_full)


class cwidget(QtWidgets.QWidget):
    def __init__(self, width):
        self.info_width = width
        super(cwidget, self).__init__()
        self.init_widget()
        self.setLayout(self.layout_central)

    def init_widget(self):
        """
        make central layout
        """
        #lay_central = self.layout()
        #if lay_central == None:
        #    lay_central = QtWidgets.QHBoxLayout()
        #else:
        #    del_layout(lay_central)
        lay_central = QtWidgets.QHBoxLayout()


        width_m_1 = self.info_width/10
        width_m_2 = (4*self.info_width)/10
        width_m_3 = self.info_width-width_m_1-width_m_2-20

        self.width_cpanel = width_m_1
        self.width_left = width_m_2
        self.width_right = width_m_3


        self.lay_cpanel = QtWidgets.QVBoxLayout()
        self.lay_left = QtWidgets.QVBoxLayout()
        self.lay_right = QtWidgets.QVBoxLayout()


        lay_central.addLayout(self.lay_cpanel)
        lay_central.addLayout(self.lay_left)
        lay_central.addLayout(self.lay_right)

        self.layout_central = lay_central

        self.form_widg_cpanel()
        self.lay_left.addWidget(QtWidgets.QFrame())
        self.form_widg_right()
        


    def form_widg_cpanel(self):


        #del_layout(self.lay_cpanel)

        width_m_1 = self.width_cpanel


        #self.widg_ph = widg_cb_ph.cwidget(model)
        #self.widg_exp = widg_cb_exp.cwidget(model)
        

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)

        #b_spgr = QtWidgets.QPushButton("unit cell")
        #b_spgr.setMaximumSize(width_m_1, 50)
        #b_spgr.setSizePolicy(sizePolicy)
        ##b_spgr.clicked.connect(lambda: self.form_widg_left(0, model))
        #b_nucl = QtWidgets.QPushButton("nucl. structure")
        ##b_nucl.clicked.connect(lambda: self.form_widg_left(1, model))
        #b_adp = QtWidgets.QPushButton("anis.displacement")
        ##b_adp.clicked.connect(lambda: self.form_widg_left(2, model))
        #b_ms = QtWidgets.QPushButton("magnetic structure")
        ##b_ms.clicked.connect(lambda: self.form_widg_left(3, model))
        #b_ms_i = QtWidgets.QPushButton("Ising model")
        ##b_ms_i.clicked.connect(lambda: self.form_widg_left(4, model))
        #b_exp = QtWidgets.QPushButton("experiment")
        ##b_exp.clicked.connect(lambda: self.form_widg_left(5, model))
        #b_profile = QtWidgets.QPushButton("profile")
        ##b_profile.clicked.connect(lambda: self.form_widg_left(6, model))
        #b_scale = QtWidgets.QPushButton("scale")
        ##b_scale.clicked.connect(lambda: self.form_widg_left(7, model))
        #b_extinction = QtWidgets.QPushButton("extinction")
        ##b_extinction.clicked.connect(lambda: self.form_widg_left(8, model))
        #b_output = QtWidgets.QPushButton("output")
        ##b_output.clicked.connect(lambda: self.form_widg_left(9, model))
        #b_sconstr = QtWidgets.QPushButton("set constraints")
        ##b_sconstr.clicked.connect(lambda: self.form_widg_left(10, model))
        #b_hide = QtWidgets.QPushButton("replot")
        ##b_hide.clicked.connect(lambda: self.form_widg_left(11, model))
        #qframe_b = QtWidgets.QFrame()
        #qframe_b.setMinimumSize(40,40)
        #qframe_b.setSizePolicy(sizePolicy)
        #qframe_s = QtWidgets.QFrame()
        #qframe_s.setMinimumSize(40,2)
        #qframe_s.setSizePolicy(sizePolicy)
        ##self.lay_cpanel.addWidget(self.widg_ph)
        #self.lay_cpanel.addWidget(qframe_s)
        #self.lay_cpanel.addWidget(b_spgr)
        #self.lay_cpanel.addWidget(b_nucl)
        #self.lay_cpanel.addWidget(b_adp)
        #self.lay_cpanel.addWidget(b_ms)
        #self.lay_cpanel.addWidget(b_ms_i)
        #self.lay_cpanel.addWidget(qframe_b)
        ##self.lay_cpanel.addWidget(self.widg_exp)
        #self.lay_cpanel.addWidget(qframe_s)
        #self.lay_cpanel.addWidget(b_exp)
        #self.lay_cpanel.addWidget(b_profile)
        #self.lay_cpanel.addWidget(b_scale)
        #self.lay_cpanel.addWidget(b_extinction)
        #self.lay_cpanel.addWidget(qframe_b)
        #self.lay_cpanel.addWidget(b_output)
        #self.lay_cpanel.addWidget(b_sconstr)
        #self.lay_cpanel.addStretch(1)
        #self.lay_cpanel.addWidget(b_hide)


    def form_widg_right(self):

        width_m_3 = self.width_right
        #qframe_r = QtWidgets.QFrame()
        qframe_r = f_rhochi_viewer.interactive_graph_mod_pwd.cwidg_central()
        qframe_r.setMinimumSize(width_m_3, 1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        qframe_r.setSizePolicy(sizePolicy)
        self.lay_right.addWidget(qframe_r)
        self._p_widg_graph = qframe_r

if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)

    f_dir_prog = os.getcwd()
    mainwind1 = cbuilder(f_dir_prog)

    sys.exit(app.exec_())