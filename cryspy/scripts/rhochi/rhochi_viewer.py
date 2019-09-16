# -*- coding: utf-8 -*-
"""
Simple RhoChi viewer
"""
__author__ = 'ikibalin'
__version__ = "2019_09_10"

import os
import os.path
import sys
import numpy

import pycifstar

f_dir = os.path.dirname(__file__)
os.chdir(f_dir)
sys.path.insert(0, f_dir)

from cryspy import RhoChi, rhochi_refinement
from cryspy.scripts.rhochi.cl_rhochi import create_temporary

from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore


import interactive_graph_mod_mono
import interactive_graph_mod_pwd
import interactive_graph_mod_pwd_2d_as_1d



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

def multistring_to_array(string):
    l_1 = string.strip().split("\n")
    l_ttheta = [float(_) for _ in l_1[0].strip().split()[1:]]
    l_phi, ll_intensity = [], []
    for line in l_1[1:]:
        l_1 = line.strip().split()
        l_phi.append(float(l_1[0]))
        ll_intensity.append([float(_) if _ != "None" else None for _ in l_1[1:]])
    x = numpy.array(l_ttheta, dtype=float)
    y = numpy.array(l_phi, dtype=float)
    zz = (numpy.array(ll_intensity, dtype=float)).transpose()
    return x, y, zz


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
        f_file_prog_setup = os.path.join(f_dir_prog, 'setup.dat')

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
        f_file_prog_setup = os.path.join(f_dir_prog, 'setup.dat')
        fid = open(f_file_prog_setup, "w")
        fid.write("\n".join(ls_out))
        fid.close()


    def def_actions(self):

        f_dir_prog = self._p_f_dir_prog
        f_dir_prog_icon = os.path.join(f_dir_prog, 'f_icon')

        create_rcif_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'create_new.png')),'Create &new', self)
        create_rcif_action.setShortcut('Ctrl+N')
        create_rcif_action.setStatusTip('Create new rcif')
        create_rcif_action.triggered.connect(self.create_rcif)

        open_rcif_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'open.png')),'&Open .rcif', self)
        open_rcif_action.setShortcut('Ctrl+O')
        open_rcif_action.setStatusTip('Open data (1D)')
        open_rcif_action.triggered.connect(self.open_rcif)


        run_rhochi = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'calc_rcif.png')),'&Run RhoChi', self)
        run_rhochi.setShortcut('Ctrl+R')
        run_rhochi.setStatusTip('Model optimization of marked parameters')
        run_rhochi.triggered.connect(self.run_rhochi)

        read_rcif_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(f_dir_prog_icon,'read_rcif.png')),'Read .rcif', self)
        read_rcif_action.setStatusTip('Read .rcif')
        read_rcif_action.triggered.connect(self.read_rcif)



        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(create_rcif_action)
        fileMenu.addAction(open_rcif_action)

        calcMenu = menubar.addMenu('&Calculations')
        calcMenu.addAction(run_rhochi)
        calcMenu.addAction(read_rcif_action)


        self.toolbar = self.addToolBar("Open")
        self.toolbar.addAction(create_rcif_action)
        self.toolbar.addAction(open_rcif_action)
        self.toolbar.addAction(run_rhochi)
        self.toolbar.addAction(read_rcif_action)

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
            if f_name_data == "main.rcif":
                l_file = [_ for _ in os.listdir(f_dir_data) if os.path.isfile(os.path.join(f_dir_data, _))]
                l_file_calc = [_ for _ in l_file if _.endswith("_calc.rcif")]
                if len(l_file_calc) != 0:
                    f_name_data = l_file_calc[0]
            f_name_full = os.path.join(f_dir_data, f_name_data)
            cif_global = pycifstar.to_global(f_name_full)

            flag_single = cif_global.is_prefix("_diffrn_refln")
            flag_powder_1d = cif_global.is_prefix("_pd_proc")
            flag_powder_2d = cif_global.is_prefix("_pd2d_proc")
            if flag_single:
                x = numpy.array(cif_global["_diffrn_refln_fr_calc"], dtype=float)
                y_exp = numpy.array(cif_global["_diffrn_refln_fr"], dtype=float)
                y_sig = numpy.array(cif_global["_diffrn_refln_fr_sigma"], dtype=float)
                self.plot_data_mono(x, y_exp, y_sig)
            if flag_powder_1d:
                x = numpy.array(cif_global["_pd_proc_2theta"], dtype=float)
                y_1_mod = numpy.array(cif_global["_pd_proc_intensity_up_total"], dtype=float)
                y_2_mod = numpy.array(cif_global["_pd_proc_intensity_down_total"], dtype=float)
                y_1_exp = numpy.array(cif_global["_pd_proc_intensity_up"], dtype=float)
                y_1_sig = numpy.array(cif_global["_pd_proc_intensity_up_sigma"], dtype=float)
                y_2_exp = numpy.array(cif_global["_pd_proc_intensity_down"], dtype=float)
                y_2_sig = numpy.array(cif_global["_pd_proc_intensity_down_sigma"], dtype=float)

                y_3_mod = y_1_mod + y_2_mod
                y_4_mod = y_1_mod - y_2_mod
                y_3_exp = y_1_exp + y_2_exp
                y_4_exp = y_1_exp - y_2_exp
                y_3_sig = (y_1_sig**2 + y_2_sig**2)**0.5
                y_4_sig = (y_1_sig**2 + y_2_sig**2)**0.5

                l_y_mod = [y_3_mod, y_4_mod]
                l_y_exp = [y_3_exp, y_4_exp]
                l_y_sig = [y_3_sig, y_4_sig]
                self.plot_data_powder_1d(x, l_y_mod, l_y_exp, l_y_sig)
            if flag_powder_2d:
                s_1_mod = cif_global["_pd2d_proc_2theta_phi_intensity_up_total"].value
                s_2_mod = cif_global["_pd2d_proc_2theta_phi_intensity_down_total"].value
                s_1_exp = cif_global["_pd2d_proc_2theta_phi_intensity_up"].value
                s_1_sig = cif_global["_pd2d_proc_2theta_phi_intensity_up_sigma"].value
                s_2_exp = cif_global["_pd2d_proc_2theta_phi_intensity_down"].value
                s_2_sig = cif_global["_pd2d_proc_2theta_phi_intensity_down_sigma"].value
                x, _2, yy_1_mod = multistring_to_array(s_1_mod)
                x, _2, yy_2_mod = multistring_to_array(s_2_mod)
                x, _2, yy_1_exp = multistring_to_array(s_1_exp)
                x, _2, yy_2_exp = multistring_to_array(s_2_exp)
                x, _2, yy_1_sig = multistring_to_array(s_1_sig)
                x, _2, yy_2_sig = multistring_to_array(s_2_sig)
                y_1_exp = numpy.where(numpy.isnan(yy_1_exp), 0., yy_1_exp).sum(axis=1)
                y_1_mod = numpy.where(numpy.isnan(yy_1_exp), 0., yy_1_mod).sum(axis=1)
                y_1_sig = ((numpy.where(numpy.isnan(yy_1_exp), 0., yy_1_sig)**2).sum(axis=1))**0.5

                y_2_exp = numpy.where(numpy.isnan(yy_2_exp), 0., yy_2_exp).sum(axis=1)
                y_2_mod = numpy.where(numpy.isnan(yy_2_exp), 0., yy_2_mod).sum(axis=1)
                y_2_sig = ((numpy.where(numpy.isnan(yy_2_exp), 0., yy_2_sig)**2).sum(axis=1))**0.5

                y_3_mod = y_1_mod + y_2_mod
                y_4_mod = y_1_mod - y_2_mod
                y_3_exp = y_1_exp + y_2_exp
                y_4_exp = y_1_exp - y_2_exp
                y_3_sig = (y_1_sig**2 + y_2_sig**2)**0.5
                y_4_sig = (y_1_sig**2 + y_2_sig**2)**0.5

                l_y_mod = [y_3_mod, y_4_mod]
                l_y_exp = [y_3_exp, y_4_exp]
                l_y_sig = [y_3_sig, y_4_sig]
                self.plot_data_powder_1d(x, l_y_mod, l_y_exp, l_y_sig)
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

    def create_rcif(self):
        if "f_dir_data" in self._p_d_setup.keys():
            f_dir_data = self._p_d_setup["f_dir_data"]
            if not(os.path.isdir(f_dir_data)):
                f_dir_data = self._p_f_dir_prog
        else:
            f_dir_data = self._p_f_dir_prog
        #f_file_data_new = QtWidgets.QFileDialog.getExistingDirectory(self, "Select a folder:", f_dir_data)
        f_file_data_new, okPressed = QtWidgets.QFileDialog.getSaveFileName(self, "Select a folder:", "full.rcif")
        if not okPressed:
            return
        if f_file_data_new != "":
            f_dir_data = os.path.dirname(f_file_data_new)
            f_name_data = os.path.basename(f_file_data_new)
            self._p_d_setup["f_dir_data"] = f_dir_data
            self._p_d_setup["f_name_data"] = f_name_data
            self.write_setup()
            f_name_full = os.path.join(f_dir_data, f_name_data)

            i, okPressed = QtWidgets.QInputDialog.getInt(self, "Type of data","Type of experiment:\n\n 1. Single crystal\n 2. Powder 1D\n 3 .Powder 2D\n 4 .Powder 2D (one-axial texture)", 1, 1, 4, 1)
            if okPressed:
                exp_type = str(i)
                create_temporary(f_name_full, exp_type)

    def open_rcif(self):
        self.open_file(s_constr="Rcif files (*.rcif)")
        flag_out = self.ask_file(ext=".rcif")
        if not flag_out:
            return None
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
        f_file_data_new, okPressed = QtWidgets.QFileDialog.getOpenFileName(self, 'Select a file:', f_dir_data, s_constr)
        if not okPressed:
            return False
        if f_file_data_new != "":
            f_dir_data = os.path.dirname(f_file_data_new)
            f_name_data = os.path.basename(f_file_data_new)
            self._p_d_setup["f_dir_data"] = f_dir_data
            self._p_d_setup["f_name_data"] = f_name_data
            self.write_setup()
        return True
    def ask_file(self, ext=".rcif"):
        cond_ask = False
        cond_cycle = True

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
            cond_cycle = self.open_file()
        return cond_cycle

    def run_rhochi(self):
        flag_out = self.ask_file(".rcif")
        if not flag_out:
            return None
        f_name_data = self._p_d_setup["f_name_data"]
        f_dir_data = self._p_d_setup["f_dir_data"]
        f_name_full = os.path.join(f_dir_data, f_name_data)
        f_name_out = os.path.join(f_dir_data, "out.rcif")
        
        #f_dir_prog = self._p_f_dir_prog
        #f_prog_full = os.path.join(f_dir_prog, "rhochi.py")
        #line = "python {:} {:} {:}".format(f_prog_full, f_name_full, f_name_full)
        f_dir = os.path.dirname(f_name_full)
        print("Step 1: Input files in directory: ", f_dir)
        rhochi = RhoChi(file_input=f_name_full)
        print("Step 2: File reading")
        rhochi.read_files(f_dir=f_dir)
        print("Step 3: Refienment")
        rhochi.refine()
        print("Step 4: Saving to files")
        rhochi.save_to_files()
        #os.system(line)
        self.try_to_plot_from_rcif()

    def read_rcif(self):
        flag_out = self.ask_file(".rcif")
        if not flag_out:
            return None
        f_name_data = self._p_d_setup["f_name_data"]
        f_dir_data = self._p_d_setup["f_dir_data"]
        f_name_full = os.path.join(f_dir_data, f_name_data)
        os.startfile(f_name_full)

    def plot_data_mono(self, x, y_exp, y_sig):
        widg_central = self._p_widg_central
        stack_widg = widg_central._p_widg_graph
        index = 0
        stack_widg.setCurrentIndex(index)
        widg_graph = stack_widg.widget(index)
        widg_graph.plot_file(x, y_exp, y_sig)

    def plot_data_powder_1d(self, x, l_y_mod, l_y_exp, l_y_sig):
        widg_central = self._p_widg_central
        stack_widg = widg_central._p_widg_graph
        index = 1
        stack_widg.setCurrentIndex(index)
        widg_graph = stack_widg.widget(index)
        widg_graph.plot_file(x, l_y_mod, l_y_exp, l_y_sig)

    def plot_data_powder_2d(self, f_name_full):
        widg_central = self._p_widg_central
        stack_widg = widg_central._p_widg_graph
        index = 2
        stack_widg.setCurrentIndex(index)
        widg_graph = stack_widg.widget(index)
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


        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)



    def form_widg_right(self):

        width_m_3 = self.width_right
        #widg_powder_1d = QtWidgets.QFrame()
        stack_widg = QtWidgets.QStackedWidget(self)
        stack_widg.setMinimumSize(width_m_3, 1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        stack_widg.setSizePolicy(sizePolicy)


        widg_mono = interactive_graph_mod_mono.cwidg_central()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        widg_mono.setSizePolicy(sizePolicy)
        stack_widg.addWidget(widg_mono)
        
        widg_powder_1d = interactive_graph_mod_pwd.cwidg_central()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        widg_powder_1d.setSizePolicy(sizePolicy)

        stack_widg.addWidget(widg_powder_1d)

        widg_powder_2d = interactive_graph_mod_pwd_2d_as_1d.cwidg_central()
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        widg_powder_2d.setSizePolicy(sizePolicy)

        stack_widg.addWidget(widg_powder_2d)
        
        self.lay_right.addWidget(stack_widg)
        self._p_widg_graph = stack_widg
        stack_widg.setCurrentIndex(0)

def main(l_arg= []):
    app = QtWidgets.QApplication(l_arg)
    f_dir_prog = os.path.dirname(__file__)
    mainwind1 = cbuilder(f_dir_prog)
    sys.exit(app.exec_())


if __name__ == '__main__':
    l_arg = sys.argv
    main(l_arg)


