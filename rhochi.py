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


from cwidgets import ccore
from cwidgets import cmodel
from cwidgets import widg_min
from cwidgets import widg_cb_ph
from cwidgets import widg_cb_exp
from cwidgets import widg_spgr
from cwidgets import widg_nucl
from cwidgets import widg_adp
from cwidgets import widg_magn
from cwidgets import widg_magn_ising
from cwidgets import widg_exptype
from cwidgets import widg_profile
from cwidgets import widg_scale
from cwidgets import widg_extinction
from cwidgets import widg_output
from cwidgets import widg_sconstr


from cwidgets import read_rcif

# =============================================================================
# #rewrite it
# =============================================================================
from cwidgets import wind_graph_inp
from cwidgets import interactive_graph_mod_pwd

from cwidgets import look_matrix
from cwidgets import interactive_graph_mod_mono



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
        #try:
        lmessage = ccore.run_refinement()
        model = cmodel.cmodel()
        model.take_from_core(ccore)
        self.signal.model = model
        #except:
        #    lmessage = ["Undefined mistake during the refinement procedure.\nSomething wrong."]
        print "Calculations are stopped"
        self.signal.message = lmessage
        self.signal.rename_signal.emit()

class cthread_signal(QtCore.QObject):
    rename_signal = QtCore.pyqtSignal()

class cbuilder(QtWidgets.QMainWindow):
    def __init__(self, model, fdirprog):
        super(cbuilder, self).__init__()
        self.thread = mythread()
        self.signal = cthread_signal()
        self.thread.signal = self.signal
        self.thread.signal.rename_signal.connect(lambda: self.reload_model(self.thread.signal))
        self.model = model
        self.fdirprog = fdirprog
        self.def_actions()
        self.init_widget()

        #self.setStyleSheet("background-color:lightgray;")
        self.show()

    def def_actions(self):
        model = self.model
        fdirprog = self.fdirprog
        fdirprog_icon = os.path.join(fdirprog,'icon_gui')

        openxmlAction = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'open.png')),'&Open XML', self)
        openxmlAction.setShortcut('Ctrl+O')
        openxmlAction.setStatusTip('Open XML file')
        openxmlAction.triggered.connect(self.open_xml)

        savexmlAction = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'save.png')),'&Save', self)
        savexmlAction.setShortcut('Ctrl+S')
        savexmlAction.setStatusTip('save information from GUI to variable and xml file')
        savexmlAction.triggered.connect(lambda: self.save_to_xml(False))

        saveasxmlAction = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'saveas.png')),'&Save as ..', self)
        saveasxmlAction.setStatusTip('save information from GUI to variable and xml file')
        saveasxmlAction.triggered.connect(lambda: self.save_to_xml(True))

        refineAction = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'calc.png')),'&Calculation', self)
        refineAction.setShortcut('Ctrl+R')
        refineAction.setStatusTip('Model optimization of marked parameters')
        refineAction.triggered.connect(self.calc_model)

        del_ph_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'del_ph.png')),'delete &phase', self)
        del_ph_action.setStatusTip('Choose phase to delete')
        del_ph_action.triggered.connect(self.del_ph)

        del_exp_action = QtWidgets.QAction(QtGui.QIcon(os.path.join(fdirprog_icon,'del_exp.png')),'delete &experiment', self)
        del_exp_action.setStatusTip('Choose experiment to delete')
        del_exp_action.triggered.connect(self.del_exp)

        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openxmlAction)
        fileMenu.addAction(savexmlAction)
        fileMenu.addAction(saveasxmlAction)

        calcMenu = menubar.addMenu('&Calculations')
        calcMenu.addAction(refineAction)
        calcMenu.addAction(del_ph_action)
        calcMenu.addAction(del_exp_action)


        self.toolbar = self.addToolBar("Open")
        self.toolbar.addAction(openxmlAction)
        self.toolbar.addAction(savexmlAction)
        self.toolbar.addAction(saveasxmlAction)
        self.toolbar.addAction(refineAction)
        self.toolbar.addAction(del_ph_action)
        self.toolbar.addAction(del_exp_action)

    def del_ph(self):
        model = self.model
        lph_name = [ph.name for ph in model.ph]
        question = "Which phase should be deleted: {:}".format(", ".join(lph_name))
        text, ok = QtWidgets.QInputDialog.getText(self, 'Question',question)
        if (not ok):
            return
        lname = [ph_name.strip() for ph_name in text.split(",")]
        for iph, ph_name in enumerate(lph_name):
            if ph_name in lname:
                #totaly incorrect solution but it is simple !!!!!!!!!!!!!!!!!!!
                ph_del = model.ph[iph]
                model.ph.remove(ph_del)
        model.mod_changed()
        
    def del_exp(self):
        model = self.model
        lexp_name = [exp.name for exp in model.exp]
        question = "Which experiment should be deleted: {:}".format(", ".join(lexp_name))
        text, ok = QtWidgets.QInputDialog.getText(self, 'Question',question)
        if (not ok):
            return
        lname = [exp_name.strip() for exp_name in text.split(",")]
        for iexp, exp_name in enumerate(lexp_name):
            if exp_name in lname:
                #totaly incorrect solution but it is simple !!!!!!!!!!!!!!!!!!!
                exp_del = model.exp[iexp]
                model.exp.remove(exp_del)
        model.mod_changed()
    def init_widget(self): 
        """
        module is reloaded for cwidget
        """
        model = self.model
        self.location_on_the_screen()
        widg_central = cwidget(model, self.info_width)
        self.setCentralWidget(widg_central)
        self.setWindowTitle('RhoChi')

    def open_xml(self):
        model = self.model
        fdirprog = self.fdirprog
        fdirxml = model.fdirxml
        if fdirxml == None:
            fdirxml = os.getcwd()
        lhelp, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdirxml,
                                                      "Extension file (*.rcif)")
        if (not ok):
            return
        fdirxml = os.path.dirname("{}".format(lhelp))
        
        fxml = os.path.basename("{}".format(lhelp))

        self.setWindowTitle('RhoChi - {:}'.format(lhelp))

        core = read_rcif.get_core_from_file(os.path.join(fdirxml, fxml))
        core.set_fdirxml(fdirxml)
        core.set_fdirprog(fdirprog) 
        core.set_fxml(fxml)        
        
        model = cmodel.cmodel()
        model.take_from_core(core)

        self.model = model
        self.init_widget()

    def save_to_xml(self, saveas = False):
        """
        save information to xml file
        """
        model = self.model
        fxml = model.fxml
        fdirxml = model.fdirxml
        fdirprog = model.fdirprog
        if ((fxml == None) | (saveas)):
            lhelp, ok =QtWidgets.QFileDialog.getSaveFileName()
            print ok
            """
            if (not (ok)):
                return
            """
            fxml = os.path.basename("{}".format(lhelp))
            fdirxml = os.path.dirname("{}".format(lhelp))
            self.setWindowTitle('RhoChi - {:}'.format(lhelp))
            
        
        model.set_fdirprog(fdirprog)
        model.set_fdirxml(fdirxml)
        model.set_fxml(fxml)
        #redo it: core should be independent from model
        core = ccore.ccore(model)
        read_rcif.save_core_to_file(core, os.path.join(fdirxml, fxml))

    def calc_model(self):
        #self.setStyleSheet("background-color:lightyellow;")
        model = self.model
        model.check_phaseexp()
        core = ccore.ccore(model)
        self.thread.core = core
        self.setWindowTitle('CALCULATIONS ARE RUNNING............')
        self.thread.start()

    def reload_model(self, obj):
        #self.setStyleSheet("background-color:lightgray;")
        try:
            lmessage = obj.message
            if lmessage!= []:
                print lmessage, obj.message
                QtWidgets.QMessageBox.question(self, 'MISATAKE', "The cacluclations are stopped because of:\n\n"+"\n".join(lmessage))
                obj.message = []
                return
        except:
            return
        self.setWindowTitle('RhoChi')
        try:
            model = obj.model
            self.model = model
            self.init_widget()
        except:
            return

    def location_on_the_screen(self):
        """
        position and size of main window
        """
        screen = QtWidgets.QDesktopWidget().screenGeometry()
        self.setMinimumSize(screen.width()*3/4, screen.height()*3/4)
        #self.setMaximumSize(screen.width()*3/4, screen.height()*3/4)
        self.info_width = screen.width()*3/4
        self.move(screen.width()/8 , screen.height()/16)






class cwidget(widg_min.cwidget_min):
    def __init__(self, model, width):
        self.info_width = width
        super(cwidget, self).__init__(model)
    def init_widget(self, model):
        """
        make central layout
        """
        lay_central = self.layout()
        if lay_central == None:
            lay_central = QtWidgets.QHBoxLayout()
        else:
            del_layout(lay_central)
        self.layout_central = lay_central


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


        self.form_widg_cpanel(model)
        self.lay_left.addWidget(QtWidgets.QFrame())
        self.form_widg_right(model)
        


    def form_widg_cpanel(self, model):

        for i in range(self.lay_cpanel.count()):
            widg = self.lay_cpanel.itemAt(i).widget()
            model.del_obs(widg)
            model.del_builder(widg)

        del_layout(self.lay_cpanel)

        width_m_1 = self.width_cpanel
        """
        self.widg_ph_exp = widg_cb_ph_exp.cwidget(model)


        self.cb_ph = QtWidgets.QComboBox()
        lname = [ph.name for ph in model.ph]
        lname.append("new phase")
        self.cb_ph.addItems(lname)



        for iph, ph in enumerate(model.ph):
            ph.add_obs_labs_pars(self, "cb_ph", iph, "name", None)
        """

        self.widg_ph = widg_cb_ph.cwidget(model)
        self.widg_exp = widg_cb_exp.cwidget(model)
        

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)

        b_spgr = QtWidgets.QPushButton("unit cell")
        b_spgr.setMaximumSize(width_m_1, 50)
        b_spgr.setSizePolicy(sizePolicy)
        b_spgr.clicked.connect(lambda: self.form_widg_left(0, model))
        b_nucl = QtWidgets.QPushButton("nucl. structure")
        b_nucl.clicked.connect(lambda: self.form_widg_left(1, model))
        b_adp = QtWidgets.QPushButton("anis.displacement")
        b_adp.clicked.connect(lambda: self.form_widg_left(2, model))
        b_ms = QtWidgets.QPushButton("magnetic structure")
        b_ms.clicked.connect(lambda: self.form_widg_left(3, model))
        
        b_ms_i = QtWidgets.QPushButton("Ising model")
        b_ms_i.clicked.connect(lambda: self.form_widg_left(4, model))


        b_exp = QtWidgets.QPushButton("experiment")
        b_exp.clicked.connect(lambda: self.form_widg_left(5, model))
        b_profile = QtWidgets.QPushButton("profile")
        b_profile.clicked.connect(lambda: self.form_widg_left(6, model))
        b_scale = QtWidgets.QPushButton("scale")
        b_scale.clicked.connect(lambda: self.form_widg_left(7, model))
        b_extinction = QtWidgets.QPushButton("extinction")
        b_extinction.clicked.connect(lambda: self.form_widg_left(8, model))
        b_output = QtWidgets.QPushButton("output")
        b_output.clicked.connect(lambda: self.form_widg_left(9, model))

        b_sconstr = QtWidgets.QPushButton("set constraints")
        b_sconstr.clicked.connect(lambda: self.form_widg_left(10, model))
        b_hide = QtWidgets.QPushButton("replot")
        b_hide.clicked.connect(lambda: self.form_widg_left(11, model))

        qframe_b = QtWidgets.QFrame()
        qframe_b.setMinimumSize(40,40)
        qframe_b.setSizePolicy(sizePolicy)

        qframe_s = QtWidgets.QFrame()
        qframe_s.setMinimumSize(40,2)
        qframe_s.setSizePolicy(sizePolicy)

        self.lay_cpanel.addWidget(self.widg_ph)
        self.lay_cpanel.addWidget(qframe_s)
        self.lay_cpanel.addWidget(b_spgr)
        self.lay_cpanel.addWidget(b_nucl)
        self.lay_cpanel.addWidget(b_adp)
        self.lay_cpanel.addWidget(b_ms)
        self.lay_cpanel.addWidget(b_ms_i)

        self.lay_cpanel.addWidget(qframe_b)

        self.lay_cpanel.addWidget(self.widg_exp)
        self.lay_cpanel.addWidget(qframe_s)
        self.lay_cpanel.addWidget(b_exp)
        self.lay_cpanel.addWidget(b_profile)
        self.lay_cpanel.addWidget(b_scale)
        self.lay_cpanel.addWidget(b_extinction)

        self.lay_cpanel.addWidget(qframe_b)


        self.lay_cpanel.addWidget(b_output)
        self.lay_cpanel.addWidget(b_sconstr)

        self.lay_cpanel.addStretch(1)
        self.lay_cpanel.addWidget(b_hide)


    def form_widg_left(self, row, model):
        #form the left widget
        ind_ph = self.widg_ph.cb_ph.currentIndex()
        ind_exp = self.widg_exp.cb_exp.currentIndex()
        ind_ref = 0


        lclass = [widg_spgr.cwidget, widg_nucl.cwidget, widg_adp.cwidget, widg_magn.cwidget,
                  widg_magn_ising.cwidget,
                  widg_exptype.cwidget, widg_profile.cwidget, widg_scale.cwidget,
                  widg_extinction.cwidget, widg_output.cwidget, widg_sconstr.cwidget]
        lparam = [("ph", ind_ph), ("ph", ind_ph), ("ph", ind_ph), ("ph", ind_ph),
                  ("ph", ind_ph),
                  ("exp", ind_exp), ("exp", ind_exp), ("exp", ind_exp),
                  ("exp", ind_exp), ("", None), ("", None)]

        for i in range(self.lay_left.count()):
            widg = self.lay_left.itemAt(i).widget()
            model.del_obs(widg)
            model.del_builder(widg)

        del_layout(self.lay_left)

        if row < len(lclass):
            (name, ind) = lparam[row]
            if name == "":
                widg = lclass[row](model)
            else:
                if ind >= len(model.__dict__[name]):
                    ind = len(model.__dict__[name])
                    if name == "ph":
                        model.add_ph()
                        self.widg_ph.cb_ph.setCurrentIndex(ind_ph)
                    elif name == "exp":
                        model.add_exp()
                        self.widg_exp.cb_exp.setCurrentIndex(ind_exp)
                widg = lclass[row](model.__dict__[name][ind])

        else:
            widg = QtWidgets.QFrame()
            if row == 11:
                self.form_widg_right(model)
            

        self.lay_left.addWidget(widg)
        

    def form_widg_right(self, model):
        width_m_3 = self.width_right
        ind_exp = self.widg_exp.cb_exp.currentIndex()

        
        for i in range(self.lay_right.count()):
            if self.lay_left.itemAt(i) != None:
                widg = self.lay_left.itemAt(i).widget()
                model.del_obs(widg)
                model.del_builder(widg)
        del_layout(self.lay_right)
        qframe_r = QtWidgets.QFrame()
        qframe_r.setMinimumSize(width_m_3, 1)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        qframe_r.setSizePolicy(sizePolicy)
        self.lay_right.addWidget(qframe_r)
        try:
            if (model.fdirxml != None)&(len(model.exp) > ind_exp):
                fdir = model.fdirxml
                finp = model.exp[ind_exp].input
                fout = model.exp[ind_exp].output
                ffull_inp = os.path.join(fdir, finp)
                ffull_out = os.path.join(fdir, fout)
                if os.path.exists(ffull_out):
                    if (model.exp[ind_exp].powder):
                        if model.exp[ind_exp].mode_2dpd:
                            fname_e = os.path.join(fdir, "e_"+fout)
                            fname_m = os.path.join(fdir, "m_"+fout)
                            widg = look_matrix.cwidg_central(fname_e, fname_m)
                        else:
                            widg = interactive_graph_mod_pwd.cwidg_central(ffull_out)
                    else:
                        widg = interactive_graph_mod_mono.cwidg_central(ffull_out)
                    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
                    widg.setSizePolicy(sizePolicy)
                    self.lay_right.addWidget(widg)
                elif os.path.exists(ffull_inp):
                    if (model.exp[ind_exp].powder):
                        widg = wind_graph_inp.cwidg_central(ffull_inp)
                    sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
                    widg.setSizePolicy(sizePolicy)
                    self.lay_right.addWidget(widg)
        except:
            pass


if __name__ == '__main__':
    larg = sys.argv

    fdirprog = os.path.dirname(os.path.abspath(larg[0]))
    print fdirprog 
    if len(larg) >= 2:
        fdirxml = os.path.abspath(os.path.dirname(larg[1]))
        fxml = os.path.basename(larg[1])
    else:
        fdirxml = fdirprog 
        fxml = None
        
    flag_gui = True
    if len(larg) >= 3:
        #consol version
        flag_gui = (larg[2] != "consol")
        print "Program is runned in the console mode."

    if fxml != None:
        core = read_rcif.get_core_from_file(os.path.join(fdirxml, fxml))
    else:
        core = ccore.ccore()
        
    core.set_fdirxml(fdirxml)
    core.set_fdirprog(fdirprog) 
    core.set_fxml(fxml)
    
    if flag_gui:
        print 5*"\n"
        model = cmodel.cmodel()
        model.take_from_core(core)
        
        app = QtWidgets.QApplication(larg)
        mainwind1 = cbuilder(model, fdirprog)
        #model.add_builder(mainwind1)
        sys.exit(app.exec_())
    else:
        print "Calculations are started\n"
        lmessage = core.run_refinement()
        print 80*"*"
        for hh in lmessage:
            print hh
        if lmessage == []:
            fxml_out = "o_{:}".format(fxml)
            print "the refined parameters are saved in the file: '{:}'".format(fxml_out)
            read_rcif.save_core_to_file(core, os.path.join(fdirxml, fxml_out))
        print 80*"*"
        print "That is all."
        print 80*"*"
    