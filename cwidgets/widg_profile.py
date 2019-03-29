# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 09:52:24 2019

Создает окно, которое используется в главной программе как виджит.

Для общности определные три класса:

    компоновщик (формирует локальное главное окно, если только этот файл запускается)
    модель (для передачи данных представителю и изменения данных при сигнале от представителя)
    представление (дает представление от модели)

@author: ikibalin
"""

import sys

from PyQt5 import QtWidgets
#from PyQt5 import QtGui
#from PyQt5 import QtCore

import widg_min
import cmodel 

import os
import wind_create_file


class cbuilder(widg_min.cbuilder_min):
    def __init__(self, model):
        super(cbuilder, self).__init__(model)
    def init_widget(self, model):
        widg_central = cwidget(model)
        self.setCentralWidget(widg_central)


class cwidget(widg_min.cwidget_min):
    def __init__(self, model):
        super(cwidget, self).__init__(model)

    def init_widget(self, model):
        """
        make central layout
        """
        super(cwidget, self).init_widget(model)
        lay_central = self.layout_central
        
        if (not model.powder):
            q_label = QtWidgets.QLabel("""
The profile can be defined only for powder type of experiment. 
For measurements with single crystal it is not neaded.""")
            lay_central.addWidget(q_label)
            return
            

        q_label = QtWidgets.QLabel("""
Introduce the name of file which contents information about background.""")
        lay_central.addWidget(q_label)


        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_bkgr", model, "bkgr", "text")
        b_open = QtWidgets.QPushButton("open")
        b_open.clicked.connect(lambda: self.open_file(model))
        lay_h.addWidget(b_open)
        lay_central.addLayout(lay_h)
        
        if model.mode_2dpd:
            q_label = QtWidgets.QLabel("""
Minimal and maximal values for the diffraction angle (ttheta, in degrees):""")
            lay_central.addWidget(q_label)
            self.put_le(lay_central, "le_tth_min", model, "tth_min", "text")
            self.put_le(lay_central, "le_tth_max", model, "tth_max", "text")
            q_label = QtWidgets.QLabel("""and for the polar angle (phi, in degrees):""")
            lay_central.addWidget(q_label)
            self.put_le(lay_central, "le_phi_min", model, "phi_min", "text")
            self.put_le(lay_central, "le_phi_max", model, "phi_max", "text")
        else:
            q_label = QtWidgets.QLabel("""
You can exclude one or several regions (example: '(0,10)(15,25)').""")
            lay_central.addWidget(q_label)
            self.put_le(lay_central, "le_excl2theta", model, "excl2theta", "text")
        
        q_label = QtWidgets.QLabel("""

Which diffraction profiles should be taken for refinement:""")
        lay_central.addWidget(q_label)
        lay_h = QtWidgets.QHBoxLayout()
        self.put_checkb(lay_h, "cb_modechi2_up", model, "modechi2_up", "up")
        self.put_checkb(lay_h, "cb_modechi2_down", model, "modechi2_down", "down")
        self.put_checkb(lay_h, "cb_modechi2_diff", model, "modechi2_diff", "diff")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)


        q_label = QtWidgets.QLabel("""
Resolution of diffractomete is defined by U, V, W, x and y.
(Pearson function 7)
HG**2 = U*tan(theta)**2 + V*tan(theta) + W + Ig/???,
HL = X*??? + Y/????.""")
        lay_central.addWidget(q_label)

        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_U", model, "U", "val")
        self.put_le(lay_h, "le_V", model, "V", "val")
        self.put_le(lay_h, "le_W", model, "W", "val")
        self.put_le(lay_h, "le_x", model, "x", "val")
        self.put_le(lay_h, "le_y", model, "y", "val")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        q_label = QtWidgets.QLabel("""
Zeroshifting is given by two parameters.
shifting = p1 + p2*???.""")
        lay_central.addWidget(q_label)

        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_zeroshift", model, "zeroshift", "val")
        self.put_le(lay_h, "le_zshift_a", model, "zshift_a", "val")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        q_label = QtWidgets.QLabel("""
Assymetry of reflection is given by four parameters.
assymetry = p1 * ??? + p2*??? + p3*??? + p4*???.""")
        lay_central.addWidget(q_label)

        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_p_asym_1", model, "p_asym_1", "val")
        self.put_le(lay_h, "le_p_asym_2", model, "p_asym_2", "val")
        self.put_le(lay_h, "le_p_asym_3", model, "p_asym_3", "val")
        self.put_le(lay_h, "le_p_asym_4", model, "p_asym_4", "val")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        lay_central.addStretch(1)
        

    def open_file(self, model):
        ffile = model.bkgr
        fdir = model.fdirxml
        
        if ((ffile == None)|(ffile == "")|(fdir == None)):
            if fdir == None:
                fdir = "./"
            lhelp, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Input file', fdir,
                                                      "Extension file (*.*)")
            if (not ok):
                return
            fdir = os.path.dirname("{}".format(lhelp))
            #fdir = fdir.replace("/","\\")
            ffile = os.path.basename("{}".format(lhelp))
            model.fdirxml = fdir
            model.input = ffile
        dtype = "bkgr"
        try:
            widg_create_file = wind_create_file.cwidg_central(fdir, ffile, dtype)
            self.redef_widget(widg_create_file, model)
        except:
            self.dialog_info("The problem with directory {:} or file {:}".format(fdir, ffile))
        
    def redef_widget(self, widget, model):
        super(cwidget, self).init_widget(model)
        lay_central = self.layout_central
        lay_central.addWidget(widget)


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_exp()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())