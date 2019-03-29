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
import os
import sys

from PyQt5 import QtWidgets
#from PyQt5 import QtGui
#from PyQt5 import QtCore

import widg_min
import cmodel 

#it should be deleted
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

        q_label = QtWidgets.QLabel("""
Introduce the name of experiment (it should be unique) and choose its type
(powder or mono crystall).""")
        lay_central.addWidget(q_label)

        self.put_le(lay_central, "le_name", model, "name", "text")

        self.put_checkb(lay_central, "cb_powder", model, "powder", "it is marked for powder and unmarked for single crystal")
        self.put_checkb(lay_central, "cb_2dpd", model, "mode_2dpd", "2d Rietveld refinement (only for powder)")

        q_label = QtWidgets.QLabel("""
Introduce the name of input and output (experimental and model measured values) data file
(input data file should be in the same directory ).""")
        lay_central.addWidget(q_label)
        
        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_input", model, "input", "text")
        b_open = QtWidgets.QPushButton("open")
        b_open.clicked.connect(lambda: self.open_file(model))
        lay_h.addWidget(b_open)
        lay_central.addLayout(lay_h)
        
        q_label = QtWidgets.QLabel("""(to create or open input file press the button 'open')""")
        lay_central.addWidget(q_label)

        self.put_le(lay_central, "le_output", model, "output", "text")
        
        q_label = QtWidgets.QLabel("""
Polarization of the beam up and down:""")
        lay_central.addWidget(q_label)
        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_pup", model, "pup", "val")
        self.put_le(lay_h, "le_pdwon", model, "pdown", "val")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)
        

    
        lay_central.addStretch(1)
        
    def open_file(self, model):
        ffile = model.input
        fdir = model.fdirxml
        
        if ((ffile == None)|(ffile == "")|(fdir == None)):
            if fdir == None:
                fdir = "./"
            lhelp, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Input file', fdir,
                                                      "Extension file (*.*)")
            if (not ok):
                return
            fdir = os.path.dirname("{}".format(lhelp))
            fdir = fdir.replace("/","\\")
            ffile = os.path.basename("{}".format(lhelp))
            model.fdirxml = fdir
            model.input = ffile
        if model.powder:
            dtype = "powder"
        else:
            dtype = "mono"

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