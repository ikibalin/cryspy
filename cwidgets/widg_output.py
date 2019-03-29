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




class cbuilder(widg_min.cbuilder_min):
    def __init__(self, model):
        super(cbuilder, self).__init__(model)
    def init_widget(self, model):
        widg_central = cwidget(model)
        self.setCentralWidget(widg_central)


class cwidget(widg_min.cwidget_min):
    def __init__(self, model):
        super(cwidget, self).__init__(model)

    def init_widget(self, model_full):
        """
        make central layout
        """
        model = model_full.ref[0]
        super(cwidget, self).init_widget(model)
        lay_central = self.layout_central
        

        q_label = QtWidgets.QLabel("""
Introduce the output file where detailed listing will be written.""")
        lay_central.addWidget(q_label)


        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_output", model, "output", "text")
        b_open = QtWidgets.QPushButton("open")
        b_open.clicked.connect(lambda: self.open_out(model))
        lay_h.addWidget(b_open)
        lay_central.addLayout(lay_h)

        self.put_checkb(lay_central, "cb_refin", model, "refin", "make refinement")
        self.put_checkb(lay_central, "cb_sigmas", model, "sigmas", "calculate errorbars")
        

        paramref = model.paramref 
        paramrefval = model.paramrefval 
        paramreferrors = model.paramreferrors 
        if ((paramref != None)&(paramrefval != None)):
            lsout = []
            for exp in model_full.exp:
                lsout.append("'{:}' chi2/n_points is {:.3f}".format(exp.name, exp.valchi2*1.0/exp.n_total))
            qlabel = QtWidgets.QLabel("\nChi2 for experiment \n\n {:}".format("\n ".join(lsout)))
            lay_central.addWidget(qlabel)
            
            tab = QtWidgets.QTableWidget()

            tab.setRowCount(len(paramref))
            tab.setColumnCount(3)
            tab.setHorizontalHeaderLabels(["parameter", "value", "errorbar"])
            if paramreferrors == None:
                paramreferrors = len(paramref)*["--"]
            else:
                paramreferrors  = ["{:.5f}".format(hh) for hh in paramreferrors]
            paramrefval = ["{:.5f}".format(hh) for hh in paramrefval]

            for ih, ref, val, error in zip(range(len(paramref)), paramref, paramrefval, paramreferrors):
                twi = QtWidgets.QTableWidgetItem()
                twi.setText(str(ref))
                tab.setItem(ih, 0, twi)

                twi = QtWidgets.QTableWidgetItem()
                twi.setText(str(val))
                tab.setItem(ih, 1, twi)

                twi = QtWidgets.QTableWidgetItem()
                twi.setText(str(error))
                tab.setItem(ih, 2, twi)

            tab.horizontalHeader().setStretchLastSection(True);
            tab.resizeColumnsToContents()
            lay_central.addWidget(tab)
        else:
            lay_central.addStretch(1)
        
    def open_out(self, model):
        fname = model.output
        fdir = model.fdirxml
        if (fdir != None)|(fname != None)|(fname != ""):
            os.startfile("{:} ".format(os.path.join(fdir,fname)))

if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_ref()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())