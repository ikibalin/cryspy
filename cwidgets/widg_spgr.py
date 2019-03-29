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
        q_label = QtWidgets.QLabel("""On the page you can introduce information about crystall strucutre 
by loading it from the '.cif' file """)
        lay_central.addWidget(q_label)

        b_lfc = QtWidgets.QPushButton("press to point the '.cif' file")
        b_lfc.clicked.connect(lambda: self.load_from_cif(model))
        lay_central.addWidget(b_lfc)

        q_label = QtWidgets.QLabel("""

or by hands. Then introduce the name of phase (it should be unique) and the space groupe 
(example: 'P2(1)/n' or 'Fd-3m 2', '2' means second choise in last example).""")
        lay_central.addWidget(q_label)

        self.put_le(lay_central, "le_name", model, "name", "text")
        
        self.put_le(lay_central, "le_spgr", model, "spgr", "text")
        
        self.put_le(lay_central, "le_spgr_code", model, "spgr_code", "text")

        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addStretch(1)
        b_les = QtWidgets.QPushButton("looks on elements of symmetry for given space groupe")
        b_les.clicked.connect(self.look_elsymm)
        lay_h.addWidget(b_les)
        lay_central.addLayout(lay_h)
        
        q_label = QtWidgets.QLabel("""
Unit cell parameters: 
       a, b, c (in angstrems)       //       alpha, beta, gamma (in degrees)""")

        lay_central.addWidget(q_label)

        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_cella", model, "cella", "val")
        self.put_le(lay_h, "le_cellb", model, "cellb", "val")
        self.put_le(lay_h, "le_cellc", model, "cellc", "val")

        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        lay_h = QtWidgets.QHBoxLayout()
        self.put_le(lay_h, "le_cellalpha", model, "cellalpha", "val")
        self.put_le(lay_h, "le_cellbeta", model, "cellbeta", "val")
        self.put_le(lay_h, "le_cellgamma", model, "cellgamma", "val")
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        q_label = QtWidgets.QLabel("""(here and further mark the corresponding check boxex to find the optimal parameters
during the refinement procedure)""")
        lay_central.addWidget(q_label)
    
        lay_central.addStretch(1)
        model.add_builder(self)
    
    def load_from_cif(self, model):
        lhelp, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', "./",
                                                      "Extension file (*.cif)")
        if (not ok):
            return
        model.load_from_cif(lhelp)
        
    def look_elsymm(self):
        pass


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_ph()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())