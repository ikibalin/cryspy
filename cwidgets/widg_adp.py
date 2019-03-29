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
        
        q_label = QtWidgets.QLabel("""
The anisotropic displacment parameters are shown in the table:
beta11*h**2 + beta22*k**2 + beta33*l**2 + 2*beta12*h*k + 2*beta13*h*l + 2*beta23*k*l""")
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        q_label.setSizePolicy(sizePolicy)            
        lay_central.addWidget(q_label)


        llab_m = ["name", "beta11", "beta22", "beta33", "beta12", "beta13", "beta23"]
        lval_type= ["text", "val", "val", "val", "val", "val", "val"]
        lab_o = "t_adp"
        lheader = ["name", "beta11", "beta22", "beta33", "beta12", "beta13", "beta23"]
        
        self.put_tab(lay_central, lab_o, model.atom, llab_m, lval_type, lheader)
        model.add_builder(self)
        
        


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_ph()
    model.add_at()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())