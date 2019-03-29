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
Information about magnetic atoms are give on the page
chi is matrix with elements 
(chi11, chi12, chi13 // chi12, chi22, chi23 // chi13, chi23, chi33 )

kappa describe the expansion (<1) or contraction (>1) of radial function R(kappa * r)

l-factor is Lande factor which defines form factor: ff = j0 + (1-2/l-factor)*j2
""")
        lay_central.addWidget(q_label)
        
        lay_h = QtWidgets.QHBoxLayout()
        b_add_at_mag = QtWidgets.QPushButton("add magnetic atom")
        b_add_at_mag.clicked.connect(lambda: self.add_at_magn(model))
        lay_h.addWidget(b_add_at_mag)
        b_del_at_mag = QtWidgets.QPushButton("delete magnetic atom")
        b_del_at_mag.clicked.connect(lambda: self.del_at_magn(model))
        lay_h.addWidget(b_del_at_mag)
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)
        q_label = QtWidgets.QLabel("""(by the following buttons the magnetic atoms can be added or exlucded from the table )""")
        lay_central.addWidget(q_label)

        llab_m = ["name", "type", "chi11", "chi22", "chi33", "chi12", "chi13", "chi23", "kappa", "lfactor"]
        lval_type= ["text", "text", "val", "val", "val", "val", "val", "val", "val", "val"]
        lab_o = "t_mag"
        lheader = ["name", "type", "chi11", "chi22", "chi33", "chi12", "chi13", "chi23", "kappa", "l-factor"]
        
        l_at_mag = [atom for atom in model.atom if atom.modemagn]
        self.put_tab(lay_central, lab_o, l_at_mag, llab_m, lval_type, lheader)
        model.add_builder(self)
        
    def add_at_magn(self, model):
        lat_name = [at.name for at in model.atom]
        text, ok = self.dialog_ask("Which atoms should be setted as magnetic: {:}".format(", ".join(lat_name)))
        if (not ok):
            return
        lname = [at_name.strip() for at_name in text.split(",")]
        for iat, at_name in enumerate(lat_name):
            if at_name in lname:
                #totaly incorrect solution but it is simple !!!!!!!!!!!!!!!!!!!
                model.atom[iat].modemagn = True
        model.mod_changed()
        
    def del_at_magn(self, model):
        lat_name = [at.name for at in model.atom]
        text, ok = self.dialog_ask("Which atoms should be setted as non-magnetic: {:}".format(", ".join(lat_name)))
        if (not ok):
            return
        lname = [at_name.strip() for at_name in text.split(",")]
        for iat, at_name in enumerate(lat_name):
            if at_name in lname:
                #totaly incorrect solution but it is simple !!!!!!!!!!!!!!!!!!!
                model.atom[iat].modemagn = False
        model.mod_changed()
        
        


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_ph()
    model.add_at()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())