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
from PyQt5 import QtCore

class cmodel(object):
    """
    minimal model for cwidget
    """
    def __init__(self):
        self.name = "phase"
        self.spgr = "P1"
        self.cella = [1., False, ""]
        self.cellb = [1., False, ""]
        self.cellc = [1., False, ""]
        self.cellalpha = [90., False, ""]
        self.cellbeta = [90., False, ""]
        self.cellgamma = [90., False, ""]
        self.observer_label_ind = []
    def set_value(self, label, value, ind):
        if ind == -1:
            self.__dict__[label] = value
        else:
            self.__dict__[label][ind] = value
        self.value_changed(label, ind)
    def add_observer_label_ind(self, observer, lab_obs, lab_mod, ind):
        self.observer_label_ind.append((observer, lab_obs, lab_mod, ind))
    def del_observer_label_ind(self, observer, lab_obs, lab_mod, ind):
        self.observer_label_ind.remove((observer, lab_obs, lab_mod, ind))
    def value_changed(self, label, index):
        for observer_lab_obs_lab_mod_ind in self.observer_label_ind:
            (obs, lab_o, lab_m, ind) = observer_lab_obs_lab_mod_ind
            if (lab_m == label) & (ind == index):
                if ind == -1:
                    obs.set_string(lab_o, self.__dict__[lab_m])
                elif ((ind == 0)|(ind == 2)):
                    obs.set_string(lab_o, self.__dict__[lab_m][ind])
                elif ind == 1:
                    obs.set_logic(lab_o, self.__dict__[lab_m][ind])
    def load_from_cif(self):
        pass




class cbuilder(QtWidgets.QMainWindow):
    def __init__(self):
        super(cbuilder, self).__init__()
        model = cmodel()
        widg_central = cwidget(model)
        #self.init_layout_central()
        #widget_central = QtWidgets.QWidget()
        #widget_central.setLayout(self.layout_central)
        self.setCentralWidget(widg_central)
        self.show()


class cwidget(QtWidgets.QWidget):
    def __init__(self, model):
        super(cwidget, self).__init__()
        self.init_layout_central(model)
        self.setLayout(self.layout_central)

    def init_layout_central(self, model):
        """
        make central layout
        """
        lay_central = QtWidgets.QVBoxLayout()

        lay_h = QtWidgets.QHBoxLayout()
        b_lfc = QtWidgets.QPushButton("load from cif")
        b_lfc.clicked.connect(model.load_from_cif)
        lay_h.addWidget(b_lfc)
        self.show_text(lay_h, "name", "le_name", "name", model, -1)

        lay_central.addLayout(lay_h)

        self.show_text(lay_central, "space groupe", "le_spgr", "spgr", model, -1)

        q_label = QtWidgets.QLabel("Unit cell:")
        lay_central.addWidget(q_label)
        lay_h = QtWidgets.QHBoxLayout()
        self.show_text(lay_h, "", "le_cella", "cella", model, 0)
        self.show_text(lay_h, "", "le_cellb", "cellb", model, 0)
        self.show_text(lay_h, "", "le_cellc", "cellc", model, 0)
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        lay_h = QtWidgets.QHBoxLayout()
        self.show_text(lay_h, "", "le_cellalpha", "cellalpha", model, 0)
        self.show_text(lay_h, "", "le_cellbeta", "cellbeta", model, 0)
        self.show_text(lay_h, "", "le_cellgamma", "cellgamma", model, 0)
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)
        lay_central.addStretch(1)

        self.layout_central = lay_central
    def set_string(self, lab_o, value):
        self.__dict__[lab_o].setText("{:}".format(value))
    def set_logic(self, lab_o, value):
        self.__dict__[lab_o].setCheckState(2*value)
    def show_text(self, lay, name_for_label, name_for_observer, name_for_model, model, ind):
        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addStretch(1)
        if name_for_label != "":
            q_label = QtWidgets.QLabel(name_for_label)
            lay_h.addWidget(q_label)

        self.__dict__[name_for_observer] = QtWidgets.QLineEdit()
        if ind == -1:
            self.set_string(name_for_observer, "{:}".format(model.__dict__[name_for_model]))
        else:
            self.set_string(name_for_observer, "{:}".format(model.__dict__[name_for_model][ind]))
            if model.__dict__[name_for_model][2] != "":
                self.__dict__[name_for_observer].setReadOnly(True)
                self.__dict__[name_for_observer].setStyleSheet("color: rgb(156, 156, 156);")
            else:
                self.__dict__[name_for_observer].setReadOnly(False)
                self.__dict__[name_for_observer].setStyleSheet("color: rgb(0, 0, 0);")

            name_for_observer_cb = name_for_observer + "_cb"
            self.__dict__[name_for_observer_cb] = QtWidgets.QCheckBox()
            self.set_logic(name_for_observer_cb, model.__dict__[name_for_model][1])
        self.__dict__[name_for_observer].setFrame(False)
        self.__dict__[name_for_observer].setAlignment(QtCore.Qt.AlignRight)
        self.__dict__[name_for_observer].setTextMargins(0, 0, 0, 0)

        model.add_observer_label_ind(self, name_for_observer, name_for_model, ind)
        self.__dict__[name_for_observer].editingFinished.connect(lambda:
            model.set_value(name_for_model, self.__dict__[name_for_observer].text(), ind))
        lay_h.addWidget(self.__dict__[name_for_observer])
        if ind != -1:
            model.add_observer_label_ind(self, name_for_observer_cb, name_for_model, 1)
            self.__dict__[name_for_observer_cb].stateChanged.connect(lambda:
                model.set_value(name_for_model, self.__dict__[name_for_observer_cb].checkState() != 0, 1))
            lay_h.addWidget(self.__dict__[name_for_observer_cb])

        lay.addLayout(lay_h)
    def remove_observers(self, model):
        """
        delete observers from model
        """
        lmod_name = ["name", "spgr", "cella", "cellb", "cellc",
                     "cellalpha", "cellbeta", "cellgamma"]
        lind = [-1, -1, 0, 0, 0, 0, 0, 0]
        for mod_name, ind in zip(lmod_name, lind):
            model.del_observer_label_ind(self, "le_"+mod_name, mod_name, ind)
            if ind == 0:
                model.del_observer_label_ind(self, "le_"+mod_name+"_cb", mod_name, 1)

if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    mainwind1 = cbuilder()

    sys.exit(app.exec_())