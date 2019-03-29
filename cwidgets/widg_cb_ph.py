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

import widg_min



class cbuilder(widg_min.cbuilder_min):
    def __init__(self, model):
        super(cbuilder, self).__init__(model)
    def init_widget(self, model):
        """
        module is reloaded for cwidget
        """
        widg_central = cwidget(model)
        self.setCentralWidget(widg_central)

class cmodel(widg_min.cmodel_min):
    """
    minimal model for cwidget
    """
    def __init__(self):
        super(cmodel, self).__init__()
        self.ph = []
        self.exp = []
    def add_ph(self):
        ph = cmodel_ph()
        self.ph.append(ph)
        self.mod_changed()

    def add_exp(self):
        exp = cmodel_exp()
        self.exp.append(exp)
        self.mod_changed()
    def del_obs(self, obs):
        """
        to add to definition of observer
        """
        super(cmodel, self).del_obs(obs)
        for ph in self.ph:
            ph.del_obs(obs)
        for exp in self.exp:
            exp.del_obs(obs)

class cmodel_ph(widg_min.cmodel_min):
    """
    minimal model for cwidget phase
    """
    def __init__(self):
        super(cmodel_ph, self).__init__()
        self.name = "phase"

class cmodel_exp(widg_min.cmodel_min):
    """
    minimal model for cwidget experiment
    """
    def __init__(self):
        super(cmodel_exp, self).__init__()
        self.name = "exp"


class cwidget(widg_min.cwidget_min):
    def __init__(self, model):
        super(cwidget, self).__init__(model)

    def init_widget(self, model):
        """
        redefine initial 'init_widget'
        """
        super(cwidget, self).init_widget(model)
        lay_central = self.layout_central

        self.cb_ph = QtWidgets.QComboBox()
        lname = [ph.name for ph in model.ph]
        lname.append("new phase")
        self.cb_ph.addItems(lname)

        for iph, ph in enumerate(model.ph):
            ph.add_obs_labs_pars(self, "cb_ph", iph, "name", None)

        model.add_builder(self)
        
        lay_central.addWidget(self.cb_ph)


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel()
    model.add_ph()
    model.add_exp()
    mainwind1 = cbuilder(model)
    """
    model.add_builder(mainwind1)
    """
    sys.exit(app.exec_())