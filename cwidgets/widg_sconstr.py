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

        q_label = QtWidgets.QLabel("""Set constraints on ...""")
        lay_central.addWidget(q_label)

        """
        self.put_le(lay_central, "le_name", model, "name", "text")
        self.put_checkb(lay_central, "cb_powder", model, "powder", "it is marked for powder and unmarked for single crystal")
        """

        self.cb_ucp = QtWidgets.QCheckBox("unit cell parameters (singony)")
        self.cb_adp = QtWidgets.QCheckBox("anisotropic beta and chi (symmetry)")
        self.cb_iso = QtWidgets.QCheckBox("spherical beta and chi for atoms: ")
        self.cb_chi_equal = QtWidgets.QCheckBox("equality of chi for atoms: ")
        self.cb_flip = QtWidgets.QCheckBox("efficiency of flipper for experiments: ")


        self.le_iso_at = QtWidgets.QLineEdit()
        self.le_iso_ph = QtWidgets.QLineEdit()

        self.le_chi_at = QtWidgets.QLineEdit()
        self.le_chi_ph = QtWidgets.QLineEdit()
        self.le_flip = QtWidgets.QLineEdit()

        q_label_1 = QtWidgets.QLabel(""" in phase: """)
        q_label_2 = QtWidgets.QLabel(""" in phase: """)

        lay_central.addWidget(self.cb_ucp)
        lay_central.addWidget(self.cb_adp)


        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addWidget(self.cb_iso)
        lay_h.addWidget(self.le_iso_at)
        lay_h.addWidget(q_label_1)
        lay_h.addWidget(self.le_iso_ph)
        lay_central.addLayout(lay_h)

        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addWidget(self.cb_chi_equal)
        lay_h.addWidget(self.le_chi_at)
        lay_h.addWidget(q_label_2)
        lay_h.addWidget(self.le_chi_ph)
        lay_central.addLayout(lay_h)

        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addWidget(self.cb_flip)
        lay_h.addWidget(self.le_flip)
        lay_central.addLayout(lay_h)

        b_gen = QtWidgets.QPushButton("generate")
        b_gen.clicked.connect(lambda: self.gen_constr(model))
        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addWidget(b_gen)
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        q_label = QtWidgets.QLabel("""Full list of constraints (apllied and generated):""")
        lay_central.addWidget(q_label)


        self.t_constr = QtWidgets.QTableWidget()
        lheader = ["p1", "coeff", "p2"]
        self.t_constr.setHorizontalHeaderLabels(lheader)
        self.t_constr.horizontalHeader().setStretchLastSection(True)
        self.t_constr.setColumnCount(3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.t_constr.setSizePolicy(sizePolicy)
        lay_central.addWidget(self.t_constr)


        b_del = QtWidgets.QPushButton("delete")
        b_iconstr = QtWidgets.QPushButton("apply")
        b_iconstr.clicked.connect(lambda: self.apply_constr(model))
        b_del.clicked.connect(lambda: self.del_constr(model))
        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addWidget(b_iconstr)
        lay_h.addStretch(1)
        lay_h.addWidget(b_del)
        lay_central.addLayout(lay_h)
        
        
        self.get_link_on_constr(model)
        
    def get_link_on_constr(self, model):
        l_link = model.get_link_on_constr()
        self.t_constr.setRowCount(len(l_link))
        for iconstr, link in enumerate(l_link):
            val, lmessage = model.set_val_by_link(link, None)
            if lmessage != []:
                self.dialog_info("ERROR\n"+"\n".join(lmessage))
                return
            lhelp = val[2].split("[")
            scoeff = lhelp[0].replace("x1","").replace("*","").strip()
            if scoeff == "":
                scoeff = "1.0"
            subj = lhelp[1].replace("]","").strip()
            constr = [link, scoeff, subj]
            for hh in range(3):
                flag = False
                twi = self.t_constr.item(iconstr, hh)
                if twi == None:
                    twi = QtWidgets.QTableWidgetItem()
                    flag = True
                twi.setText(constr[hh])
                if flag:
                    self.t_constr.setItem(iconstr, hh, twi)

    def gen_constr(self, model):
        lconstr = []
        if self.cb_ucp.checkState() != 0:
            try:
                lconstr_1 = model.get_constr_ucp()
            except:
                self.dialog_info("Problem with restriction of unit cell parameters. Try to save file first and then load it.")
                return
            lconstr.extend(lconstr_1)
        if  self.cb_adp.checkState() != 0:
            try:
                lconstr_1 = model.get_constr_adp()
            except:
                self.dialog_info("Problem with restriction of anisotropic displacement parameters. Try to save file first and then load it.")
                return
            lconstr.extend(lconstr_1)
        lconstr_small = [hh for hh in lconstr if len(hh)==3]
        self.t_constr.setRowCount(len(lconstr_small))
        for iconstr, constr in enumerate(lconstr_small):
            twi = QtWidgets.QTableWidgetItem()
            twi.setText(constr[0])
            self.t_constr.setItem(iconstr, 0, twi)
            twi = QtWidgets.QTableWidgetItem()
            twi.setText(str(constr[1]))
            self.t_constr.setItem(iconstr, 1, twi)
            twi = QtWidgets.QTableWidgetItem()
            twi.setText(constr[2])
            self.t_constr.setItem(iconstr, 2, twi)
    """
    def get_constr_ucp(self, model):
        lconstr = model.get_constr_ucp()
        return lconstr
    """
    def apply_constr(self, model):
        lmessage = []
        n_r = self.t_constr.rowCount()
        if n_r < 0:
            return
        for i_r in range(n_r):
            it_obj = self.t_constr.item(i_r, 0)
            it_coeff = self.t_constr.item(i_r, 1)
            it_sub = self.t_constr.item(i_r, 2)
            if all([it_obj != None, it_coeff != None, it_sub != None]):
                link_obj = "{:}".format(it_obj.text()).strip()
                coeff = "{:}".format(it_coeff.text()).strip()
                link_sub = "{:}".format(it_sub.text()).strip()
                link = "{:}*x1 [{}]".format(coeff,link_sub)
                val_old, lmessage = model.set_val_by_link(link_sub, None)
                if lmessage != []:
                    self.dialog_info("ERROR\n"+"\n".join(lmessage))
                    return
                val = [float(coeff)*val_old[0], False, link]
                val_new, lmessage= model.set_val_by_link(link_obj, val)
                if lmessage != []:
                    self.dialog_info("ERROR\n"+"\n".join(lmessage))
                    return
            lmessage.append("Constraint {:} is applied: {:}".format(i_r, val))
        self.dialog_info("\n".join(lmessage))
            
            
    def del_constr(self, model):
        text, ok = self.dialog_ask("Which constraints should be deleted \n(example: 1-3, 5, 7 or all)")
        n_r = self.t_constr.rowCount()
        lnumb = []
        if text.strip() == "all":
            lnumb = range(n_r)
        else:
            lhelp = text.split(",")
            for srange in lhelp:
                if srange.strip().isdigit():
                    lnumb.append(int(srange)-1)
                else:
                    lhelp2 = srange.split("-")
                    if len(lhelp2) !=2:
                        pass
                    elif all([lhelp2[0].strip().isdigit(), lhelp2[1].strip().isdigit()]):
                        n1, n2 = int(lhelp2[0]), int(lhelp2[1])
                        lnumb.extend(range(n1-1,n2))
        lnumb = [numb for numb in lnumb if numb < n_r]
        for numb in lnumb:
            twi = self.t_constr.item(numb, 0)
            slink = twi.text()
            val_old, lmessage = model.set_val_by_link(slink, None)
            if lmessage != []:
                self.dialog_info("\n".join(lmessage))
                return
            val = [val_old[0], False, ""]
            model.set_val_by_link(slink, val)
        self.dialog_info("Constraints are deleted")
        self.get_link_on_constr(model)

if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel.cmodel_exp()
    mainwind1 = cbuilder(model)

    sys.exit(app.exec_())