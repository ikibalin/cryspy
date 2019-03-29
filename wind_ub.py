# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 09:52:24 2019

@author: ikibalin
"""

import os
import sys
import traceback

import xml
import xml.etree.ElementTree
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore

#it is not correct in general case when start from any folder !!! SHOULD BE CORRECTED
larg = sys.argv
if len(larg) >= 2:
    #script mode
    fdirprog = os.path.dirname(larg[0])
else:
    fdirprog = os.getcwd()

sys.path.insert(0,os.path.join(fdirprog,'code_core_agent'))

import cfunc

class cwind_central(QtWidgets.QMainWindow):
    def __init__(self):
        super(cwind_central, self).__init__()
        widg_central = cwidg_central()
        #self.init_layout_central()
        #widget_central = QtWidgets.QWidget()
        #widget_central.setLayout(self.layout_central)
        self.setCentralWidget(widg_central)
        self.show()


class cwidg_central(QtWidgets.QWidget):
    def __init__(self):
        super(cwidg_central, self).__init__()
        self.init_layout_central()
        self.setLayout(self.layout_central)
        
    def init_layout_central(self):
        """
        make central layout
        """
        lay_central = QtWidgets.QVBoxLayout()
        
        l_ub = QtWidgets.QLabel("UB matrix:")
        lay_central.addWidget(l_ub)
        
        self.te_ub = QtWidgets.QTextEdit()
        lay_central.addWidget(self.te_ub)

        self.cb_transp = QtWidgets.QCheckBox("transposed")
        lay_central.addWidget(self.cb_transp)

        lay_h = QtWidgets.QHBoxLayout()
        lay_h.addStretch(1)
        b_calc = QtWidgets.QPushButton("calculate U and B matrices")
        b_calc.clicked.connect(self.calc_ub)
        lay_h.addWidget(b_calc)
        lay_central.addLayout(lay_h)

        l_u = QtWidgets.QLabel("U matrix (output):")
        lay_central.addWidget(l_u)
        
        self.te_u = QtWidgets.QTextEdit()
        self.te_u.setReadOnly(True)
        lay_central.addWidget(self.te_u)

        l_b = QtWidgets.QLabel("B matrix (output):")
        lay_central.addWidget(l_b)
        
        self.te_b = QtWidgets.QTextEdit()
        self.te_b.setReadOnly(True)
        lay_central.addWidget(self.te_b)

        l_cell = QtWidgets.QLabel("unit cell (output):")
        lay_central.addWidget(l_cell)
        
        lay_h = QtWidgets.QHBoxLayout()
        self.le_a = QtWidgets.QLineEdit()
        self.le_a.setReadOnly(True)
        lay_h.addWidget(self.le_a)
        
        self.le_b = QtWidgets.QLineEdit()        
        self.le_b.setReadOnly(True)
        lay_h.addWidget(self.le_b)
        
        self.le_c = QtWidgets.QLineEdit()        
        self.le_c.setReadOnly(True)
        lay_h.addWidget(self.le_c)
        
        lay_central.addLayout(lay_h)
        

        lay_h = QtWidgets.QHBoxLayout()

        self.le_al = QtWidgets.QLineEdit()
        self.le_al.setReadOnly(True)
        lay_h.addWidget(self.le_al)
        
        self.le_be = QtWidgets.QLineEdit()
        self.le_be.setReadOnly(True)
        lay_h.addWidget(self.le_be)
        
        self.le_ga = QtWidgets.QLineEdit()
        self.le_ga.setReadOnly(True)
        lay_h.addWidget(self.le_ga)
        
        lay_central.addLayout(lay_h)
        self.layout_central = lay_central
        
    def calc_ub(self):
        s_ub = str(self.te_ub.toPlainText())
        lltext = [hh.split() for hh in s_ub.split("\n")]
        lflags = [len(hh) != 3 for hh in lltext]
        if ((len(lltext) != 3)&(any(lflags))):
            print "Text is not correspond dimension matrix 3x3"
            return
        try:
            m_ub_input = [[float(hh1) for hh1 in hh] for hh in lltext]
        except:
            print "Can not convert text to float"
            return
        if self.cb_transp.checkState() == 2:
            m_ub = cfunc.transposeM(m_ub_input)
        else:
            m_ub = m_ub_input
        m_u, m_b, ucp = cfunc.calc_u_b_ucp_from_ub(m_ub)
        text_u = "\n".join([" ".join(["{:12.7f}".format(hh1) for hh1 in hh]) for hh in m_u])
        text_b = "\n".join([" ".join(["{:12.7f}".format(hh1) for hh1 in hh]) for hh in m_b])
        self.te_u.setPlainText(text_u)
        self.te_b.setPlainText(text_b)
        self.le_a.setText("{:.7f}".format(ucp[0]))
        self.le_b.setText("{:.7f}".format(ucp[1]))
        self.le_c.setText("{:.7f}".format(ucp[2]))
        self.le_al.setText("{:.3f}".format(ucp[3]))
        self.le_be.setText("{:.3f}".format(ucp[4]))
        self.le_ga.setText("{:.3f}".format(ucp[5]))
        
    
if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    mainwind1 = cwind_central()

    sys.exit(app.exec_())
