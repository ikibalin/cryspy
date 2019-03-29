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
    def __init__(self, fdir = None, ffile = None, data_type = "powder"):
        super(cwidg_central, self).__init__()
        self.data_type = data_type
        self.fdir = fdir
        self.ffile = ffile
        self.init_layout_central()
        self.setLayout(self.layout_central)
        if (not((ffile == None)|(fdir == None)|(ffile == ""))):
            if os.path.exists(os.path.join(fdir,ffile)):
                self.open_file()

    def init_layout_central(self):
        """
        make central layout
        """
        lay_central = QtWidgets.QVBoxLayout()

        data_type = self.data_type
        
        if data_type != "bkgr":
            lay_h = QtWidgets.QHBoxLayout()
            l_wl = QtWidgets.QLabel("wavelength (in angstrems):")
            lay_h.addWidget(l_wl)
            lay_h.addStretch(1)
            self.le_wl = QtWidgets.QLineEdit("")
            lay_h.addWidget(self.le_wl)
            lay_central.addLayout(lay_h)

            lay_h = QtWidgets.QHBoxLayout()
            l_field = QtWidgets.QLabel("magnetic field (in teslas):")
            lay_h.addWidget(l_field)
            lay_h.addStretch(1)
            self.le_field = QtWidgets.QLineEdit("")
            lay_h.addWidget(self.le_field)
            lay_central.addLayout(lay_h)
            
        if data_type == "mono":
            l_u = QtWidgets.QLabel("orientation matrix (U):")
            lay_central.addWidget(l_u)
            self.te_u = QtWidgets.QTextEdit()
            self.te_u.setPlainText(" 1.000    0.000    0.000\n 0.000    1.000    0.000\n 0.000    0.000    1.000\n")
            lay_central.addWidget(self.te_u)

        if data_type == "powder":
            l_refl = QtWidgets.QLabel("\n\nSeparated by spaces: \nttheta (in degrees), intensity_up, sigma, intensity_down, sigma\n")
            lay_central.addWidget(l_refl)

            lay_h = QtWidgets.QHBoxLayout()
            l_sim = QtWidgets.QLabel("(in case of simmulation just write first diffraction angle, step and final one through spaces)")
            lay_h.addWidget(l_sim)
            le_sim = QtWidgets.QLabel("")
            lay_h.addWidget(le_sim)
            lay_central.addLayout(lay_h)
        elif data_type == "mono":
            l_refl = QtWidgets.QLabel("\n\nSeparated by spaces: \n h, k, l, flip ratio, sigma\n")
            lay_central.addWidget(l_refl)

            lay_h = QtWidgets.QHBoxLayout()
            l_sim = QtWidgets.QLabel("(in case of simmulation write h_min hmax k_min k_max l_min l_max through spaces)")
            lay_h.addWidget(l_sim)
            le_sim = QtWidgets.QLabel("")
            lay_h.addWidget(le_sim)
            lay_central.addLayout(lay_h)
        elif data_type == "bkgr":
            l_refl = QtWidgets.QLabel("\n\nSeparated by spaces: \n ttheta (in degrees), intensity\n")
            lay_central.addWidget(l_refl)

        self.te_refl = QtWidgets.QTextEdit()
        lay_central.addWidget(self.te_refl)
        
        if data_type == "powder":
            self.cb_sigma = QtWidgets.QCheckBox("calculate sigmas as root of intensity (if true, sigmas are skipped)")
            lay_central.addWidget(self.cb_sigma)
            self.cb_down = QtWidgets.QCheckBox("intensity down is equal to intensity up (if true, intensity_down and its sigma are skipped)")
            lay_central.addWidget(self.cb_down)

            lay_central.addLayout(lay_h)



        lay_h = QtWidgets.QHBoxLayout()
        #b_open = QtWidgets.QPushButton("open")
        #b_open.clicked.connect(self.open_new_file)
        #lay_h.addWidget(b_open)
        b_save = QtWidgets.QPushButton("save")
        b_save.clicked.connect(self.create_file)
        lay_h.addWidget(b_save)
        lay_h.addStretch(1)
        lay_central.addLayout(lay_h)

        self.layout_central = lay_central

    def open_new_file(self):
        fdir = self.fdir
        if fdir == None:
            fdir = os.getcwd()
        lhelp, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdir,
                                                      "Extension file (*.dat)")
        if (not ok):
            return
        fdir = os.path.dirname("{}".format(lhelp))
        fdir = fdir.replace("/","\\")
        ffile = os.path.basename("{}".format(lhelp))
        self.ffile = ffile
        self.fdir = fdir
        self.open_file()
        
    def open_file(self):
        fdir = self.fdir
        ffile = self.ffile

        data_type = self.data_type
        ddata = cfunc.readexpdata(os.path.join(fdir, ffile))
        
        if data_type == "powder": 
            self.le_wl.setText("{:}".format(ddata["wavelength"]))
            self.le_field.setText(" {0[2]:}".format(ddata["field"]))
            lang = ddata["ttheta"]
            lint_u = ddata["IntUP"]
            lsint_u = ddata["sIntUP"]
            lint_d = ddata["IntDOWN"]
            lsint_d = ddata["sIntDOWN"]
            ltext=[]
            for ang, int_u, sint_u, int_d, sint_d in zip(lang, lint_u, lsint_u, lint_d, lsint_d):
                ltext.append(" {:7.2f} {:15.5f} {:15.5f} {:15.5f} {:15.5f}".format(ang, int_u, sint_u, int_d, sint_d))
        elif data_type == "mono": 
            self.le_wl.setText("{:}".format(ddata["wavelength"]))
            self.le_field.setText(" {0[2]:}".format(ddata["field"]))
            mu = [ddata["orientation"][0:3], ddata["orientation"][3:6], ddata["orientation"][6:9]]
            self.te_u.setPlainText("\n".join([" {0[0]:} {0[1]:} {0[2]:}".format(hh) for hh in mu]))
            lh = ddata["h"]
            lk = ddata["k"]
            ll = ddata["l"]
            lfr = ddata["FR"]
            lsfr = ddata["sFR"]
            ltext=[]
            for h, k, l, fr, sfr in zip(lh, lk, ll, lfr, lsfr):
                ltext.append(" {:4} {:4} {:4} {:10.5f} {:10.5f}".format(h, k, l, fr, sfr))
        elif data_type == "bkgr": 
            lttheta = ddata["ttheta"]
            lbkgr = ddata["IntBKGR"]
            ltext=[]
            for ttheta, bkgr in zip(lttheta, lbkgr):
                ltext.append(" {:7.2f} {:15.5f}".format(ttheta, bkgr))
        self.te_refl.setPlainText("\n".join(ltext))

    def create_file(self):
        data_type = self.data_type
        if data_type == "powder":
            lsout = self.create_file_powder()
        elif data_type == "mono":
            lsout = self.create_file_mono()
        elif data_type == "bkgr":
            lsout = self.create_file_bkgr()
        else:
            print "Unknwon type of data"
            return
        ffile = self.ffile
        fdir = self.fdir
        if fdir == None:
            fdir = os.getcwd()
        if ffile == "":
            lhelp, ok = QtWidgets.QFileDialog.getSaveFileName(self, 'Save file', fdir,
                                                      "Extension file (*.dat)")
            if (not ok):
                return
            fdir = os.path.dirname("{}".format(lhelp))
            fdir = fdir.replace("/","\\")
            ffile = os.path.basename("{}".format(lhelp))
        self.ffile = ffile
        self.fdir = fdir

        ffull = os.path.join(fdir, ffile)
        fid = open(ffull,'w')
        fid.write("\n".join(lsout))
        fid.close()
        QtWidgets.QMessageBox.question(self, "Info", "File has been saved.")

    def create_file_powder(self):
        lsout=[]
        lsout.append("#wavelength {}".format(self.le_wl.text()))
        lsout.append("#field 0. 0. {}".format(self.le_field.text()))
        lsout.append("#   ttheta     IntUP    sIntUP   IntDOWN  sIntDOWN")
        ltext = ["{:}".format(hh.strip()) for hh in self.te_refl.toPlainText().split("\n") if hh.strip() != ""]
        if ((len(ltext) == 1)&(len(ltext[0].split()) == 3)):
            #simmulation mode
            [sangle1, sstep, sangle2] = ltext[0].split()
            angle1, step, angle2 = float(sangle1), float(sstep), float(sangle2)
            npoint = int((angle2 - angle1)//step)+1
            langle = ["{:7.2f}".format(angle1+step*ind) for ind in range(npoint)]
            lint_u = ["0.00000" for hh in langle]
            lsint_u = ["1.00000" for hh in langle]
            lint_d = ["0.00000" for hh in langle]
            lsint_d = ["1.00000" for hh in langle]
        else:
            langle, lint_u, lsint_u, lint_d, lsint_d = [], [], [], [], []
            ind_tth, ind_int_u, ind_sint_u, ind_int_d, ind_sint_d  = 0, 1, 2, 3, 4
            if (self.cb_sigma.checkState()!=0):
                ind_sint_u, ind_sint_d = -1, -1
                if (self.cb_down.checkState()!=0):
                    ind_int_d = ind_int_u
                else:
                    ind_int_d = ind_int_u + 1
            else:
                if (self.cb_down.checkState()!=0):
                    ind_int_d = ind_int_u
                    ind_sint_d = ind_sint_u
            for line in ltext:
                langle.append(line.split()[ind_tth])
                lint_u.append(line.split()[ind_int_u])
                lint_d.append(line.split()[ind_int_d])
                if ind_sint_u == -1:
                    lsint_u.append(" {:.3f} ".format((float(line.split()[ind_int_u]))**0.5))
                else:
                    lsint_u.append(line.split()[ind_sint_u])
                if ind_sint_d == -1:
                    lsint_d.append(" {:.3f} ".format((float(line.split()[ind_int_d]))**0.5))
                else:
                    lsint_d.append(line.split()[ind_sint_d])
        for angle, int_u, sint_u, int_d, sint_d in zip(langle, lint_u, lsint_u, lint_d, lsint_d):
            lsout.append(" {:}   {:}   {:}   {:}   {:}".format(angle, int_u, sint_u, int_d, sint_d))
        return lsout
    
    def create_file_mono(self):
        lsout=[]
        lsout.append("#wavelength {}".format(self.le_wl.text()))
        lsout.append("#field 0. 0. {}".format(self.le_field.text()))
        lsout.append("#orientation {}".format(" ".join(self.te_u.toPlainText().split("\n"))))
        lsout.append("#     h      k       l      FR     sFR")
        ltext = ["{:}".format(hh.strip()) for hh in self.te_refl.toPlainText().split("\n") if hh.strip() != ""]
        lh, lk, ll, lfr, lsfr = [], [], [], [], []
        if len(ltext) == 1:
            [h_min, h_max, k_min, k_max, l_min, l_max] = ltext[0].split()[:6]
            for ih in range(int(h_min),int(h_max)+1):
                for ik in range(int(k_min),int(k_max)+1):
                    for il in range(int(l_min),int(l_max)+1):
                        lh.append(" {:4}".format(ih))
                        lk.append(" {:4}".format(ik))
                        ll.append(" {:4}".format(il))
            lfr = len(lh)*[ "1.00000" ]
            lsfr = len(lh)*[ "0.10000" ]
        for line in ltext:
            [h, k, l, fr, sfr] = line.split()[:5]
            lh.append(h)
            lk.append(k)
            ll.append(l)
            lfr.append(fr)
            lsfr.append(sfr)
        for h, k, l, fr, sf in zip(lh, lk, ll, lfr, lsfr):
            lsout.append(" {:} {:} {:} {:} {:}".format(h, k, l, fr, sf))
        return lsout

    def create_file_bkgr(self):
        lsout=[]
        lsout.append("#   ttheta     IntBKGR")
        ltext = ["{:}".format(hh.strip()) for hh in self.te_refl.toPlainText().split("\n") if hh.strip() != ""]
        for line in ltext:
            ttheta = line.split()[0]
            int_bkgr = line.split()[1]
            lsout.append(" {:}   {:}".format(ttheta, int_bkgr))
        return lsout

if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    mainwind1 = cwind_central()

    sys.exit(app.exec_())
