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
from PyQt5 import QtGui
from PyQt5 import QtCore

def del_layout(layout):
    """
    delete all elements from layouts
    """
    for i in reversed(range(layout.count())):
        if layout.itemAt(i).widget() != None:
           layout.itemAt(i).widget().setParent(None)
        elif layout.itemAt(i).layout() != None:
            del_layout(layout.itemAt(i).layout())
            layout.itemAt(i).layout().setParent(None)
        else:
            layout.removeItem(layout.itemAt(i))
    return

class cmodel_min(object):
    """
    minimal model for cwidget
    """
    def __init__(self):
        super(cmodel_min, self).__init__()
        self.obs_lab_o_par_o_lab_m_par_m = []
        self.builder = []

    def set_val(self, lab_m, par_m, val, par_val):
        """
        par_m can be None
        par_val can be None
        """
        #it is temporary solution, it is not good for logical values
        try:
            val_new = float(val)
        except:
            val_new = val
        if par_m == None:
            self.__dict__[lab_m] = val_new
        else:
            self.__dict__[lab_m][par_m] = val_new
        self.val_changed(lab_m, par_m, val, par_val)

    def add_obs_labs_pars(self, obs, lab_o, par_o, lab_m, par_m):
        """
        par_o, par_m can be None
        """
        if (not ((obs, lab_o, par_o, lab_m, par_m) in self.obs_lab_o_par_o_lab_m_par_m)):
            self.obs_lab_o_par_o_lab_m_par_m.append((obs, lab_o, par_o, lab_m, par_m))

    def del_obs_labs_pars(self, obs, lab_o, par_o, lab_m, par_m):
        """
        par_o, par_m can be None
        """
        self.obs_lab_o_par_o_lab_m_par_m.remove((obs, lab_o, par_o, lab_m, par_m))

    def del_obs(self, obs):
        """
        delete observer
        """
        l_to_remove = [hh for hh in self.obs_lab_o_par_o_lab_m_par_m if hh[0] == obs]
        for to_remove in l_to_remove:
            self.obs_lab_o_par_o_lab_m_par_m.remove(to_remove)

    def val_changed(self, lab_m, par_m, val, par_val):
        for obs_lab_o_par_o_lab_m_par_m in self.obs_lab_o_par_o_lab_m_par_m:
            (obs, l_o, p_o, l_m, p_m) = obs_lab_o_par_o_lab_m_par_m
            if (l_m, p_m) == (lab_m, par_m):
                val_old = obs.get_val(l_o, p_o, par_val)
                if val_old != val:
                    obs.set_val(l_o, p_o, val, par_val)

    def mod_changed(self):
        for build in self.builder:
            self.del_obs(build)
            build.init_widget(self)

    def add_builder(self, build):
        if (not (build in self.builder)):
            self.builder.append(build)

    def del_builder(self, build):
        l_to_remove = [hh for hh in self.builder if hh == build]
        for to_remove in l_to_remove:
            self.builder.remove(to_remove)

    #the changing of the model should be inserted
    def set_val_by_link(self, slink, val):
        """take value of parameter from object
        obj is object
        slink is a link on parameter (attribute of object,name of value)(attribute of object)
        val -value to change. If None than just output of the value over the link without replacement
        el_obj should have a name attribute by defifnition
        """
        val_new, lmessage = None, []
        sname_el_obj = "name"
    
        num1, num2 = slink.find("("), slink.find(")")
        sh1 = slink[(num1+1):num2]
        lhelp = sh1.split(",")
        sname_atr_obj = lhelp[0]
        lmessage_new = []
        if (sname_atr_obj in self.__dict__.keys()):
            if (len(lhelp)==2):
                sval_el_obj = lhelp[1]
                for el_obj in self.__dict__[sname_atr_obj]:
                    if (sname_el_obj in el_obj.__dict__.keys()):
                        if el_obj.__dict__[sname_el_obj] == sval_el_obj:
                            val_new, lmessage_new = el_obj.set_val_by_link(slink[(num2+1):], val)
                            lmessage.extend(lmessage_new)
                            break
                    else:
                        lmessage.append("there is no parameter '{:}' in the object".format(sname_el_obj))
                if (val_new == None)&(lmessage_new == []):
                    lmessage.append("the name '{:}' for '{:}' not found".format(sval_el_obj, sname_atr_obj))
            elif (len(lhelp)==1):
                val_old = self.__dict__[sname_atr_obj]
                if (val == None):
                    val_new = val_old
                elif (not(isinstance(val, float))):
                    val_new = [hh for hh in val]
                    self.__dict__[sname_atr_obj] = val_new
                else:
                    val_new = [val, val_old[1], val_old[2]]
                    self.__dict__[sname_atr_obj] = val_new
            else:
                lmessage.append("unknown link {:} for the object".format(slink))
        else:
            lmessage.append("there is no parameter '{:}' in the object".format(sname_atr_obj))
        return val_new, lmessage


class cwidget_min(QtWidgets.QWidget):
    """
    minimal widget which are compatible with model
    """
    def __init__(self, model):
        super(cwidget_min, self).__init__()
        self.init_widget(model)
        self.setLayout(self.layout_central)

    def set_val(self, lab_o, par_o, val, par_val):
        """
        par_val can be None

        For all:
            if par_val is string, then val = val.__dict__[par_val] if val is object and
            val = [hh.__dict__[par_val] for hh in val] if val is object

        For TableWidget
        par_val = None - change text
        par_val = 1 - change checkbox in tableitem
        par_val = 2 - change editability in tableitem
        """
        if isinstance(par_val, str):
            if isinstance(val, list):
                val = [hh.__dict__[par_val] for hh in val]
            else:
                val = val.__dict__[par_val]

        if par_o == None:
            if (isinstance(self.__dict__[lab_o], QtWidgets.QLineEdit)|
                isinstance(self.__dict__[lab_o], QtWidgets.QLabel)):
                if not(isinstance(val, str)):
                    try:
                        val = "{:.5f}".format(val)
                    except:
                        val = str(val)
                self.__dict__[lab_o].setText(val)
            elif isinstance(self.__dict__[lab_o], QtWidgets.QCheckBox):
                self.__dict__[lab_o].setCheckState(2*int(val))
            elif isinstance(self.__dict__[lab_o], QtWidgets.QComboBox):
                #here val is list
                self.__dict__[lab_o].addItems(val)
            else:
                print "Class '{:}' is not defined for cwidget_min".format(type(self.__dict__[lab_o]))
                self.__dict__[lab_o] = val
        else:
            if isinstance(self.__dict__[lab_o], QtWidgets.QComboBox):
                #c_ind = self.__dict__[lab_o].currentIndex()
                self.__dict__[lab_o].setItemText(par_o, val)
            elif isinstance(self.__dict__[lab_o], QtWidgets.QTableWidget):
                flag = False
                item_old = self.__dict__[lab_o].item(par_o[0], par_o[1])
                if item_old == None:
                    flag = True
                    item_old = QtWidgets.QTableWidgetItem()
                    
                if par_val == None:
                    if not(isinstance(val, str)):
                        try:
                            val = "{:.5f}".format(val)
                        except:
                            val = str(val)
                    item_old.setText(val)
                elif par_val == 1:
                    #refinement
                    item_old.setCheckState(2*int(val))
                elif par_val == 2:
                    #constraints
                    if val != "":
                        qcolor = QtGui.QColor(255,100,0,50)
                        qbrush = QtGui.QBrush(qcolor) 
                        item_old.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled)
                        item_old.setBackground(qbrush)
                    else:
                        qcolor = QtGui.QColor(255,255,255,255)
                        qbrush = QtGui.QBrush(qcolor) 
                        item_old.setFlags(QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsEditable|QtCore.Qt.ItemIsUserCheckable)
                        item_old.setBackground(qbrush)
                else:
                    print "parameter '{:}' for table widget is unknown".format(par_val)
                if flag:
                    self.__dict__[lab_o].setItem(par_o[0], par_o[1], item_old)
                    
    def get_val(self, lab_o, par_o, par_val):
        """
        take val from the widget
        par_val can be None

        For all:
            if par_val is string, then val = val.__dict__[par_val] if val is object and
            val = [hh.__dict__[par_val] for hh in val] if val is object

        For TableWidget
        par_val = None - change text
        par_val = 1 - change checkbox in tableitem
        par_val = 2 - change editability in tableitem
        """

        if par_o == None:
            if (isinstance(self.__dict__[lab_o], QtWidgets.QLineEdit)|
                isinstance(self.__dict__[lab_o], QtWidgets.QLabel)):
                val = "{:}".format(self.__dict__[lab_o].text())
            elif isinstance(self.__dict__[lab_o], QtWidgets.QCheckBox):
                val  = (self.__dict__[lab_o].checkState() != 0)
            elif isinstance(self.__dict__[lab_o], QtWidgets.QComboBox):
                #here val is list
                val = ["{:}".format(self.__dict__[lab_o].itemText(hh)) for hh in range(self.__dict__[lab_o].count())]
            else:
                print "Class '{:}' is not defined for cwidget_min".format(type(self.__dict__[lab_o]))
                val = self.__dict__[lab_o]
        else:
            if isinstance(self.__dict__[lab_o], QtWidgets.QComboBox):
                #c_ind = self.__dict__[lab_o].currentIndex()
                val = self.__dict__[lab_o].itemText(par_o)
            elif isinstance(self.__dict__[lab_o], QtWidgets.QTableWidget):
                item_old = self.__dict__[lab_o].item(par_o[0], par_o[1])
                if item_old == None:
                    val = None
                    return val
                if par_val == None:
                    val = "{:}".format(item_old.text())
                elif par_val == 1:
                    #refinement
                    val = (item_old.checkState() != 0)
                elif par_val == 2:
                    #constraints
                    val = None
                else:
                    print "parameter '{:}' for table widget is unknown".format(par_val)
        return val
    
    def init_widget(self, model):
        lay_central = self.layout()
        if lay_central == None:
            lay_central = QtWidgets.QVBoxLayout()
        else:
            del_layout(lay_central)
        self.layout_central = lay_central
        
    def put_le(self, lay, lab_o, model, lab_m, val_type):
        """
        put line edit for variable of model such as value (line edit + checkbox) or text (line edit)
        val_type = 'text' or 'val'
        """
        lay_h = QtWidgets.QHBoxLayout()
        self.__dict__[lab_o] = QtWidgets.QLineEdit()
        if val_type == 'text':
            val = "{:}".format(model.__dict__[lab_m])
            self.set_val(lab_o, None, val, None)
            par_m = None
            self.__dict__[lab_o].editingFinished.connect(lambda:
                model.set_val(lab_m, par_m, self.__dict__[lab_o].text(), None))
            self.__dict__[lab_o].setAlignment(QtCore.Qt.AlignLeft)
        elif val_type == 'val':
            val = model.__dict__[lab_m]
            self.set_val(lab_o, None, val[0], None)
            
            if val[2] != "":
                self.__dict__[lab_o].setReadOnly(True)
                self.__dict__[lab_o].setStyleSheet("color: rgb(156, 156, 156);")
            else:
                self.__dict__[lab_o].setReadOnly(False)
                self.__dict__[lab_o].setStyleSheet("color: rgb(0, 0, 0);")

            lab_o_cb = lab_o + "_cb"
            self.__dict__[lab_o_cb] = QtWidgets.QCheckBox()
            self.set_val(lab_o_cb, None, val[1], None)
            par_m = 0
            self.__dict__[lab_o].editingFinished.connect(lambda:
                model.set_val(lab_m, par_m, self.__dict__[lab_o].text(), None))
            self.__dict__[lab_o].setAlignment(QtCore.Qt.AlignRight)
            
        self.__dict__[lab_o].setFrame(False)
        self.__dict__[lab_o].setTextMargins(0, 0, 0, 0)
        model.add_obs_labs_pars(self, lab_o, None, lab_m, par_m)
    
            
        lay_h.addWidget(self.__dict__[lab_o])
        if val_type == 'val':
            model.add_obs_labs_pars(self, lab_o_cb, None, lab_m, 1)
            self.__dict__[lab_o_cb].stateChanged.connect(lambda:
                model.set_val(lab_m, 1, self.__dict__[lab_o_cb].checkState() != 0, None))
            lay_h.addWidget(self.__dict__[lab_o_cb])
            
        lay.addLayout(lay_h)
        
    def put_checkb(self, lay, lab_o, model, lab_m, label):
        """
        model.__dict__[lab_m] should be logical
        """
        self.__dict__[lab_o] = QtWidgets.QCheckBox(label)
        self.set_val(lab_o, None, model.__dict__[lab_m], None)
        model.add_obs_labs_pars(self, lab_o, None, lab_m, None)
        self.__dict__[lab_o].stateChanged.connect(lambda n1: model.set_val(lab_m, None, n1!=0, None))
        lay.addWidget(self.__dict__[lab_o])


    def put_tab(self, lay, lab_o, lmodel, llab_m, lval_type, lheader):
        """
        put table for variable of model such as value (table item with checkbox) or text (table item)
        val_type = 'text' or 'val'
        lmodels is a list of model
        llab_m is a list of lab_m which are in the model: model.lab_m
        lval_type is a list of val_type for each model.lab_m
        
        par_val = None    ==>  val is text
        par_val = 1    ==>  val is for checkbox
        par_val = 2    ==>  val is for constraints
        """        
        n_row = len(lmodel)
        n_count = len(llab_m)
        self.__dict__[lab_o] = QtWidgets.QTableWidget(n_row, n_count)
        """
        self.__dict__[lab_o].setRowCount(n_row)
        self.__dict__[lab_o].setColumnCount(n_count)
        """
        if lheader!=[]:
            self.__dict__[lab_o].setHorizontalHeaderLabels(lheader)
            self.__dict__[lab_o].horizontalHeader().setStretchLastSection(True)
            
        for iatom, atom in enumerate(lmodel):
            ilab_m = -1
            for lab_m, val_type in zip(llab_m, lval_type):
                ilab_m += 1
                par_o = (iatom, ilab_m)
                val = atom.__dict__[lab_m]
                if val_type == "text":
                    self.set_val(lab_o, par_o, val, None)
                    atom.add_obs_labs_pars(self, lab_o, par_o, lab_m, None)
                elif val_type == "val":
                    self.set_val(lab_o, par_o, val[0], None)
                    self.set_val(lab_o, par_o, val[1], 1)
                    self.set_val(lab_o, par_o, val[2], 2)
                    atom.add_obs_labs_pars(self, lab_o, par_o, lab_m, 0)
                    atom.add_obs_labs_pars(self, lab_o, par_o, lab_m, 1)
                    atom.add_obs_labs_pars(self, lab_o, par_o, lab_m, 2)
        lay.addWidget(self.__dict__[lab_o])
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.__dict__[lab_o].setSizePolicy(sizePolicy)
        self.__dict__[lab_o].resizeColumnsToContents()
        #actually the used signal is bad because it is emitted  when cell is changed, 
        #but it should be emitted when user changes the values
        self.__dict__[lab_o].cellChanged.connect(lambda nr, nc: self.signal_for_tab(lab_o, (nr,nc), lmodel, llab_m, lval_type))

    def signal_for_tab(self, lab_o, par_o, lmodel, llab_m, lval_type):
        """
        It is additional function to self.put_tab
        """
        n_at = par_o[0]
        n_par = par_o[1]
        model = lmodel[n_at]
        lab_m = llab_m[n_par]
        val_type = lval_type[n_par]
        if val_type == 'text':
            par_m, par_val = None, None
            val = self.get_val(lab_o, par_o, par_val)
            model.set_val(lab_m, par_m, val, par_val)
        elif val_type == 'val':
            par_m, par_val = 0, None
            val = self.get_val(lab_o, par_o, par_val)
            model.set_val(lab_m, par_m, val, par_val)
        
            par_m, par_val = 1, 1
            val = self.get_val(lab_o, par_o, par_val)
            model.set_val(lab_m, par_m, val, par_val)

    def dialog_ask(self,question):
        text, ok = QtWidgets.QInputDialog.getText(self, 'Question',question)
        return text, ok
    
    def dialog_info(self,info):
        QtWidgets.QMessageBox.information(self, "Information", info)
        

class cbuilder_min(QtWidgets.QMainWindow):
    """
    minimal builder
    """
    def __init__(self, model):
        super(cbuilder_min, self).__init__()
        self.init_widget(model)
        self.show()
    def init_widget(self, model):
        widg_central = cwidget_min(model)
        self.setCentralWidget(widg_central)


if __name__ == '__main__':
    larg = sys.argv
    app = QtWidgets.QApplication(larg)
    model = cmodel_min()
    mainwind = cbuilder_min(model)

    sys.exit(app.exec_())