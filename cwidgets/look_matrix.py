import os
import sys
 
from PyQt5 import QtWidgets
from PyQt5 import QtGui
from PyQt5 import QtCore
 
import matplotlib
import matplotlib.backends.backend_qt5agg
import matplotlib.figure 
import matplotlib.pyplot

import numpy
import math
import scipy
import scipy.optimize

class ccomunicate(QtCore.QObject):
    sig_m = QtCore.pyqtSignal()
    sig_l_tth = QtCore.pyqtSignal()
    sig_l_phi = QtCore.pyqtSignal()
    sig_l_redraw = QtCore.pyqtSignal()


def func_g(x, pos, width, iint):
    val_1 = iint/(width*(2.*math.pi)**0.5)
    val_2 = numpy.exp(-0.5*((x-pos)/width)**2)
    return val_1*val_2


def func_g_2d(np_tth, np_phi, p_iint, p_tth0, p_sig):
    """
    calc gauss function in 2D space:
        p_iint: iint = p1 + p2*sin(phi)**2 + p3*sin(phi)**4 
        p_tth: tth = p1 + p2*phi + p3*phi**2 
        p_sig: sig = p1 + p2*phi + p3*phi**2
    """
    np_tth_rad, np_phi_rad = np_tth*math.pi/180., np_phi*math.pi/180.
    np_tth_2d, np_phi_2d = numpy.meshgrid(np_tth, np_phi_rad)
    np_iint_2d = p_iint[0] + p_iint[1] * numpy.cos(np_phi_2d)**2 + p_iint[2] * numpy.cos(np_phi_2d)**4
    np_tth0_2d = p_tth0[0] + p_tth0[1] * np_phi_2d + p_tth0[2] * np_phi_2d**2
    np_sig_2d = p_sig[0] + p_sig[1] * np_phi_2d + p_sig[2] * np_phi_2d**2
    
    res_2d = func_g(np_tth_2d, np_tth0_2d, np_sig_2d, np_iint_2d)
    return res_2d


def det_2(m11, m12, m21, m22):
    return m11*m22 - m12*m21

def plane_by_3_points(x, y, x1, y1, z1, x2, y2, z2, x3, y3, z3):
    d_x = det_2(y2-y1, z2-z1, y3-y1, z3-z1)
    d_y = det_2(x2-x1, z2-z1, x3-x1, z3-z1)
    d_z = det_2(x2-x1, y2-y1, x3-x1, y3-y1)
    z = z1 + (- d_x*(x-x1) + d_y*(y-y1))/d_z
    return z

xyz1, xyz2, xyz3 = (1,2,3), (2,7,0), (6,-7,4) 



def func_gausses_2d_bkg(np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3):
    """
    
    """
    np_g_2d = numpy.zeros((np_phi.size, np_tth.size), dtype = float)
    for p_iint, p_tth0, p_sig in zip(lp_iint, lp_tth0, lp_sig):
        np_g_2d += func_g_2d(np_tth, np_phi, p_iint, p_tth0, p_sig)
        
    np_tth_2d, np_phi_2d = numpy.meshgrid(np_tth, np_phi)
    np_b = plane_by_3_points(np_tth_2d, np_phi_2d, *(b1+b2+b3))            
    res = np_g_2d + np_b
    return res




def refine_by_gauss(inp_x, inp_y, inp_sy, int_0_x, int_0_y, bl_x, bl_0_y,
                    br_x, br_0_y):
    width_0 = abs(br_x-bl_x)*0.2
    iint_0 = (int_0_y - 0.5*(bl_0_y+br_0_y))*width_0*(2.*math.pi)**0.5
    def func_gl(x, pos, width, iint, b1_x, b1_y, b2_x, b2_y):
        val1 = func_g(x, pos, width, iint)
        val2 = (b2_y-b1_y)*(x-b1_x)/(b2_x-b1_x)+b1_y
        return val1 + val2
    
    if bl_x < br_x:
        b1_x = bl_x 
        b1_y = bl_0_y 
        b2_x = br_x 
        b2_y = br_0_y 
    else:
        b2_x = bl_x 
        b2_y = bl_0_y
        b1_x = br_x 
        b1_y = br_0_y 
        

    flag_1 = inp_x >= b1_x
    flag_2 = inp_x <= b2_x
    
    np_x = inp_x[numpy.logical_and(flag_1, flag_2)]
    np_y_e = inp_y[numpy.logical_and(flag_1, flag_2)]
    np_sy_e = inp_sy[numpy.logical_and(flag_1, flag_2)]
    
    np_y = func_gl(np_x, int_0_x, width_0, iint_0, b1_x, b1_y, b2_x, b2_y)
    lpar_0 = [int_0_x, iint_0, width_0]
    
    def tempfunc(par):
        pos = par[0]
        iint = par[1]
        width = par[2]
        np_y = func_gl(np_x, pos, width, iint, b1_x, b1_y, b2_x, b2_y)
        chi_sq = ((np_y_e-np_y)/np_sy_e)**2
        return chi_sq.sum()
    
    res = scipy.optimize.minimize(tempfunc, lpar_0, method='Nelder-Mead')
    lpar = [float(hh) for hh in res.x]
    pos, iint, width = lpar[0], lpar[1], lpar[2]
    np_y = func_gl(np_x, pos, width, iint, b1_x, b1_y, b2_x, b2_y)
    return np_x, np_y, b1_y, b2_y, pos, iint, width 

class cmodel(QtCore.QObject):
    """
    class to control changings in the widget and to interact with data
    """
    def __init__(self):
        super(cmodel, self).__init__()
        self.data_matrix_e_u = numpy.nan
        self.data_matrix_e_su = numpy.nan
        self.data_matrix_e_d = numpy.nan
        self.data_matrix_e_sd = numpy.nan
        self.data_matrix_e_sum = numpy.nan
        self.data_matrix_e_ssum = numpy.nan
        self.data_matrix_e_diff = numpy.nan
        self.data_matrix_e_sdiff = numpy.nan
        self.data_matrix_m_u = numpy.nan
        self.data_matrix_m_d = numpy.nan
        self.data_matrix_m_sum = numpy.nan
        self.data_matrix_m_diff = numpy.nan
        self.data_range_x = numpy.nan
        self.data_range_y = numpy.nan
        
        self.data_lines_xz_e_u = [None]
        self.data_lines_xz_e_su = [None]
        self.data_lines_xz_e_d = [None]
        self.data_lines_xz_e_sd = [None]
        self.data_lines_xz_e_sum = [None]
        self.data_lines_xz_e_ssum = [None]
        self.data_lines_xz_e_diff = [None]
        self.data_lines_xz_e_sdiff = [None]
        
        self.data_lines_yz_e_u = [None]
        self.data_lines_yz_e_su = [None]
        self.data_lines_yz_e_d = [None]
        self.data_lines_yz_e_sd = [None]
        self.data_lines_yz_e_sum = [None]
        self.data_lines_yz_e_ssum = [None]
        self.data_lines_yz_e_diff = [None]
        self.data_lines_yz_e_sdiff = [None]
        
        self.data_lines_xz_m_u = [None]
        self.data_lines_xz_m_d = [None]
        self.data_lines_xz_m_sum = [None]
        self.data_lines_xz_m_diff = [None]
        
        self.data_lines_yz_m_u = [None]
        self.data_lines_yz_m_d = [None]
        self.data_lines_yz_m_sum = [None]
        self.data_lines_yz_m_diff = [None]
        
        self.min_x = numpy.nan
        self.max_x = numpy.nan
        self.min_y = numpy.nan
        self.max_y = numpy.nan
        self.min_z_u = numpy.nan
        self.max_z_u = numpy.nan
        self.min_z_d = numpy.nan
        self.max_z_d = numpy.nan
        self.min_z_sum = numpy.nan
        self.max_z_sum = numpy.nan
        self.min_z_diff = numpy.nan
        self.max_z_diff = numpy.nan
        self.observer_x = []
        self.observer_y = []
        self.observer_z = []
        
    def read_matrices(self, finp, m_type):
        """
        read matrix
        m_type is 'exp' or 'mod' (experimentally measured (with sigmas) or modelled one (without sigmas))
        """
        ll_int, l_ang1, l_ang2 = [], [], [] 

        fid = open(finp, 'r')
        lcont = fid.readlines()
        fid.close()
               
        lcont = [hh for hh in lcont if (not(hh.startswith("#")))]
        lmat, mat = [], []
        for hh in lcont:
            if hh.strip() == "":
                if mat != []:
                    lmat.append(mat)
                    mat = []
            else:
                mat.append(hh)
        if mat != []:
            lmat.append(mat)
            mat = []
        
        for mat in lmat:
            ang1 = [float(hh) for hh in mat[0].split()[1:]]
            l_int, ang2 = [], []
            for hh in mat[1:]:
                lhelp = hh.split()
                ang2.append(float(lhelp[0]))
                l_int.append([float(hh2) if hh2 != "None" else None for hh2 in lhelp[1:]])
            ll_int.append(l_int)
            l_ang1.append(ang1)
            l_ang2.append(ang2)
        self.data_range_x = numpy.array(l_ang1[0], float)
        self.data_range_y = numpy.array(l_ang2[0], float)

        if m_type == "exp":
            self.data_matrix_e_u = numpy.array(ll_int[0], float)
            self.data_matrix_e_su = numpy.array(ll_int[1], float)
            self.data_matrix_e_d = numpy.array(ll_int[2], float)
            self.data_matrix_e_sd = numpy.array(ll_int[3], float)
            self.data_matrix_e_sum = self.data_matrix_e_u + self.data_matrix_e_d
            self.data_matrix_e_ssum = (self.data_matrix_e_su**2 + self.data_matrix_e_sd**2)**0.5
            self.data_matrix_e_diff = self.data_matrix_e_u - self.data_matrix_e_d
            self.data_matrix_e_sdiff = (self.data_matrix_e_su**2 + self.data_matrix_e_sd**2)**0.5
        elif m_type == "mod":
            self.data_matrix_m_u = numpy.array(ll_int[0], float)
            self.data_matrix_m_d = numpy.array(ll_int[2], float)
            self.data_matrix_m_sum = self.data_matrix_m_u + self.data_matrix_m_d 
            self.data_matrix_m_diff = self.data_matrix_m_u - self.data_matrix_m_d 
        else:
            print "m_type '{:}' is unknown.\nIt should be 'exp' or 'mod'.".format(m_type)
        self.reset_limit()
        
    def reset_limit(self):
        self.min_x = numpy.nanmin(self.data_range_x)
        self.max_x = numpy.nanmax(self.data_range_x)
        self.min_y = numpy.nanmin(self.data_range_y)
        self.max_y = numpy.nanmax(self.data_range_y)        
        
        
        if isinstance(self.data_matrix_e_u, numpy.ndarray):
            self.min_z_u = numpy.nanmin(self.data_matrix_e_u)
            self.max_z_u = numpy.nanmax(self.data_matrix_e_u)
        else:
            self.min_z_u = numpy.nanmin(self.data_matrix_m_u)
            self.max_z_u = numpy.nanmax(self.data_matrix_m_u)
            
        if isinstance(self.data_matrix_e_d, numpy.ndarray):
            self.min_z_d = numpy.nanmin(self.data_matrix_e_d)
            self.max_z_d = numpy.nanmax(self.data_matrix_e_d)
        else:
            self.min_z_d = numpy.nanmin(self.data_matrix_m_d)
            self.max_z_d = numpy.nanmax(self.data_matrix_m_d)
            
        if isinstance(self.data_matrix_e_sum, numpy.ndarray):
            self.min_z_sum = numpy.nanmin(self.data_matrix_e_sum)
            self.max_z_sum = numpy.nanmax(self.data_matrix_e_sum)
        else:
            self.min_z_sum = numpy.nanmin(self.data_matrix_m_sum)
            self.max_z_sum = numpy.nanmax(self.data_matrix_m_sum)
            
        if isinstance(self.data_matrix_e_diff, numpy.ndarray):
            self.min_z_diff = numpy.nanmin(self.data_matrix_e_diff)
            self.max_z_diff = numpy.nanmax(self.data_matrix_e_diff)
        else:
            self.min_z_diff = numpy.nanmin(self.data_matrix_m_diff)
            self.max_z_diff = numpy.nanmax(self.data_matrix_m_diff)

            
    def add_observer_x(self, obs):
        self.observer_x.append(obs)
        
    def add_observer_y(self, obs):
        self.observer_y.append(obs)
        
    def add_observer_z(self, obs):
        self.observer_z.append(obs)
        
    def del_observer_x(self, obs):
        self.observer_x.remove(obs)
        
    def del_observer_y(self, obs):
        self.observer_y.remove(obs)
        
    def del_observer_z(self, obs):
        self.observer_z.remove(obs)
        
    def get_lines_by_index(self, ind_x, ind_y):
        cond_e = isinstance(self.data_matrix_e_u, numpy.ndarray)
        cond_m = isinstance(self.data_matrix_m_u, numpy.ndarray)
        if (cond_e|cond_m):
            cond_x = ind_x != None
            cond_y = ind_y != None
            
        if cond_e:
            if cond_y:
                if (ind_y in range(self.data_matrix_e_u.shape[0])):
                    self.data_lines_xz_e_u[0] = self.data_matrix_e_u[ind_y,:]
                    self.data_lines_xz_e_su[0] = self.data_matrix_e_su[ind_y, :]
                    self.data_lines_xz_e_d[0] = self.data_matrix_e_d[ind_y, :]
                    self.data_lines_xz_e_sd[0] = self.data_matrix_e_sd[ind_y, :]
                    self.data_lines_xz_e_sum[0] = self.data_matrix_e_sum[ind_y, :]
                    self.data_lines_xz_e_ssum[0] = self.data_matrix_e_ssum[ind_y, :]
                    self.data_lines_xz_e_diff[0] = self.data_matrix_e_diff[ind_y, :]
                    self.data_lines_xz_e_sdiff[0] = self.data_matrix_e_sdiff[ind_y, :]
            if cond_x:
                if (ind_x in range(self.data_matrix_e_u.shape[1])):
                    self.data_lines_yz_e_u[0] = self.data_matrix_e_u[:,ind_x]
                    self.data_lines_yz_e_su[0] = self.data_matrix_e_su[:, ind_x]
                    self.data_lines_yz_e_d[0] = self.data_matrix_e_d[:, ind_x]
                    self.data_lines_yz_e_sd[0] = self.data_matrix_e_sd[:, ind_x]
                    self.data_lines_yz_e_sum[0] = self.data_matrix_e_sum[:, ind_x]
                    self.data_lines_yz_e_ssum[0] = self.data_matrix_e_ssum[:, ind_x]
                    self.data_lines_yz_e_diff[0] = self.data_matrix_e_diff[:, ind_x]
                    self.data_lines_yz_e_sdiff[0] = self.data_matrix_e_sdiff[:, ind_x]
                    
        if cond_m:
            if cond_y:
                if (ind_y in range(self.data_matrix_m_u.shape[0])):
                    self.data_lines_xz_m_u[0] = self.data_matrix_m_u[ind_y, :]
                    self.data_lines_xz_m_d[0] = self.data_matrix_m_d[ind_y, :]
                    self.data_lines_xz_m_sum[0] = self.data_matrix_m_sum[ind_y, :]
                    self.data_lines_xz_m_diff[0] = self.data_matrix_m_diff[ind_y, :]
            if cond_x:
                if (ind_x in range(self.data_matrix_m_u.shape[1])):
                    self.data_lines_yz_m_u[0] = self.data_matrix_m_u[:, ind_x]
                    self.data_lines_yz_m_d[0] = self.data_matrix_m_d[:, ind_x]
                    self.data_lines_yz_m_sum[0] = self.data_matrix_m_sum[:, ind_x]
                    self.data_lines_yz_m_diff[0] = self.data_matrix_m_diff[:, ind_x]
        
    def get_lines_by_values(self, val_x, val_y):
        ind_x, ind_y = self.value_to_index(val_x, val_y)
        self.get_lines_by_index(ind_x, ind_y)
        
    def index_to_value(self, ind_x, ind_y):
        iint_x = int(round(ind_x))
        iint_y = int(round(ind_y))
        val_x = self.data_range_x[iint_x]
        val_y = self.data_range_y[iint_y]
        return val_x, val_y
        
    def value_to_index(self, val_x, val_y):
        min_x = numpy.nanmin(self.data_range_x)
        max_x = numpy.nanmax(self.data_range_x)
        min_y = numpy.nanmin(self.data_range_y)
        max_y = numpy.nanmax(self.data_range_y)
        ind_x = int(round((self.data_range_x.size-1)*(val_x-min_x)*1./(max_x-min_x)))
        ind_y = int(round((self.data_range_y.size-1)*(val_y-min_y)*1./(max_y-min_y)))
        return ind_x, ind_y
    
class cgraph(matplotlib.backends.backend_qt5agg.FigureCanvasQTAgg):
 
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = matplotlib.figure.Figure(figsize=(width, height), dpi = dpi)
        self.fig.subplots_adjust(left = 0.12,
                            right = 0.90,
                            top = 0.90,
                            bottom = 0.12,
                            wspace = 0.0,
                            hspace = 0.0)
        super(cgraph, self).__init__(self.fig)
        self.xlim_1 = None
        self.xlim_2 = None
        self.ylim_1 = None
        self.ylim_2 = None
        self.zlim_1 = None
        self.zlim_2 = None
        self.lim_flag = True
        
        self.val_x = None
        self.val_y = None
        self.val_z = None
        
        self.info_press = (False, False)
        self.control = parent
        self.figure = self.fig
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        self.ax_pri = self.fig.add_subplot(111)
        #fig.canvas.mpl_connect("button_press_event", self.onclick)
        
        
    def onclick(self, event):
        print "x position: ", event.xdata
        print "y position: ", event.ydata
        
    def set_x(self, np_x):
        self.val_x = np_x
        
    def set_y(self, np_y):
        self.val_y = np_y
        
    def set_z(self, np_z):
        self.val_z = np_z
    
    def set_lim(self, xlim_1, xlim_2, ax):
        
        if xlim_1 != None:
            if ax == "x":
                self.xlim_1 = xlim_1
            elif ax == "y":
                self.ylim_1 = xlim_1
            elif ax == "z":
                self.zlim_1 = xlim_1
        if xlim_2 != None:
            if ax == "x":
                self.xlim_2 = xlim_2
            elif ax == "y":
                self.ylim_2 = xlim_2
            elif ax == "z":
                self.zlim_2 = xlim_2
        if ax == "x":
            if (self.xlim_1 != None) & (self.xlim_2 != None):
                x_min = min([self.xlim_1, self.xlim_2])
                x_max = max([self.xlim_1, self.xlim_2])
                self.ax_pri.set_xlim([x_min, x_max])
                self.draw()        
        elif ax == "y":
            if (self.ylim_1 != None) & (self.ylim_2 != None):
                x_min = min([self.ylim_1, self.ylim_2])
                x_max = max([self.ylim_1, self.ylim_2])
                self.ax_pri.set_ylim([x_min, x_max])
                self.draw()        
        elif ax == "z":
            if (self.zlim_1 != None) & (self.zlim_2 != None):
                x_min = min([self.zlim_1, self.zlim_2])
                x_max = max([self.zlim_1, self.zlim_2])
                self.ax_pri.set_zlim([x_min, x_max])
                self.draw()        


class cwidg_ref_2d(QtWidgets.QWidget):
    def __init__(self, parent = None, np_x = None, np_y = None, np_z_2d = None, lim_x = None, lim_y = None):    
        super(cwidg_ref_2d, self).__init__(parent, QtCore.Qt.Window)
        self.init_layout_central(np_x, np_y, np_z_2d, lim_x, lim_y)
        self.setLayout(self.layout_central)
        self.mode = ["None"]# it can be '["bkg", 0], ["bkg", 1], ["bkg", 2], ["gauss", 0, "pos"], ["gauss", 0, "sig"]
        self.ctemp1 = 0.0 #temporary variable
        self.ctemp2 = 0.0 #temporary variable
        
    def init_layout_central(self, np_x, np_y, np_z_2d, lim_x, lim_y):
        lay_main = QtWidgets.QVBoxLayout()
        self.graph_m = cgraph_matrix(self)
        self.graph_m.figure.canvas.mpl_connect("button_press_event", self.onclick_matrix)
        
        self.graph_m_2 = cgraph_matrix(self)
        self.graph_m_2.hide()
        
        self.graph_l_1 = cgraph_line(self)
        self.graph_l_1.hide()
        
        self.graph_m.set_x(np_x)
        self.graph_m.set_y(np_y)
        self.graph_m.set_z(np_z_2d)
        self.graph_m.set_lim(lim_x[0], lim_x[1], "x")
        self.graph_m.set_lim(lim_y[0], lim_y[1], "y")
        self.graph_m.set_data_to_graph()
        
        
        splitter_h = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter_v = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        
        splitter_v.addWidget(self.graph_m)
        splitter_v.addWidget(self.graph_m_2)
        splitter_v.addWidget(self.graph_l_1)
        splitter_h.addWidget(splitter_v)

        lay_param = self.init_lay_param()
        widg_param = QtWidgets.QWidget()
        widg_param.setLayout(lay_param)
        splitter_h.addWidget(widg_param)
        
        
        lay_main.addWidget(splitter_h)
        
        
        self.layout_central = lay_main
        
    def init_lay_param(self):
        lay = QtWidgets.QVBoxLayout()
        
        q_pb = QtWidgets.QPushButton("set background by 3 points")
        q_pb.clicked.connect(lambda: self.set_mode(["bkg",0]))
        lay.addWidget(q_pb)
            
        lay_b_g =QtWidgets.QGridLayout()
        q_l = QtWidgets.QLabel("coord x")
        lay_b_g.addWidget(q_l, 0, 1)

        q_l = QtWidgets.QLabel("coord y")
        lay_b_g.addWidget(q_l, 0, 2)

        q_l = QtWidgets.QLabel("coord z")
        lay_b_g.addWidget(q_l, 0, 3)
        
        q_l = QtWidgets.QLabel("point 1")
        lay_b_g.addWidget(q_l, 1, 0)
        self.q_le_b1_x = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b1_x, 1, 1)
        self.q_le_b1_y = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b1_y, 1, 2)
        self.q_le_b1_z = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b1_z, 1, 3)

        q_l = QtWidgets.QLabel("point 2")
        lay_b_g.addWidget(q_l, 2, 0)
        self.q_le_b2_x = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b2_x, 2, 1)
        self.q_le_b2_y = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b2_y, 2, 2)
        self.q_le_b2_z = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b2_z, 2, 3)

        q_l = QtWidgets.QLabel("point 3")
        lay_b_g.addWidget(q_l, 3, 0)
        self.q_le_b3_x = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b3_x, 3, 1)
        self.q_le_b3_y = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b3_y, 3, 2)
        self.q_le_b3_z = QtWidgets.QLineEdit()
        lay_b_g.addWidget(self.q_le_b3_z, 3, 3)

        lay.addLayout(lay_b_g)
        q_pb = QtWidgets.QPushButton("refine background")
        q_pb.clicked.connect(self.fit_bkg)
        lay.addWidget(q_pb)
        
        
        
        q_l = QtWidgets.QLabel("number of reflections:")
        lay.addWidget(q_l)

        self.q_sb = QtWidgets.QSpinBox()
        self.q_sb.setMinimum(1)
        self.q_sb.setMaximum(4)
        self.q_sb.setValue(1)
        self.q_sb.valueChanged.connect(self.nr_changed)
        lay.addWidget(self.q_sb)

        q_pb = QtWidgets.QPushButton("define starting position and width for each intensity")
        q_pb.clicked.connect(lambda: self.set_mode(["gauss",0,"pos"]))
        lay.addWidget(q_pb)
        
        q_tab = QtWidgets.QTableWidget()
        q_tab.setRowCount(1)
        q_tab.setColumnCount(3)
        q_tab.setHorizontalHeaderLabels(('position', 'p_1', 'p_2'))
        self.q_tab_pos = q_tab
        lay.addWidget(q_tab)

        q_tab = QtWidgets.QTableWidget()
        q_tab.setRowCount(1)
        q_tab.setColumnCount(3)
        q_tab.setHorizontalHeaderLabels(('FWHM', 'p_1', 'p_2'))
        self.q_tab_sig = q_tab
        lay.addWidget(q_tab)

        q_tab = QtWidgets.QTableWidget()
        q_tab.setRowCount(1)
        q_tab.setColumnCount(3)
        q_tab.setHorizontalHeaderLabels(('int. intensity', 'p_1', 'p_2'))
        self.q_tab_iint = q_tab
        lay.addWidget(q_tab)
        
        q_pb = QtWidgets.QPushButton("fit model to experiment")
        q_pb.clicked.connect(self.fit_model)
        lay.addWidget(q_pb)
        lay.addStretch(1)
        q_pb = QtWidgets.QPushButton("plot model")
        q_pb.clicked.connect(self.plot_model)
        lay.addWidget(q_pb)
        q_pb = QtWidgets.QPushButton("plot exp.-model")
        q_pb.clicked.connect(self.plot_diff)
        lay.addWidget(q_pb)
        
        q_pb = QtWidgets.QPushButton("plot iint(phi)")
        q_pb.clicked.connect(self.plot_iint)
        lay.addWidget(q_pb)

        q_pb = QtWidgets.QPushButton("plot pos(phi)")
        q_pb.clicked.connect(self.plot_pos)
        lay.addWidget(q_pb)

        q_pb = QtWidgets.QPushButton("plot FWHM(phi)")
        q_pb.clicked.connect(self.plot_sig)
        lay.addWidget(q_pb)
        
        
        
        
        return lay

    def fit_bkg(self):
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        lpar_0 = [b1[2], b2[2], b3[2]]
            
        def temp_func(lpar):
            b1_t = (b1[0], b1[1], lpar[0])
            b2_t = (b2[0], b2[1], lpar[1])
            b3_t = (b3[0], b3[1], lpar[2])
            res_2d = func_gausses_2d_bkg(np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1_t, b2_t, b3_t)
            chi_sq = ((np_int_2d_exp - res_2d)/np_int_2d_sexp)**2
            return chi_sq.sum()
        
        res = scipy.optimize.minimize(temp_func, lpar_0, method='Nelder-Mead')
        lpar = [float(hh) for hh in res.x]
        
        b1_t = (b1[0], b1[1], lpar[0])
        b2_t = (b2[0], b2[1], lpar[1])
        b3_t = (b3[0], b3[1], lpar[2])
        
        self.load_param_model(lp_iint, lp_tth0, lp_sig, b1_t, b2_t, b3_t)
    
    def fit_model(self):
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        nr = len(lp_iint)
        lpar_0 = []
        for hh1, hh2, hh3 in zip(lp_iint, lp_tth0, lp_sig):
            lpar_0.append(hh1[0])
            lpar_0.append(hh1[1])
            lpar_0.append(hh1[2])
            lpar_0.append(hh2[0])
            lpar_0.append(hh2[1])
            lpar_0.append(hh2[2])
            lpar_0.append(hh3[0])
            lpar_0.append(hh3[1])
            lpar_0.append(hh3[2])
            
        def temp_func(lpar):
            lp_iint_t, lp_tth0_t, lp_sig_t = [], [], []
            for ir in range(nr):
                lp_iint_t.append((lpar[ir*9+0], lpar[ir*9+1], lpar[ir*9+2]))
                lp_tth0_t.append((lpar[ir*9+3], lpar[ir*9+4], lpar[ir*9+5]))
                lp_sig_t.append((lpar[ir*9+6], lpar[ir*9+7], lpar[ir*9+8]))
            res_2d = func_gausses_2d_bkg(np_tth, np_phi, lp_iint_t, lp_tth0_t, lp_sig_t, b1, b2, b3)
            chi_sq = ((np_int_2d_exp - res_2d)/np_int_2d_sexp)**2
            return chi_sq.sum()
        
        res = scipy.optimize.minimize(temp_func, lpar_0, method='Nelder-Mead')
        lpar = [float(hh) for hh in res.x]
        
        lp_iint_t, lp_tth0_t, lp_sig_t = [], [], []
        for ir in range(nr):
            lp_iint_t.append((lpar[ir*9+0], lpar[ir*9+1], lpar[ir*9+2]))
            lp_tth0_t.append((lpar[ir*9+3], lpar[ir*9+4], lpar[ir*9+5]))
            lp_sig_t.append((lpar[ir*9+6], lpar[ir*9+7], lpar[ir*9+8]))
        self.load_param_model(lp_iint_t, lp_tth0_t, lp_sig_t, b1, b2, b3)

    def load_param_model(self, lp_iint_t, lp_tth0_t, lp_sig_t, b1, b2, b3):
        nr = len(lp_iint_t)
        for ir, p_iint, p_tth0, p_sig in zip(range(nr), lp_iint_t, lp_tth0_t, lp_sig_t):
            self.q_tab_pos.setItem(ir, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(p_tth0[0])))
            self.q_tab_pos.setItem(ir, 1, QtWidgets.QTableWidgetItem("{:.3f}".format(p_tth0[1])))
            self.q_tab_pos.setItem(ir, 2, QtWidgets.QTableWidgetItem("{:.3f}".format(p_tth0[2])))
            self.q_tab_sig.setItem(ir, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(p_sig[0])))
            self.q_tab_sig.setItem(ir, 1, QtWidgets.QTableWidgetItem("{:.3f}".format(p_sig[1])))
            self.q_tab_sig.setItem(ir, 2, QtWidgets.QTableWidgetItem("{:.3f}".format(p_sig[2])))
            self.q_tab_iint.setItem(ir, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(p_iint[0])))
            self.q_tab_iint.setItem(ir, 1, QtWidgets.QTableWidgetItem("{:.3f}".format(p_iint[1])))
            self.q_tab_iint.setItem(ir, 2, QtWidgets.QTableWidgetItem("{:.3f}".format(p_iint[2])))
            
        self.q_le_b1_z.setText("{:.3f}".format(b1[2]))
        self.q_le_b2_z.setText("{:.3f}".format(b2[2]))
        self.q_le_b3_z.setText("{:.3f}".format(b3[2]))
        return
        
    def plot_model(self):
        self.graph_m_2.show()
        self.graph_l_1.hide()
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        res_2d = func_gausses_2d_bkg(np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3)
        self.graph_m_2.set_x(np_tth)
        self.graph_m_2.set_y(np_phi)
        self.graph_m_2.set_z(res_2d)
        self.graph_m_2.xlim_1 = None
        self.graph_m_2.xlim_2 = None
        self.graph_m_2.ylim_1 = None
        self.graph_m_2.ylim_2 = None
        self.graph_m_2.set_data_to_graph()

    def plot_iint(self):
        self.graph_m_2.hide()
        self.graph_l_1.show()
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        np_phi_rad = np_phi * math.pi/180.
        np_iint = lp_iint[0][0] + lp_iint[0][1] * numpy.cos(np_phi_rad)**2 + lp_iint[0][2]*numpy.cos(np_phi_rad)**4
        self.graph_l_1.set_x(np_phi)
        self.graph_l_1.set_y(np_iint)
        self.graph_l_1.set_data_to_graph()

    def plot_pos(self):
        self.graph_m_2.hide()
        self.graph_l_1.show()
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        np_phi_rad = np_phi * math.pi/180.
        np_tth0 = lp_tth0[0][0] + lp_tth0[0][1] * np_phi_rad + lp_tth0[0][2]*np_phi_rad**2
        self.graph_l_1.set_x(np_phi)
        self.graph_l_1.set_y(np_tth0)
        self.graph_l_1.set_data_to_graph()

    def plot_sig(self):
        self.graph_m_2.hide()
        self.graph_l_1.show()
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        np_phi_rad = np_phi * math.pi/180.
        np_sig = lp_sig[0][0] + lp_sig[0][1] * np_phi_rad + lp_sig[0][2]*np_phi_rad**2
        self.graph_l_1.set_x(np_phi)
        self.graph_l_1.set_y(np_sig)
        self.graph_l_1.set_data_to_graph()
        
    def plot_diff(self):
        self.graph_m_2.show()
        self.graph_l_1.hide()
        np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3 = self.take_param_model()
        res_2d = func_gausses_2d_bkg(np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3)
        self.graph_m_2.set_x(np_tth)
        self.graph_m_2.set_y(np_phi)
        self.graph_m_2.set_z(np_int_2d_exp-res_2d)
        self.graph_m_2.xlim_1 = None
        self.graph_m_2.xlim_2 = None
        self.graph_m_2.ylim_1 = None
        self.graph_m_2.ylim_2 = None
        self.graph_m_2.set_data_to_graph()        

    def take_param_model(self):
        xlim_min = min([self.graph_m.xlim_1, self.graph_m.xlim_2])
        xlim_max = max([self.graph_m.xlim_1, self.graph_m.xlim_2])
        ylim_min = min([self.graph_m.ylim_1, self.graph_m.ylim_2])
        ylim_max = max([self.graph_m.ylim_1, self.graph_m.ylim_2])
        np_tth = self.graph_m.val_x[xlim_min:xlim_max]
        np_phi = self.graph_m.val_y[ylim_min:ylim_max]
        np_int_2d_exp = self.graph_m.val_z[ylim_min:ylim_max, xlim_min:xlim_max]
        np_int_2d_sexp = numpy.ones(np_int_2d_exp.shape, float)
        
        b1 = (float(self.q_le_b1_x.text()), float(self.q_le_b1_y.text()), float(self.q_le_b1_z.text()))
        b2 = (float(self.q_le_b2_x.text()), float(self.q_le_b2_y.text()), float(self.q_le_b2_z.text())) 
        b3 = (float(self.q_le_b3_x.text()), float(self.q_le_b3_y.text()), float(self.q_le_b3_z.text()))
        
        nr = self.q_sb.value()
        lp_iint, lp_tth0, lp_sig = [], [], []
        for ir in range(nr):
            pos = (float(self.q_tab_pos.item(ir, 0).text()), float(self.q_tab_pos.item(ir, 1).text()), float(self.q_tab_pos.item(ir, 2).text()))
            sig = (float(self.q_tab_sig.item(ir, 0).text()), float(self.q_tab_sig.item(ir, 1).text()), float(self.q_tab_sig.item(ir, 2).text()))
            iint = (float(self.q_tab_iint.item(ir, 0).text()), float(self.q_tab_iint.item(ir, 1).text()), float(self.q_tab_iint.item(ir, 2).text()))
            lp_tth0.append(pos)
            lp_sig.append(sig)
            lp_iint.append(iint)
        

        return np_int_2d_exp, np_int_2d_sexp, np_tth, np_phi, lp_iint, lp_tth0, lp_sig, b1, b2, b3
            
        
    def nr_changed(self, nr):
        self.q_tab_pos.setRowCount(nr)
        self.q_tab_sig.setRowCount(nr)
        self.q_tab_iint.setRowCount(nr)
        return
    
    def onclick_matrix(self, event):
        ind_x = int(round(event.xdata))
        ind_y = int(round(event.ydata))
        val_x = self.graph_m.val_x[ind_x]
        val_y = self.graph_m.val_y[ind_y]
        val_z = self.graph_m.val_z[ind_y, ind_x]
        mode = self.mode
        if mode[0] == "bkg":
            if mode[1] == 0:
                self.q_le_b1_x.setText("{:.3f}".format(val_x))
                self.q_le_b1_y.setText("{:.3f}".format(val_y))
                self.q_le_b1_z.setText("{:.3f}".format(val_z))
                mode[1] = 1
            elif mode[1] == 1:
                self.q_le_b2_x.setText("{:.3f}".format(val_x))
                self.q_le_b2_y.setText("{:.3f}".format(val_y))
                self.q_le_b2_z.setText("{:.3f}".format(val_z))
                mode[1] = 2
            elif mode[1] == 2:
                self.q_le_b3_x.setText("{:.3f}".format(val_x))
                self.q_le_b3_y.setText("{:.3f}".format(val_y))
                self.q_le_b3_z.setText("{:.3f}".format(val_z))
                mode[0] = "None"
        elif mode[0] == "gauss":
            nr = mode[1]
            nr_max = self.q_sb.value()
            if nr > nr_max:
                mode[0] = "None"
                return
            else:
                if mode[2] == "pos":
                    self.q_tab_pos.setItem(nr, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(val_x)))
                    self.q_tab_pos.setItem(nr, 1, QtWidgets.QTableWidgetItem("0.000"))
                    self.q_tab_pos.setItem(nr, 2, QtWidgets.QTableWidgetItem("0.000"))
                    self.ctemp1 = val_x
                    self.ctemp2 = val_z
                    mode[2] = "sig"
                elif mode[2] == "sig":
                    sig = abs(val_x-self.ctemp1)
                    self.q_tab_sig.setItem(nr, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(sig)))
                    self.q_tab_sig.setItem(nr, 1, QtWidgets.QTableWidgetItem("0.000"))
                    self.q_tab_sig.setItem(nr, 2, QtWidgets.QTableWidgetItem("0.000"))
                    
                    iint = (self.ctemp2 - val_z)*(2.*math.pi)**0.5*sig*1./(1.-math.exp(-0.5))
                    
                    self.q_tab_iint.setItem(nr, 0, QtWidgets.QTableWidgetItem("{:.3f}".format(iint)))
                    self.q_tab_iint.setItem(nr, 1, QtWidgets.QTableWidgetItem("0.000"))
                    self.q_tab_iint.setItem(nr, 2, QtWidgets.QTableWidgetItem("0.000"))
                    mode[2] = "pos"
                    mode[1] = nr + 1
        
    def set_mode(self, mode):
        self.mode = mode
        
class cwidg_central(QtWidgets.QWidget):
    def __init__(self, fname_exp = None, fname_mod = None):
        super(cwidg_central, self).__init__()
        self.model = cmodel()
        self.init_layout_central()
        self.setLayout(self.layout_central)
        self.mode = ""#it can be 'e_u', 'e_d', 'e_sum', 'e_diff', 'm_u', 'm_d', 'm_sum', 'm_diff'
        self.fdir_xml = None
        
        if fname_mod != None:
            self.model.read_matrices(fname_mod, "mod")
            self.fdir_xml = os.path.dirname(fname_mod)
            self.plot_m_u()
            
        if fname_exp != None:
            self.model.read_matrices(fname_exp, "exp")
            self.fdir_xml = os.path.dirname(fname_exp)
            self.plot_e_u()
        
 
    def init_layout_central(self):
        
        lay_main = QtWidgets.QVBoxLayout()
        
        
        self.graph_m = cgraph_matrix(self)
        self.graph_l_tth = cgraph_line(self)
        self.graph_l_phi = cgraph_line(self)
        self.graph_m.figure.canvas.mpl_connect("button_press_event", self.onclick_matrix)
        self.graph_l_tth.figure.canvas.mpl_connect("button_press_event", self.onclick_l_tth)
        self.graph_l_phi.figure.canvas.mpl_connect("button_press_event", self.onclick_l_phi)
        """
        signal_pl = ccomunicate()
        self.graph_m.sig = signal_pl
        self.graph_l_tth.sig = signal_pl
        self.graph_l_phi.sig = signal_pl
        """

        splitter_h = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        splitter_v = QtWidgets.QSplitter(QtCore.Qt.Vertical)


        splitter_v.addWidget(self.graph_l_tth)
        splitter_v.addWidget(self.graph_l_phi)
        
        splitter_h.addWidget(self.graph_m)
        splitter_h.addWidget(splitter_v)

        lay_cp = self.control_panel()
        
        lay_1 = QtWidgets.QVBoxLayout()
        lay_1.addWidget(splitter_h)
        lay_1.addLayout(lay_cp)
        lay_main.addLayout(lay_1)

        self.layout_central = lay_main
        
    def onclick_matrix(self, event):
        ind_x = int(round(event.xdata))
        ind_y = int(round(event.ydata))
        val_x, val_y = self.model.index_to_value(ind_x, ind_y)
        mode = self.mode
        str_x = "{:.3f}".format(val_x)
        str_y = "{:.3f}".format(val_y)
    
        val_z = None
        if mode == "e_u":
            val_z = self.model.data_matrix_e_u[ind_y, ind_x]
        elif mode == "e_d":
            val_z = self.model.data_matrix_e_d[ind_y, ind_x]
        elif mode == "e_sum":
            val_z = self.model.data_matrix_e_sum[ind_y, ind_x]
        elif mode == "e_diff":
            val_z = self.model.data_matrix_e_diff[ind_y, ind_x]
        elif mode == "m_u":
            val_z = self.model.data_matrix_m_u[ind_y, ind_x]
        elif mode == "m_d":
            val_z = self.model.data_matrix_m_d[ind_y, ind_x]
        elif mode == "m_sum":
            val_z = self.model.data_matrix_m_sum[ind_y, ind_x]
        elif mode == "m_diff":
            val_z = self.model.data_matrix_m_diff[ind_y, ind_x]

        str_z = "{:.3f}".format(val_z)
        self.q_le_coord_x.setText(str_x)
        self.q_le_coord_y.setText(str_y)
        self.q_le_coord_z.setText(str_z)

        if event.button == 1:
            self.plot_lines(ind_x, ind_y, "xy")

        elif event.button == 2:
            self.graph_m.set_lim(0, len(self.graph_m.val_x)-1, "x")
            self.graph_m.set_lim(0, len(self.graph_m.val_y)-1, "y")
            self.graph_m.lim_flag = True
            self.graph_l_tth.set_lim(min(self.graph_l_tth.val_x), max(self.graph_l_tth.val_x), "x")
            self.graph_l_tth.lim_flag = True
            self.graph_l_phi.set_lim(min(self.graph_l_phi.val_x), max(self.graph_l_phi.val_x), "x")
            self.graph_l_phi.lim_flag = True
            
            
        elif event.button == 3:
            if self.graph_m.lim_flag:
                self.graph_m.lim_flag = False
                self.graph_m.xlim_1 = ind_x
                self.graph_m.ylim_1 = ind_y

                self.graph_l_tth.lim_flag = False
                self.graph_l_tth.xlim_1 = val_x

                self.graph_l_phi.lim_flag = False
                self.graph_l_phi.xlim_1 = val_y
            else:
                self.graph_m.lim_flag = True
                self.graph_m.set_lim(None, ind_x, "x")
                self.graph_m.set_lim(None, ind_y, "y")
                
                self.graph_l_tth.lim_flag = True
                self.graph_l_tth.set_lim(None, val_x, "x")

                self.graph_l_phi.lim_flag = True
                self.graph_l_phi.set_lim(None, val_y, "x")
        
    def onclick_l_tth(self, event):
        dblclick = event.dblclick
        val_x = event.xdata
        val_z = event.ydata
        str_x = "{:.3f}".format(val_x)
        str_z = "{:.3f}".format(val_z)
        self.q_le_coord_x.setText(str_x)
        self.q_le_coord_z.setText(str_z)
        try:
            val_y = float(self.q_le_coord_y.text())
        except:
            return    

        ind_x, ind_y = self.model.value_to_index(val_x, val_y)
        
        if ((event.button == 1) & (not(dblclick))):
            self.plot_lines(ind_x, ind_y, "y")
            
        elif event.button == 2:
            self.graph_l_tth.lim_flag = True
            self.graph_l_tth.set_lim(min(self.graph_l_tth.val_x), max(self.graph_l_tth.val_x), "x")
            self.graph_m.lim_flag = True
            self.graph_m.set_lim(0, len(self.graph_m.val_x)-1, "x")
            
        elif event.button == 3:
            if self.graph_l_tth.lim_flag:
                self.graph_l_tth.lim_flag = False
                self.graph_l_tth.xlim_1 = val_x

                self.graph_m.lim_flag = False
                self.graph_m.xlim_1 = ind_x
            else:
                self.graph_l_tth.lim_flag = True
                self.graph_l_tth.set_lim(None, val_x, "x")

                self.graph_m.lim_flag = True
                self.graph_m.set_lim(None, ind_x, "x")
                
        elif ((event.button == 1) & dblclick):
            if self.graph_l_tth.gaus_flag == 0:
                self.graph_l_tth.gaus_flag = 1
                self.graph_l_tth.gaus_i_max_x = val_x
                self.graph_l_tth.gaus_i_max_y = val_z
                self.graph_l_tth.add_scatter_to_graph(val_x, val_z)
            elif self.graph_l_tth.gaus_flag == 1:
                self.graph_l_tth.gaus_flag = 2
                self.graph_l_tth.gaus_bord_l_x = val_x
                self.graph_l_tth.gaus_bord_l_y = val_z
                self.graph_l_tth.add_scatter_to_graph(val_x, val_z)
            elif self.graph_l_tth.gaus_flag == 2:
                self.graph_l_tth.gaus_flag = 0
                self.graph_l_tth.gaus_bord_r_x = val_x
                self.graph_l_tth.gaus_bord_r_y = val_z
                self.graph_l_tth.add_scatter_to_graph(val_x, val_z)
                np_x, np_y, bl, br, pos, iint, width = refine_by_gauss(
                    self.graph_l_tth.val_x, self.graph_l_tth.val_y, self.graph_l_tth.val_sy,
                    self.graph_l_tth.gaus_i_max_x, self.graph_l_tth.gaus_i_max_y,
                    self.graph_l_tth.gaus_bord_l_x, self.graph_l_tth.gaus_bord_l_y,
                    self.graph_l_tth.gaus_bord_r_x,self.graph_l_tth.gaus_bord_r_y)
                self.q_le_bl.setText("{:.3f}".format(bl))
                self.q_le_br.setText("{:.3f}".format(br))
                self.q_le_pos.setText("{:.3f}".format(pos))
                self.q_le_iint.setText("{:.3f}".format(iint))
                self.q_le_width.setText("{:.3f}".format(width))
                self.graph_l_tth.ax_pri.plot(np_x, np_y, "-")
                self.graph_l_tth.draw()
            
    def onclick_l_phi(self, event):
        dblclick = event.dblclick
        val_y = event.xdata
        val_z = event.ydata
        str_y = "{:.3f}".format(val_y)
        str_z = "{:.3f}".format(val_z)
        self.q_le_coord_y.setText(str_y)
        self.q_le_coord_z.setText(str_z)
        try:
            val_x = float(self.q_le_coord_x.text())
        except:
            return

        ind_x, ind_y = self.model.value_to_index(val_x, val_y)
        if ((event.button == 1) & (not(dblclick))):
            self.plot_lines(ind_x, ind_y, "x")
        elif event.button == 2:
            self.graph_l_phi.lim_flag = True
            self.graph_l_phi.set_lim(min(self.graph_l_phi.val_x), max(self.graph_l_phi.val_x), "x")
            self.graph_m.lim_flag = True
            self.graph_m.set_lim(0, len(self.graph_m.val_y)-1, "y")
            
        elif event.button == 3:
            if self.graph_l_phi.lim_flag:
                self.graph_l_phi.lim_flag = False
                self.graph_l_phi.xlim_1 = val_y

                self.graph_m.lim_flag = False
                self.graph_m.ylim_1 = ind_y
            else:
                self.graph_l_phi.lim_flag = True
                self.graph_l_phi.set_lim(None, val_y, "x")

                self.graph_m.lim_flag = True
                self.graph_m.set_lim(None, ind_y, "y")

        elif ((event.button == 1) & dblclick):
            if self.graph_l_phi.gaus_flag == 0:
                self.graph_l_phi.gaus_flag = 1
                self.graph_l_phi.gaus_i_max_x = val_y
                self.graph_l_phi.gaus_i_max_y = val_z
                self.graph_l_phi.add_scatter_to_graph(val_y, val_z)
            elif self.graph_l_phi.gaus_flag == 1:
                self.graph_l_phi.gaus_flag = 2
                self.graph_l_phi.gaus_bord_l_x = val_y
                self.graph_l_phi.gaus_bord_l_y = val_z
                self.graph_l_phi.add_scatter_to_graph(val_y, val_z)
            elif self.graph_l_phi.gaus_flag == 2:
                self.graph_l_phi.gaus_flag = 0
                self.graph_l_phi.gaus_bord_r_x = val_y
                self.graph_l_phi.gaus_bord_r_y = val_z
                self.graph_l_phi.add_scatter_to_graph(val_y, val_z)
                np_x, np_y, bl, br, pos, iint, width = refine_by_gauss(
                    self.graph_l_phi.val_x, self.graph_l_phi.val_y, self.graph_l_phi.val_sy,
                    self.graph_l_phi.gaus_i_max_x, self.graph_l_phi.gaus_i_max_y,
                    self.graph_l_phi.gaus_bord_l_x, self.graph_l_phi.gaus_bord_l_y,
                    self.graph_l_phi.gaus_bord_r_x,self.graph_l_phi.gaus_bord_r_y)
                self.q_le_bl.setText("{:.3f}".format(bl))
                self.q_le_br.setText("{:.3f}".format(br))
                self.q_le_pos.setText("{:.3f}".format(pos))
                self.q_le_iint.setText("{:.3f}".format(iint))
                self.q_le_width.setText("{:.3f}".format(width))
                self.graph_l_phi.ax_pri.plot(np_x, np_y, "-")
                self.graph_l_phi.draw()


    def b_plot_sections(self):
        """
        button to plot sections
        """
        val_x = float(self.q_le_coord_x.text())
        val_y = float(self.q_le_coord_y.text())
        ind_x, ind_y = self.model.value_to_index(val_x, val_y)
        self.plot_lines(ind_x, ind_y, "xy")
        
    def plot_lines(self, ind_x, ind_y, plot_opt):
        """
        plot_opt is "x", "y" or "xy"
        """
        cond_x = "x" in plot_opt
        cond_y = "y" in plot_opt
        self.model.get_lines_by_index(ind_x, ind_y)
        mode = self.mode
        
        if cond_x:
            self.graph_l_tth.set_x(self.model.data_range_x)
        if cond_y:
            self.graph_l_phi.set_x(self.model.data_range_y)

        if mode[1:] == "_u":
            if isinstance(self.model.data_matrix_e_u, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y(self.model.data_lines_xz_e_u[-1])
                    self.graph_l_tth.set_sy(self.model.data_lines_xz_e_su[-1])
                if cond_y:
                    self.graph_l_phi.set_y(self.model.data_lines_yz_e_u[-1])
                    self.graph_l_phi.set_sy(self.model.data_lines_yz_e_su[-1]) 
                
            if isinstance(self.model.data_matrix_m_u, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y_m(self.model.data_lines_xz_m_u[-1])
                if cond_y:
                    self.graph_l_phi.set_y_m(self.model.data_lines_yz_m_u[-1])
               
        elif mode[1:] == "_d":
            if isinstance(self.model.data_matrix_e_d, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y(self.model.data_lines_xz_e_d[-1])
                    self.graph_l_tth.set_sy(self.model.data_lines_xz_e_sd[-1])
                if cond_y:
                    self.graph_l_phi.set_y(self.model.data_lines_yz_e_d[-1])
                    self.graph_l_phi.set_sy(self.model.data_lines_yz_e_sd[-1])
            
            if isinstance(self.model.data_matrix_m_d, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y_m(self.model.data_lines_xz_m_d[-1])
                if cond_y:
                    self.graph_l_phi.set_y_m(self.model.data_lines_yz_m_d[-1])
               
        elif mode[1:] == "_sum":
            if isinstance(self.model.data_matrix_e_sum, numpy.ndarray):
                if cond_x:  
                    self.graph_l_tth.set_y(self.model.data_lines_xz_e_sum[-1])
                    self.graph_l_tth.set_sy(self.model.data_lines_xz_e_ssum[-1])
                if cond_y:
                    self.graph_l_phi.set_y(self.model.data_lines_yz_e_sum[-1])
                    self.graph_l_phi.set_sy(self.model.data_lines_yz_e_ssum[-1])

            if isinstance(self.model.data_matrix_m_sum, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y_m(self.model.data_lines_xz_m_sum[-1])
                if cond_y:
                    self.graph_l_phi.set_y_m(self.model.data_lines_yz_m_sum[-1])
               
        elif mode[1:] == "_diff":
            if isinstance(self.model.data_matrix_e_diff, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y(self.model.data_lines_xz_e_diff[-1])
                    self.graph_l_tth.set_sy(self.model.data_lines_xz_e_sdiff[-1])
                if cond_y:
                    self.graph_l_phi.set_y(self.model.data_lines_yz_e_diff[-1])
                    self.graph_l_phi.set_sy(self.model.data_lines_yz_e_sdiff[-1])
                
            if isinstance(self.model.data_matrix_m_diff, numpy.ndarray):
                if cond_x:
                    self.graph_l_tth.set_y_m(self.model.data_lines_xz_m_diff[-1])
                if cond_y:
                    self.graph_l_phi.set_y_m(self.model.data_lines_yz_m_diff[-1])
                
        if cond_x:
            self.graph_l_tth.set_data_to_graph()
        if cond_y:
            self.graph_l_phi.set_data_to_graph()
        
    def control_panel(self):
        lay_2 = QtWidgets.QGridLayout()
        
        kk = 0
        q_pb = QtWidgets.QPushButton("open exp. data")
        q_pb.clicked.connect(self.open_exp_data)
        lay_2.addWidget(q_pb, 0, kk)
        q_pb = QtWidgets.QPushButton("open mod. data")
        q_pb.clicked.connect(self.open_mod_data)
        lay_2.addWidget(q_pb, 1, kk)

        q_pb = QtWidgets.QPushButton("clean data")
        q_pb.clicked.connect(self.clean_data)
        lay_2.addWidget(q_pb, 3, kk)

        kk += 1
        q_pb = QtWidgets.QPushButton("plot 'up' exp")
        q_pb.clicked.connect(self.plot_e_u)
        lay_2.addWidget(q_pb, 0, kk)
        q_pb = QtWidgets.QPushButton("plot 'down' exp")
        q_pb.clicked.connect(self.plot_e_d)
        lay_2.addWidget(q_pb, 1, kk)
        q_pb = QtWidgets.QPushButton("plot 'sum' exp")
        q_pb.clicked.connect(self.plot_e_sum)
        lay_2.addWidget(q_pb, 2, kk)
        q_pb = QtWidgets.QPushButton("plot 'diff' exp")
        q_pb.clicked.connect(self.plot_e_diff)
        lay_2.addWidget(q_pb, 3, kk)


        kk += 1
        q_pb = QtWidgets.QPushButton("plot 'up' mod")
        q_pb.clicked.connect(self.plot_m_u)
        lay_2.addWidget(q_pb, 0, kk)
        q_pb = QtWidgets.QPushButton("plot 'down' mod")
        q_pb.clicked.connect(self.plot_m_d)
        lay_2.addWidget(q_pb, 1, kk)
        q_pb = QtWidgets.QPushButton("plot 'sum' mod")
        q_pb.clicked.connect(self.plot_m_sum)
        lay_2.addWidget(q_pb, 2, kk)
        q_pb = QtWidgets.QPushButton("plot 'diff' mod")
        q_pb.clicked.connect(self.plot_m_diff)
        lay_2.addWidget(q_pb, 3, kk)

        
        kk += 1
        qlabel = QtWidgets.QLabel("coord x")
        lay_2.addWidget(qlabel, 0, kk)
        qlabel = QtWidgets.QLabel("coord y")
        lay_2.addWidget(qlabel, 1, kk)
        qlabel = QtWidgets.QLabel("coord z")
        lay_2.addWidget(qlabel, 2, kk)
        q_pb = QtWidgets.QPushButton("plot sections")
        q_pb.clicked.connect(self.b_plot_sections)
        lay_2.addWidget(q_pb, 3, kk)

        
        kk += 1
        self.q_le_coord_x = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_coord_x, 0, kk)
        self.q_le_coord_y = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_coord_y, 1, kk)
        self.q_le_coord_z = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_coord_z, 2, kk)
        

        kk += 1
        q_pb = QtWidgets.QPushButton("set min. and max. values")
        q_pb.clicked.connect(self.set_limit)
        lay_2.addWidget(q_pb, 0, kk)
        self.q_le_val_min = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_val_min, 1, kk)
        self.q_le_val_max = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_val_max, 2, kk)


        kk += 1
        self.q_rb_x = QtWidgets.QRadioButton("for x")
        lay_2.addWidget(self.q_rb_x, 0, kk)
        self.q_rb_y = QtWidgets.QRadioButton("for y")
        lay_2.addWidget(self.q_rb_y, 1, kk)
        self.q_rb_z = QtWidgets.QRadioButton("for z")
        self.q_rb_z.setChecked(True)
        lay_2.addWidget(self.q_rb_z, 2, kk)


        kk += 1
        qlabel = QtWidgets.QLabel("b1 and b2")
        lay_2.addWidget(qlabel, 0, kk)
        self.q_le_bl = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_bl, 1, kk)
        self.q_le_br = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_br, 2, kk)


        kk += 1
        qlabel = QtWidgets.QLabel("position")
        lay_2.addWidget(qlabel, 0, kk)
        qlabel = QtWidgets.QLabel("width")
        lay_2.addWidget(qlabel, 1, kk)
        qlabel = QtWidgets.QLabel("intensity")
        lay_2.addWidget(qlabel, 2, kk)
        q_pb = QtWidgets.QPushButton("save to file")
        q_pb.clicked.connect(self.save_to_file)
        lay_2.addWidget(q_pb, 3, kk)


        kk += 1
        self.q_le_pos = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_pos, 0, kk)
        self.q_le_width = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_width, 1, kk)
        self.q_le_iint = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_iint, 2, kk)
        self.q_le_fout = QtWidgets.QLineEdit()
        lay_2.addWidget(self.q_le_fout, 3, kk)


        kk += 1
        q_pb = QtWidgets.QPushButton("x fit by Gauss")
        q_pb.clicked.connect(self.gauss_fit_x)
        lay_2.addWidget(q_pb, 0, kk)
        q_pb = QtWidgets.QPushButton("y fit by Gauss")
        q_pb.clicked.connect(self.gauss_fit_y)
        lay_2.addWidget(q_pb, 1, kk)
        
        
        q_pb = QtWidgets.QPushButton("2d fit by Gauss")
        q_pb.clicked.connect(self.gauss_fit_2d)
        lay_2.addWidget(q_pb, 3, kk)
        
        return lay_2
            
        """
        np_tth = model.data_range_x
        np_phi = model.data_range_y
        np_int = model.data_matrix_m_u
        
            
        
        self.graph_m.set_x(np_tth)
        self.graph_m.set_y(np_phi)
        self.graph_m.set_z(np_int)
            
            
        
        self.graph_m.set_data_to_graph()

        self.graph_l_tth.set_x(np_tth)
        self.graph_l_tth.set_y(np_int[10,:])
        self.graph_l_tth.set_data_to_graph()
            
        self.graph_l_phi.set_x(np_phi)
        self.graph_l_phi.set_y(np_int[:,10])
        self.graph_l_phi.set_data_to_graph()

        
        signal_pl.sig_m.connect(self.rewdraw_sections)
        signal_pl.sig_l_redraw.connect(self.lim_change_for_l)
        """

    def clean_data(self):
        self.model = cmodel()
        
    def gauss_fit_x(self):
        bl_x = self.graph_l_tth.xlim_1
        br_x = self.graph_l_tth.xlim_2
        
        bl_0_y = float(self.q_le_bl.text())
        br_0_y = float(self.q_le_br.text())
        
        iint_0 = float(self.q_le_iint.text())
        pos_0 = float(self.q_le_pos.text())
        width_0 = float(self.q_le_width.text())
        
        int_0_x = pos_0
        int_0_y = iint_0/(width_0*(0.5*math.pi)**2) + 0.5*(bl_0_y+br_0_y)
        
        inp_x = self.graph_l_tth.val_x
        inp_y = self.graph_l_tth.val_y
        inp_sy = self.graph_l_tth.val_sy
        
        np_x, np_y, bl, br, pos, iint, width  = refine_by_gauss(
                inp_x, inp_y, inp_sy, int_0_x, int_0_y, 
                bl_x, bl_0_y, br_x, br_0_y)
        
        self.q_le_bl.setText("{:.3f}".format(bl))
        self.q_le_br.setText("{:.3f}".format(br))
        self.q_le_pos.setText("{:.3f}".format(pos))
        self.q_le_iint.setText("{:.3f}".format(iint))
        self.q_le_width.setText("{:.3f}".format(width))
        self.graph_l_tth.ax_pri.plot(np_x, np_y, "-")
        self.graph_l_tth.draw()        

    def gauss_fit_y(self):
        bl_x = self.graph_l_phi.xlim_1
        br_x = self.graph_l_phi.xlim_2
        
        bl_0_y = float(self.q_le_bl.text())
        br_0_y = float(self.q_le_br.text())
        
        iint_0 = float(self.q_le_iint.text())
        pos_0 = float(self.q_le_pos.text())
        width_0 = float(self.q_le_width.text())
        
        int_0_x = pos_0
        int_0_y = iint_0/(width_0*(0.5*math.pi)**2) + 0.5*(bl_0_y+br_0_y)
        
        inp_x = self.graph_l_phi.val_x
        inp_y = self.graph_l_phi.val_y
        inp_sy = self.graph_l_phi.val_sy
        
        np_x, np_y, bl, br, pos, iint, width  = refine_by_gauss(
                inp_x, inp_y, inp_sy, int_0_x, int_0_y, 
                bl_x, bl_0_y, br_x, br_0_y)
        
        self.q_le_bl.setText("{:.3f}".format(bl))
        self.q_le_br.setText("{:.3f}".format(br))
        self.q_le_pos.setText("{:.3f}".format(pos))
        self.q_le_iint.setText("{:.3f}".format(iint))
        self.q_le_width.setText("{:.3f}".format(width))
        self.graph_l_phi.ax_pri.plot(np_x, np_y, "-")
        self.graph_l_phi.draw()        

    def gauss_fit_2d(self):
        widg = cwidg_ref_2d(self, self.graph_m.val_x, self.graph_m.val_y, self.graph_m.val_z,
                            (self.graph_m.xlim_1, self.graph_m.xlim_2),
                            (self.graph_m.ylim_1, self.graph_m.ylim_2))
        widg.show()
        
        pass
        
    def save_to_file(self):
        fdir_xml = self.fdir_xml
        if fdir_xml == None:
            fdir_xml = os.getcwd()
        fname = os.path.join(fdir_xml, self.q_le_fout.text())
        cond = os.path.isfile(fname)
        fid = open(fname, "a")
        if (not cond):
            line = "     coordx      coordy      b_left     b_right   integ.int        pos.       width\n"
            fid.write(line)
        coord_x = float(self.q_le_coord_x.text())
        coord_y = float(self.q_le_coord_y.text())
        bl = float(self.q_le_bl.text())
        br = float(self.q_le_br.text())
        iint = float(self.q_le_iint.text())
        pos = float(self.q_le_pos.text())
        width = float(self.q_le_width.text())
        line = " {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f} \n".format(
                coord_x, coord_y, bl, br, iint, pos, width)
        
        fid.write(line)
        fid.close()
        
    def set_limit(self):
        val_min = float(self.q_le_val_min.text())
        val_max = float(self.q_le_val_max.text())
        
        if self.q_rb_x.isChecked():
            self.set_limit_x(val_min, val_max)
        elif self.q_rb_y.isChecked():
            self.set_limit_y(val_min, val_max)
        elif self.q_rb_z.isChecked():
            self.set_limit_z(val_min, val_max)
            
    def set_limit_x(self, val_min, val_max):
        ind_min, hh = self.model.value_to_index(val_min, 0.)
        ind_max, hh = self.model.value_to_index(val_max, 0.)
        self.graph_m.set_lim(ind_min, ind_max, "x")
        self.graph_l_tth.set_lim(val_min, val_max, "x")
        
    def set_limit_y(self, val_min, val_max):
        hh, ind_min = self.model.value_to_index(0., val_min)
        hh, ind_max = self.model.value_to_index(0., val_max)
        self.graph_m.set_lim(ind_min, ind_max, "y")
        self.graph_l_phi.set_lim(val_min, val_max, "x")

    def set_limit_z(self, val_min, val_max):
        self.graph_m.zlim_1 = val_min
        self.graph_m.zlim_2 = val_max
        self.graph_m.set_data_to_graph()
        self.graph_l_tth.set_lim(val_min, val_max, "y")
        self.graph_l_phi.set_lim(val_min, val_max, "y")

    def open_exp_data(self):

        
        fdirxml = self.fdir_xml
        if fdirxml == None:
            fdirxml = os.getcwd()
        fname, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdirxml,
                                                      "Extension file (*.*)")
        if (not ok):
            return        
        self.fdir_xml = os.path.dirname("{}".format(fname))
        self.model.read_matrices(fname, "exp")
        
        
    def open_mod_data(self):
        fname = "e_HoTi_bc.out"
        fname = "test_gn.dat"
        fname = "c_mat_ed_bc.out"

        fdirxml = self.fdir_xml
        if fdirxml == None:
            fdirxml = os.getcwd()
        fname, ok = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', fdirxml,
                                                      "Extension file (*.*)")
        if (not ok):
            return        
        self.fdir_xml = os.path.dirname("{}".format(fname))
        self.model.read_matrices(fname, "mod")
    
    def plot_e_u(self):
        self.mode = "e_u"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_e_u)
        self.graph_m.set_data_to_graph()

    def plot_e_d(self):
        self.mode = "e_d"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_e_d)
        self.graph_m.set_data_to_graph()        
        
    def plot_e_sum(self):
        self.mode = "e_sum"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_e_sum)
        self.graph_m.set_data_to_graph() 

    def plot_e_diff(self):
        self.mode = "e_diff"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_e_diff)
        self.graph_m.set_data_to_graph()        
        
    def plot_m_u(self):
        self.mode = "m_u"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_m_u)
        self.graph_m.set_data_to_graph()

    def plot_m_d(self):
        self.mode = "m_d"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_m_d)
        self.graph_m.set_data_to_graph()        

    def plot_m_sum(self):
        self.mode = "m_sum"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_m_sum)
        self.graph_m.set_data_to_graph()            
        
    def plot_m_diff(self):
        self.mode = "m_diff"
        self.graph_m.set_x(self.model.data_range_x)
        self.graph_m.set_y(self.model.data_range_y)
        self.graph_m.set_z(self.model.data_matrix_m_diff)
        self.graph_m.set_data_to_graph()        
        
    def change_value(self, value):
        print value
        return


        

class cgraph_line(cgraph):
    def __init__(self, parent=None):
        super(cgraph_line, self).__init__(cgraph)
        self.ylim_1 = None
        self.ylim_2 = None
        self.val_sy = None
        self.val_y_m = None
        self.gaus_flag = 0 #three type of flag 0, 1, 2 (center, left border, right border)
        self.gaus_i_max_x = None
        self.gaus_i_max_y = None
        self.gaus_bord_l_x = None
        self.gaus_bord_l_y = None
        self.gaus_bord_r_x = None
        self.gaus_bord_r_y = None

    def set_sy(self, np_y):
        self.val_sy = np_y

    def set_y_m(self, np_y):
        self.val_y_m = np_y

    """    
    def onclick(self, event):
        if event.button == 1:
            print "x point: ", event.xdata
            print "y point: ", event.ydata
            #ind_phi = int(event.ydata)
        else:
            pass            
    """    
    def set_data_to_graph(self):
        self.ax_pri.cla()
        #self.ax_pri.set_xlim([0., 90.])
        #self.ax_pri.plot(self.val_x, self.val_y, "-")
        if ((self.xlim_1 != None)&(self.xlim_2 != None)):
            x_min = min([self.xlim_1 , self.xlim_2])
            x_max = max([self.xlim_1 , self.xlim_2])
            self.ax_pri.set_xlim([x_min, x_max])
        if ((self.ylim_1 != None)&(self.ylim_2 != None)):
            y_min = min([self.ylim_1 , self.ylim_2])
            y_max = max([self.ylim_1 , self.ylim_2])
            self.ax_pri.set_ylim([y_min, y_max])        
            

        if isinstance(self.val_y_m, numpy.ndarray):
            self.ax_pri.plot(self.val_x, self.val_y_m, "-")
            
        if isinstance(self.val_sy, numpy.ndarray):
            col_1 = "#000000"
            self.ax_pri.errorbar(self.val_x, self.val_y, yerr = self.val_sy, ecolor = col_1, fmt='o', color=col_1, linewidth = 0.5)
            
        elif isinstance(self.val_y, numpy.ndarray):
            self.ax_pri.plot(self.val_x, self.val_y, ".")

        #self.ax_pri.set_xticks(self.np_tth)
        #self.ax_pri.set_yticks(self.np_phi)
        #print self.np_phi.size
        #print self.np_tth.size
        #print self.np_int.shape
        #self.ax_pri.plot(self.data_x, self.data_x, "k-", linewidth=1.0)    
        #self.ax_pri.errorbar(self.data_x, self.data_y, yerr = self.data_sy, ecolor = col_1, fmt='o', color=col_1, linewidth = 0.5)
        self.draw()
    def add_scatter_to_graph(self, val_x, val_y):
        self.ax_pri.scatter(val_x, val_y)
        self.draw()

class cgraph_matrix(cgraph):
    def __init__(self, parent=None):
        super(cgraph_matrix, self).__init__(cgraph)
        self.zlim_1 = None
        self.zlim_2 = None
        self.fig.subplots_adjust(left = 0.02,
                            right = 0.95,
                            top = 0.95,
                            bottom = 0.02,
                            wspace = 0.0,
                            hspace = 0.0)
        
        
        
        
    def onclick(self, event):
        ind_tth = int(event.xdata)
        ind_phi = int(event.ydata)
        if event.button == 1:
            self.val_x_inline = self.val_z[ind_phi, :]
            self.val_y_inline = self.val_z[:, ind_tth]
            print "ttheta: {:.2f}, phi: {:.2f} ".format(self.val_x[ind_tth], 
                           self.val_y[ind_phi])
            self.sig.sig_m.emit()
        elif event.button == 2:
            self.set_lim(0, len(self.val_x)-1, "x")
            self.set_lim(0, len(self.val_y)-1, "y")
            self.lim_flag = True
            self.sig.sig_l_redraw.emit()
            
        elif event.button == 3:
            ind_tth = int(event.xdata)
            ind_phi = int(event.ydata)
            if self.lim_flag:
                print "Step 1"
                self.lim_flag = False
                self.xlim_1 = ind_tth
                self.ylim_1 = ind_phi
            else:
                print "Step 2",[self.xlim_1, ind_tth], [self.ylim_1, ind_phi]
                self.lim_flag = True
                self.set_lim(None, ind_tth, "x")
                self.set_lim(None, ind_phi, "y")
                self.sig.sig_l_redraw.emit()

        
    def set_data_to_graph(self):
        self.ax_pri.cla()
        if (self.zlim_1 != None) and (self.zlim_2 != None):
            self.ax_pri.imshow(self.val_z, aspect = 'auto', origin = 'lower', vmin = self.zlim_1, vmax = self.zlim_2)
        else:
            self.ax_pri.imshow(self.val_z, aspect = 'auto', origin = 'lower')
        #self.ax_pri.set_axis_off()
        #numpy.arange(len(self.val_x)),self.val_x
        self.ax_pri.set_xticks([])
        
        self.ax_pri.set_yticks([])
        if ((self.xlim_1 != None)&(self.xlim_2 != None)):
            x_min = min([self.xlim_1 , self.xlim_2])
            x_max = max([self.xlim_1 , self.xlim_2])
            self.ax_pri.set_xlim([x_min, x_max])
        if ((self.ylim_1 != None)&(self.ylim_2 != None)):
            y_min = min([self.ylim_1 , self.ylim_2])
            y_max = max([self.ylim_1 , self.ylim_2])
            self.ax_pri.set_ylim([y_min, y_max])
        #self.ax_pri.plot(self.data_x, self.data_y, "-")
        #self.ax_pri.matshow(self.np_int)
        #self.ax_pri.set_xticks(self.np_tth)
        #self.ax_pri.set_yticks(self.np_phi)
        #print self.np_phi.size
        #print self.np_tth.size
        #print self.np_int.shape
        #self.ax_pri.plot(self.data_x, self.data_x, "k-", linewidth=1.0)    
        #self.ax_pri.errorbar(self.data_x, self.data_y, yerr = self.data_sy, ecolor = col_1, fmt='o', color=col_1, linewidth = 0.5)
        self.draw()

        


        
class cwind_central(QtWidgets.QMainWindow):
    def __init__(self, larg):
        super(cwind_central, self).__init__()
        self.title = "program 'matrix Graph'"
        self.setWindowTitle(self.title)
        fdir = os.getcwd()
        """
        fname_e = os.path.join(fdir, "e_Fe3O4.out")
        fname_m = os.path.join(fdir, "m_Fe3O4.out")
        widg_central = cwidg_central(fname_e, fname_m)
        """
        
        fname_e, fname_m = None, None
        if len(larg) == 2:
            fname_e = larg[1]
        elif len(larg) >= 3:
            fname_e, fname_m = larg[1], larg[2]
        
        widg_central = cwidg_central(fname_e, fname_m)
        
        self.setCentralWidget(widg_central)
        self.show()

 
if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    larg = sys.argv
    ex = cwind_central(larg)
    
    
    
    sys.exit(app.exec_())
    
    
    
"""
    def clean_fig(self):
        self.ldata = []
        self.ldata_simple = []
        self.pcanvas.figure.clf()
        self.pcanvas.make_plot()
        self.pcanvas.draw()
        
    def replot_fig(self):
        lhelp=['#005555','#BB5555','#550055','#55BB55','#555500','#5555BB']
        for iddata, ddata in enumerate(self.ldata_simple):
            lflags = [ (isinstance(hh1,(float, int))&isinstance(hh2,(float, int))&isinstance(hh3,(float, int))&isinstance(hh4,(float, int))) for hh1,hh2,hh3,hh4 in zip(ddata["IntUP"],ddata["sIntUP"],ddata["IntDOWN"],ddata["sIntDOWN"])]
            lx = [hh for hh, flag in zip(ddata["ttheta"], lflags) if flag]
            li_u_exp = [hh for hh, flag in zip(ddata["IntUP"], lflags) if flag]
            lsi_u_exp = [hh for hh, flag in zip(ddata["sIntUP"], lflags) if flag]
            li_d_exp = [hh for hh, flag in zip(ddata["IntDOWN"], lflags) if flag]
            lsi_d_exp = [hh for hh, flag in zip(ddata["sIntDOWN"], lflags) if flag]
            
            #li_u_exp = ddata["IntUP"]
            #lsi_u_exp = ddata["sIntUP"]
            #li_d_exp = ddata["IntDOWN"]
            #lsi_d_exp = ddata["sIntDOWN"]
            ldiff_exp = [hh1-hh2 for hh1, hh2 in zip(li_u_exp, li_d_exp)]
            lsdiff_exp = [(hh1**2+hh2**2)**0.5  for hh1, hh2 in zip(lsi_u_exp, lsi_d_exp)]
            lcolors=[lhelp[iddata] for hh in range(3)]
            self.pcanvas.plot(lx,[li_u_exp,lsi_u_exp,None],[li_d_exp,lsi_d_exp,None],[ldiff_exp,lsdiff_exp,None],None,lcolors)
        xmin = float(self.le_xmin.text())
        xmax = float(self.le_xmax.text())
        ymin = float(self.le_ymin.text())
        ymax = float(self.le_ymax.text())
        self.pcanvas.set_limits(xmin,xmax,ymin,ymax)          
        self.pcanvas.draw()
"""    
