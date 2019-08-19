__author__ = 'ikibalin'
__version__ = "2019_04_06"
import os
import numpy

from neupy.f_common.cl_variable import Variable
from neupy.f_interface.cl_abstract_fract import AbstractFract


class Fract(AbstractFract):
    """
    Fract of atom_site(s) in unit cell.
    """
    def __init__(self, x = 0., y = 0., z = 0.):
        super(Fract, self).__init__()
        self._p_x = None
        self._p_y = None
        self._p_z = None
        
        self._refresh(x, y, z)
        self.set_val()


    def __repr__(self):
        lsout = """Fract: \n xyz: {:} {:} {:}""".format(self._p_x, self._p_y, 
                                  self._p_z)
        return lsout


    def _refresh(self, x, y, z):
        if not(isinstance(x, type(None))):
            self._p_x = numpy.mod(x, 1.)
        if not(isinstance(y, type(None))):
            self._p_y = numpy.mod(y, 1.)
        if not(isinstance(z, type(None))):
            self._p_z = numpy.mod(z, 1.)

            
    def set_val(self, x = None, y = None, z = None):
        self._refresh(x, y, z)
        
    def get_val(self, label):
        lab = "_p_"+label
        
        if lab in self.__dict__.keys():
            val = self.__dict__[lab]
            if isinstance(val, type(None)):
                self.set_val()
                val = self.__dict__[lab]
        else:
            print("The value '{:}' is not found".format(lab))
            val = None
        return val

    def list_vals(self):
        """
        give a list of parameters with small descripition
        """
        lsout = """
Parameters:
x, y, z is atoms coordinate
        """
        print(lsout)
    
    def calc_phase(self, space_group, h, k, l):
        """
        calculate phase: exp(-2 pi i * (h*x+k*y+l*z))
        r_11, r_22, r_33, r_12, r_13, r_23 are element of symmetry 
        """
        
        x, y, z = self._p_x, self._p_y, self._p_z

        r_11, r_12 = space_group.get_val("r_11"), space_group.get_val("r_12")
        r_13, r_21 = space_group.get_val("r_13"), space_group.get_val("r_21")
        r_22, r_23 = space_group.get_val("r_22"), space_group.get_val("r_23")
        r_31, r_32 = space_group.get_val("r_31"), space_group.get_val("r_32")
        r_33 = space_group.get_val("r_33")
        b_1, b_2 = space_group.get_val("b_1"), space_group.get_val("b_2")
        b_3 = space_group.get_val("b_3")
        
        np_h, np_x, np_r_11 = numpy.meshgrid(h, x, r_11, indexing="ij")
        np_k, np_y, np_r_22 = numpy.meshgrid(k, y, r_22, indexing="ij")
        np_l, np_z, np_r_33 = numpy.meshgrid(l, z, r_33, indexing="ij")
        
        np_r_12 = numpy.meshgrid(h, x, r_12, indexing="ij")[2]
        np_r_13 = numpy.meshgrid(k, y, r_13, indexing="ij")[2]
        np_r_23 = numpy.meshgrid(l, z, r_23, indexing="ij")[2]
        np_r_21 = numpy.meshgrid(h, x, r_21, indexing="ij")[2]
        np_r_31 = numpy.meshgrid(k, y, r_31, indexing="ij")[2]
        np_r_32 = numpy.meshgrid(l, z, r_32, indexing="ij")[2]

        np_b_1 = numpy.meshgrid(l, z, b_1, indexing="ij")[2]
        np_b_2 = numpy.meshgrid(l, z, b_2, indexing="ij")[2]
        np_b_3 = numpy.meshgrid(l, z, b_3, indexing="ij")[2]
        
        
        np_x_s = np_x*np_r_11 + np_y*np_r_12 + np_z*np_r_13 + np_b_1
        np_y_s = np_x*np_r_21 + np_y*np_r_22 + np_z*np_r_23 + np_b_2
        np_z_s = np_x*np_r_31 + np_y*np_r_32 + np_z*np_r_33 + np_b_3
        
        phase = numpy.exp(2*numpy.pi*1j*(np_h*np_x_s + np_k*np_y_s+ np_l*np_z_s))
        
        return phase
        
        
    def els4pos(self, space_group):
        """
        give the lelements of symmetry which transfer atom to the same atom
        """
        
        lelsymm = space_group.get_val("el_symm")
        lorig = space_group.get_val("orig")
        centr = space_group.get_val("centr")
        pcentr = space_group.get_val("pcentr")
    
        lelsat = []
        lelsuniqat, lcoorduniqat = [], []
        [x, y, z] = xyz
        x, y, z = x%1, y%1, z%1
        for els in lelsymm:
            for orig in lorig:
                xat = (els[0] + els[1]*x + els[ 2]*y + els[ 3]*z+orig[0])%1
                yat = (els[4] + els[5]*x + els[ 6]*y + els[ 7]*z+orig[1])%1
                zat = (els[8] + els[9]*x + els[10]*y + els[11]*z+orig[2])%1
                elsn = [els[0]+orig[0],els[1],els[2],els[3],els[4]+orig[1],els[5],els[6],els[7],
                els[8]+orig[2],els[9],els[10],els[11]]
                if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): lelsat.append(elsn)
                xyzatu = (round(xat,4),round(yat,4),round(zat,4))
                if (not(xyzatu in lcoorduniqat)):
                    lcoorduniqat.append(xyzatu)
                    lelsuniqat.append(elsn)
                if (centr):
                    elsn=[2*pcentr[0]-els[0]-orig[0],-1*els[1],-1*els[2],-1*els[3],
                        2*pcentr[1]-els[4]-orig[1],-1*els[5],-1*els[6],-1*els[7],
                        2*pcentr[2]-els[8]-orig[2],-1*els[9],-1*els[10],-1*els[11]]
                    xat,yat,zat=(2*pcentr[0]-xat)%1,(2*pcentr[1]-yat)%1,(2*pcentr[2]-zat)%1
                    if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): lelsat.append(elsn)
                    xyzatu=(round(xat,4),round(yat,4),round(zat,4))
                    if (not(xyzatu in lcoorduniqat)):
                        lcoorduniqat.append(xyzatu)
                        lelsuniqat.append(elsn)
        return lelsat,lelsuniqat

