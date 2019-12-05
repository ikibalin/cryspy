"""
internal class to calculate phases
"""
__author__ = 'ikibalin'
__version__ = "2019_09_02"
import os
import numpy


class Fract(object):
    """
    Fract of atom_sites in unit cell.
    """
    def __init__(self, x=numpy.array([0], dtype=float),
                       y=numpy.array([0], dtype=float),
                       z=numpy.array([0], dtype=float)):
        super(Fract, self).__init__()
        self.__atom_site_fract_x = None
        self.__atom_site_fract_y = None
        self.__atom_site_fract_z = None
        self.x = x 
        self.y = y
        self.z = z 

    def _trans_to_float_array(self, x):
        if isinstance(x, numpy.ndarray):
            x_out = x.astype(float)
        else:
            x_out = numpy.array([x], dtype=float)
        res = numpy.mod(x_out, 1.)
        return res

    @property
    def x(self):
        """

        reference:
        """
        return self.__atom_site_fract_x
    @x.setter
    def x(self, x):
        self.__atom_site_fract_x = self._trans_to_float_array(x)
    @property
    def y(self):
        """

        reference:
        """
        return self.__atom_site_fract_y
    @y.setter
    def y(self, x):
        self.__atom_site_fract_y = self._trans_to_float_array(x)
    @property
    def z(self):
        """

        reference:
        """
        return self.__atom_site_fract_z
    @z.setter
    def z(self, x):
        self.__atom_site_fract_z = self._trans_to_float_array(x)


    def __repr__(self):
        ls_out = ["Fract:\n       x        y        z"]
        ls_out.extend(["{:8.5f} {:8.5f} {:8.5f}".format(hh_1, hh_2, hh_3) \
                       for hh_1, hh_2, hh_3 in zip(self.x, self.y, self.z)])
        return "\n".join(ls_out)

    
    def calc_phase(self, space_group, h, k, l):
        """
        calculate phase: exp(-2 pi i * (h*x+k*y+l*z))
        r_11, r_22, r_33, r_12, r_13, r_23 are element of symmetry 
        """
        
        x, y, z = self.x, self.y, self.z

        r_11, r_12 = space_group.r_11, space_group.r_12
        r_13, r_21 = space_group.r_13, space_group.r_21
        r_22, r_23 = space_group.r_22, space_group.r_23
        r_31, r_32 = space_group.r_31, space_group.r_32
        r_33 = space_group.r_33
        b_1, b_2 = space_group.b_1, space_group.b_2
        b_3 = space_group.b_3
        
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
        give the l_element of symmetry which transfer atom to the same atom
        """
        
        l_elsymm = space_group.el_symm
        l_orig = space_group.orig
        centr = space_group.centr
        pcentr = space_group.pcentr
    
        l_elsat = []
        l_els_uniq_at, l_coord_uniq_at = [], []
        [x, y, z] = xyz
        x, y, z = x%1, y%1, z%1
        for els in l_elsymm:
            for orig in l_orig:
                xat = (els[0] + els[1]*x + els[ 2]*y + els[ 3]*z+orig[0])%1
                yat = (els[4] + els[5]*x + els[ 6]*y + els[ 7]*z+orig[1])%1
                zat = (els[8] + els[9]*x + els[10]*y + els[11]*z+orig[2])%1
                elsn = [els[0]+orig[0],els[1],els[2],els[3],els[4]+orig[1],els[5],els[6],els[7],
                els[8]+orig[2],els[9],els[10],els[11]]
                if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): l_elsat.append(elsn)
                xyzatu = (round(xat,4),round(yat,4),round(zat,4))
                if (not(xyzatu in l_coord_uniq_at)):
                    l_coord_uniq_at.append(xyzatu)
                    l_els_uniq_at.append(elsn)
                if (centr):
                    elsn=[2*pcentr[0]-els[0]-orig[0],-1*els[1],-1*els[2],-1*els[3],
                        2*pcentr[1]-els[4]-orig[1],-1*els[5],-1*els[6],-1*els[7],
                        2*pcentr[2]-els[8]-orig[2],-1*els[9],-1*els[10],-1*els[11]]
                    xat,yat,zat=(2*pcentr[0]-xat)%1,(2*pcentr[1]-yat)%1,(2*pcentr[2]-zat)%1
                    if ((abs(xat-x)<10**-5)&(abs(yat-y)<10**-5)&(abs(zat-z)<10**-5)): l_elsat.append(elsn)
                    xyzatu=(round(xat,4),round(yat,4),round(zat,4))
                    if (not(xyzatu in l_coord_uniq_at)):
                        l_coord_uniq_at.append(xyzatu)
                        l_els_uniq_at.append(elsn)
        return l_elsat,l_els_uniq_at

