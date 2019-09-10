"""
define class Pd2dMeas which describes the 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_10"
import os
import numpy


from pystar import CIFglobal


class Pd2dMeas(object):
    """
    This section contains the measured diffractogram and information
    about the conditions used for the measurement of the diffraction 
    data set, prior to processing and application of correction
    terms.

    Example:

    _pd2d_meas_2theta_phi_intensity_up
    ;
    ;

    _pd2d_meas_2theta_phi_intensity_up_sigma
    ;
    ;

    _pd2d_meas_2theta_phi_intensity_down
    ;
    ;

    _pd2d_meas_2theta_phi_intensity_down_sigma
    ;
    ;
    """
    def __init__(self, ttheta=[], phi=[], up=[[]], down=[[]], up_sigma=[[]], down_sigma=[[]]
                 ):
        super(Pd2dMeas, self).__init__()
        self.__pd2d_meas_angle_2theta = None
        self.__pd2d_meas_angle_phi = None
        self.__pd2d_meas_intensity_up = None
        self.__pd2d_meas_intensity_up_sigma = None
        self.__pd2d_meas_intensity_down = None
        self.__pd2d_meas_intensity_down_sigma = None

        self.phi = phi
        self.ttheta = ttheta
        self.up = up
        self.up_sigma = up_sigma
        self.down = down
        self.down_sigma = down_sigma

    @property
    def ttheta(self):
        return self.__pd2d_meas_angle_2theta
    @ttheta.setter
    def ttheta(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_meas_angle_2theta = np_x_in

    @property
    def phi(self):
        return self.__pd2d_meas_angle_phi
    @phi.setter
    def phi(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd2d_meas_angle_phi = np_x_in

    @property
    def up(self):
        return self.__pd2d_meas_intensity_up
    @up.setter
    def up(self, ll_x):
        ll_x_in = []
        for l_x in ll_x:
            l_x_in = []
            for x in l_x:
                if (isinstance(x, float) | (x is None)):
                    x_in = x
                else:
                    x_in = float(x)
                l_x_in.append(x_in)
            ll_x_in.append(l_x_in)
        np_x_in = numpy.array(ll_x_in, dtype=float)
        self.__pd2d_meas_intensity_up = np_x_in

    @property
    def up_sigma(self):
        return self.__pd2d_meas_intensity_up_sigma
    @up_sigma.setter
    def up_sigma(self, ll_x):
        ll_x_in = []
        for l_x in ll_x:
            l_x_in = []
            for x in l_x:
                if (isinstance(x, float) | (x is None)):
                    x_in = x
                else:
                    x_in = float(x)
                l_x_in.append(x_in)
            ll_x_in.append(l_x_in)
        np_x_in = numpy.array(ll_x_in, dtype=float)
        self.__pd2d_meas_intensity_up_sigma = np_x_in

    @property
    def down(self):
        return self.__pd2d_meas_intensity_down
    @down.setter
    def down(self, ll_x):
        ll_x_in = []
        for l_x in ll_x:
            l_x_in = []
            for x in l_x:
                if (isinstance(x, float) | (x is None)):
                    x_in = x
                else:
                    x_in = float(x)
                l_x_in.append(x_in)
            ll_x_in.append(l_x_in)
        np_x_in = numpy.array(ll_x_in, dtype=float)
        self.__pd2d_meas_intensity_down = np_x_in

    @property
    def down_sigma(self):
        return self.__pd2d_meas_intensity_down_sigma
    @down_sigma.setter
    def down_sigma(self, ll_x):
        ll_x_in = []
        for l_x in ll_x:
            l_x_in = []
            for x in l_x:
                if (isinstance(x, float) | (x is None)):
                    x_in = x
                else:
                    x_in = float(x)
                l_x_in.append(x_in)
            ll_x_in.append(l_x_in)
        np_x_in = numpy.array(ll_x_in, dtype=float)
        self.__pd2d_meas_intensity_down_sigma = np_x_in
            
    def __repr__(self):
        ls_out = ["Pd2dMeas:"]
        ls_out.append("\n ttheta: {:}".format(str(self.ttheta)))
        ls_out.append(" phi: {:}".format(str(self.phi)))
        ls_out.append("\n up")
        ls_out.append(str(self.up))
        ls_out.append("\n up_sigma")
        ls_out.append(str(self.up_sigma))
        ls_out.append("\n down")
        ls_out.append(str(self.down))
        ls_out.append("\n down_sigma")
        ls_out.append(str(self.down_sigma))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("_pd2d_meas_2theta_phi_intensity_up")
            ls_out.append(";")
            ls_out.append("{:12} ".format(len(self.ttheta)) + " ".join(["{:6.2f}      ".format(_) for _ in self.phi]))
            for ttheta, l_intensity in zip(self.ttheta, self.up):
                ls_out.append("{:12.2f} ".format(ttheta) + " ".join(["{:12}".format(_) if _ is not numpy.nan else "        None" for _ in l_intensity]))
            ls_out.append(";")

            ls_out.append("\n_pd2d_meas_2theta_phi_intensity_up_sigma")
            ls_out.append(";")
            ls_out.append("{:12} ".format(len(self.ttheta)) + " ".join(["{:6.2f}      ".format(_) for _ in self.phi]))
            for ttheta, l_intensity in zip(self.ttheta, self.up_sigma):
                ls_out.append("{:12.2f} ".format(ttheta) + " ".join(["{:12}".format(_) if _ is not numpy.nan else "        None" for _ in l_intensity]))
            ls_out.append(";")

            ls_out.append("\n_pd2d_meas_2theta_phi_intensity_down")
            ls_out.append(";")
            ls_out.append("{:12} ".format(len(self.ttheta)) + " ".join(["{:6.2f}      ".format(_) for _ in self.phi]))
            for ttheta, l_intensity in zip(self.ttheta, self.down):
                ls_out.append("{:12.2f} ".format(ttheta) + " ".join(["{:12}".format(_) if _ is not numpy.nan else "        None" for _ in l_intensity]))
            ls_out.append(";")

            ls_out.append("\n_pd2d_meas_2theta_phi_intensity_down_sigma")
            ls_out.append(";")
            ls_out.append("{:12} ".format(len(self.ttheta)) + " ".join(["{:6.2f}      ".format(_) for _ in self.phi]))
            for ttheta, l_intensity in zip(self.ttheta, self.down_sigma):
                ls_out.append("{:12.2f} ".format(ttheta) + " ".join(["{:12}".format(_) if _ is not numpy.nan else "        None" for _ in l_intensity]))
            ls_out.append(";")
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = CIFglobal()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = cif_global.is_prefix("_pd2d_meas_2theta_phi_intensity_up")
        if flag:
            cif_value = cif_global["_pd2d_meas_2theta_phi_intensity_up"]
            string = cif_value.value
            l_1 = string.strip().split("\n")
            l_phi = [float(_) for _ in l_1[0].strip().split()[1:]]
            l_ttheta, ll_intensity = [], []
            for line in l_1[1:]:
                l_1 = line.strip().split()
                l_ttheta.append(float(l_1[0]))
                ll_intensity.append([float(_) if _ != "None" else None for _ in l_1[1:]])
            self.phi = l_phi
            self.ttheta = l_ttheta
            self.up = ll_intensity
        flag = cif_global.is_prefix("_pd2d_meas_2theta_phi_intensity_up_sigma")
        if flag:
            cif_value = cif_global["_pd2d_meas_2theta_phi_intensity_up_sigma"]
            string = cif_value.value
            l_1 = string.strip().split("\n")
            l_phi = [float(_) for _ in l_1[0].strip().split()[1:]]
            l_ttheta, ll_intensity = [], []
            for line in l_1[1:]:
                l_1 = line.strip().split()
                l_ttheta.append(float(l_1[0]))
                ll_intensity.append([float(_) if _ != "None" else None for _ in l_1[1:]])
            self.phi = l_phi
            self.ttheta = l_ttheta
            self.up_sigma = ll_intensity
        flag = cif_global.is_prefix("_pd2d_meas_2theta_phi_intensity_down")
        if flag:
            cif_value = cif_global["_pd2d_meas_2theta_phi_intensity_down"]
            string = cif_value.value
            l_1 = string.strip().split("\n")
            l_phi = [float(_) for _ in l_1[0].strip().split()[1:]]
            l_ttheta, ll_intensity = [], []
            for line in l_1[1:]:
                l_1 = line.strip().split()
                l_ttheta.append(float(l_1[0]))
                ll_intensity.append([float(_) if _ != "None" else None for _ in l_1[1:]])
            self.phi = l_phi
            self.ttheta = l_ttheta
            self.down = ll_intensity
        flag = cif_global.is_prefix("_pd2d_meas_2theta_phi_intensity_down_sigma")
        if flag:
            cif_value = cif_global["_pd2d_meas_2theta_phi_intensity_down_sigma"]
            string = cif_value.value
            l_1 = string.strip().split("\n")
            l_phi = [float(_) for _ in l_1[0].strip().split()[1:]]
            l_ttheta, ll_intensity = [], []
            for line in l_1[1:]:
                l_1 = line.strip().split()
                l_ttheta.append(float(l_1[0]))
                ll_intensity.append([float(_) if _ != "None" else None for _ in l_1[1:]])
            self.phi = l_phi
            self.ttheta = l_ttheta
            self.down_sigma = ll_intensity
        return True

    @property
    def is_defined(self):
        cond = all([self.phi is not None, self.ttheta is not None, 
                    self.up is not None, self.up_sigma is not None, 
                    self.down is not None, self.down_sigma is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
