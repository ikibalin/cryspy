"""
define classe PdProc for 1d powder diffraction experiment
"""
__author__ = 'ikibalin'
__version__ = "2019_09_06"
import os
import numpy


from pystar import Global


class PdProc(object):
    """
    This section contains the diffraction data set after processing
    and application of correction terms. If the data set is
    reprocessed, this section may be replaced (with the addition of
    a new _pd_block_id entry).

    Example:

    loop_
    _pd_proc_2theta
    _pd_proc_2theta_corrected
    _pd_proc_d_spacing
    _pd_proc_intensity_up_net
    _pd_proc_intensity_down_net
    _pd_proc_intensity_up_total
    _pd_proc_intensity_down_total
    _pd_proc_intensity_bkg_calc
    _pd_proc_intensity_up
    _pd_proc_intensity_up_sigma
    _pd_proc_intensity_down
    _pd_proc_intensity_down_sigma
     4.00  3.9  11.2   400.00   317.00   460.00   377.00   60.00    465.80000   128.97000   301.88000   129.30000

    reference: https://www.iucr.org/__data/iucr/cifdic_html/1/cif_pd.dic/Cpd_proc.html
    """
    def __init__(self, ttheta=[], ttheta_corrected=[], d=[], 
                up_net=[], down_net=[], up_total=[], down_total=[], bkg_calc=[], 
                up=[], down=[], up_sigma=[], down_sigma=[]):
        super(PdProc, self).__init__()
        self.__pd_proc_2theta = None
        self.__pd_proc_2theta_corrected = None
        self.__pd_proc_d_spacing = None
        self.__pd_proc_intensity_up_net = None
        self.__pd_proc_intensity_down_net = None
        self.__pd_proc_intensity_up_total = None
        self.__pd_proc_intensity_down_total = None
        self.__pd_proc_intensity_bkg_calc = None
        self.__pd_proc_intensity_up = None
        self.__pd_proc_intensity_up_sigma = None
        self.__pd_proc_intensity_down = None
        self.__pd_proc_intensity_down_sigma = None



        self.ttheta = ttheta
        self.ttheta_corrected = ttheta_corrected
        self.d = d
        self.up_net = up_net
        self.down_net = down_net
        self.up_total = up_total
        self.down_total = down_total
        self.bkg_calc = bkg_calc
        self.up = up
        self.up_sigma = up_sigma
        self.down = down
        self.down_sigma = down_sigma

    @property
    def ttheta(self):
        return self.__pd_proc_2theta
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
        self.__pd_proc_2theta = np_x_in

    @property
    def ttheta_corrected(self):
        return self.__pd_proc_2theta_corrected
    @ttheta_corrected.setter
    def ttheta_corrected(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_2theta_corrected = np_x_in


    @property
    def d(self):
        return self.__pd_proc_d_spacing
    @d.setter
    def d(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_d_spacing = np_x_in

    @property
    def up_net(self):
        return self.__pd_proc_intensity_up_net
    @up_net.setter
    def up_net(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_up_net = np_x_in


    @property
    def down_net(self):
        return self.__pd_proc_intensity_down_net
    @down_net.setter
    def down_net(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_down_net = np_x_in


    @property
    def up_total(self):
        return self.__pd_proc_intensity_up_total
    @up_total.setter
    def up_total(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_up_total = np_x_in


    @property
    def down_total(self):
        return self.__pd_proc_intensity_down_total
    @down_total.setter
    def down_total(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_down_total = np_x_in


    @property
    def bkg_calc(self):
        return self.__pd_proc_intensity_bkg_calc
    @bkg_calc.setter
    def bkg_calc(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_bkg_calc = np_x_in


    @property
    def up(self):
        return self.__pd_proc_intensity_up
    @up.setter
    def up(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_up = np_x_in

    @property
    def up_sigma(self):
        return self.__pd_proc_intensity_up_sigma
    @up_sigma.setter
    def up_sigma(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_up_sigma = np_x_in

    @property
    def down(self):
        return self.__pd_proc_intensity_down
    @down.setter
    def down(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_down = np_x_in

    @property
    def down_sigma(self):
        return self.__pd_proc_intensity_down_sigma
    @down_sigma.setter
    def down_sigma(self, l_x):
        l_x_in = []
        for x in l_x:
            if isinstance(x, float):
                x_in = x
            else:
                x_in = float(x)
            l_x_in.append(x_in)
        np_x_in = numpy.array(l_x_in, dtype=float)
        self.__pd_proc_intensity_down_sigma = np_x_in


            
    def __repr__(self):
        ls_out = ["PdProc:"]
        ls_out.append(" ttheta  ttheta_cor. up_total           up   down_total         down")
        for _1, _2, _3, _4, _5, _6 in zip(self.ttheta, self.ttheta_corrected, self.up_total, self.up, self.down_total, self.down):
            ls_out.append("{:7.2f} {:7.2f} {:12.5f} {:12.5f} {:12.5f} {:12.5f}".format(_1, _2, _3, _4, _5, _6))
        return "\n".join(ls_out)

    @property
    def to_cif(self):
        ls_out = []
        if self.is_defined:
            ls_out.append("loop_")
            ls_out.append("_pd_proc_2theta")
            ls_out.append("_pd_proc_2theta_corrected")
            #ls_out.append("_pd_proc_d_spacing")
            ls_out.append("_pd_proc_intensity_up_net")
            ls_out.append("_pd_proc_intensity_down_net")
            ls_out.append("_pd_proc_intensity_up_total")
            ls_out.append("_pd_proc_intensity_down_total")
            ls_out.append("_pd_proc_intensity_bkg_calc")
            ls_out.append("_pd_proc_intensity_up")
            ls_out.append("_pd_proc_intensity_up_sigma")
            ls_out.append("_pd_proc_intensity_down")
            ls_out.append("_pd_proc_intensity_down_sigma")
            #for _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12 in zip(self.ttheta, self.ttheta_corrected, self.d, 
            #    self.up_net, self.down_net, self.up_total, self.down_total, self.bkg_calc, 
            #    self.up, self.up_sigma, self.down, self.down_sigma):
            #    ls_out.append("{:} {:} {:} {:} {:} {:} {:} {:} {:} {:} {:} {:}".format(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12))
            for _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11 in zip(self.ttheta, self.ttheta_corrected, 
                self.up_net, self.down_net, self.up_total, self.down_total, self.bkg_calc, 
                self.up, self.up_sigma, self.down, self.down_sigma):
                ls_out.append("{:.3f} {:.3f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f} {:.5f}".format(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11))
        return "\n".join(ls_out)

    def from_cif(self, string: str):
        cif_global = Global()
        flag = cif_global.take_from_string(string)
        if not flag:
            return False
        flag = False
        flag = cif_global.is_prefix("_pd_proc")
        if flag:
            cif_loop = cif_global["_pd_proc"]
            l_name = cif_loop.names
            if "_pd_proc_2theta" in l_name:
                self.ttheta = [float(_1) for _1 in cif_loop["_pd_proc_2theta"]]
            if "_pd_proc_2theta_corrected" in l_name:
                self.ttheta_corrected = [float(_1) for _1 in cif_loop["_pd_proc_2theta_corrected"]]
            if "_pd_proc_d_spacing" in l_name:
                self.d = [float(_1) for _1 in cif_loop["_pd_proc_d_spacing"]]
            if "_pd_proc_intensity_up_net" in l_name:
                self.up_net = [float(_1) for _1 in cif_loop["_pd_proc_intensity_up_net"]]
            if "_pd_proc_intensity_down_net" in l_name:
                self.down_net = [float(_1) for _1 in cif_loop["_pd_proc_intensity_down_net"]]
            if "_pd_proc_intensity_up_total" in l_name:
                self.up_total = [float(_1) for _1 in cif_loop["_pd_proc_intensity_up_total"]]
            if "_pd_proc_intensity_down" in l_name:
                self.down_total = [float(_1) for _1 in cif_loop["_pd_proc_intensity_down_total"]]
            if "_pd_proc_intensity_bkg_calc" in l_name:
                self.bkg_calc = [float(_1) for _1 in cif_loop["_pd_proc_intensity_bkg_calc"]]
            if "_pd_proc_intensity_up" in l_name:
                self.up = [float(_1) for _1 in cif_loop["_pd_proc_intensity_up"]]
            if "_pd_proc_intensity_up_sigma" in l_name:
                self.up_sigma = [float(_1) for _1 in cif_loop["_pd_proc_intensity_up_sigma"]]
            if "_pd_proc_intensity_down" in l_name:
                self.down = [float(_1) for _1 in cif_loop["_pd_proc_intensity_down"]]
            if "_pd_proc_intensity_down_sigma" in l_name:
                self.down_sigma = [float(_1) for _1 in cif_loop["_pd_proc_intensity_down_sigma"]]
        else:
            self.ttheta, self.ttheta_corrected, self.d = [], [], []
            self.up_net, self.down_net, self.up_total, self.down_total, self.bkg_calc = [], [], [], [], []
            self.up, self.up_sigma, self.down, self.down_sigma = [], [], [], []
        return True

    @property
    def is_defined(self):
        cond = all([self.ttheta is not None, self.ttheta_corrected is not None, self.d is not None, 
                    self.up_net is not None, self.down_net is not None, self.up_total is not None, self.down_total is not None,
                    self.bkg_calc is not None,
                    self.up is not None, self.up_sigma is not None, self.down is not None, self.down_sigma is not None])
        return cond

    @property
    def is_variable(self):
        return False
    
    def get_variables(self):
        return []

    def _show_message(self, s_out: str):
        print("***  Error ***")
        print(s_out)
