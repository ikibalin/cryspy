"""
The program performs integration over Debye ring with different limit conditions
"""
__author__ = 'ikibalin'
__version__ = "2019_04_19"
import os
import numpy
import sys



import matplotlib
import matplotlib.pyplot

 
def in_tag(lcont, lab):
    numb_b = lcont.find("<"+lab)
    numb_e = lcont.find("/"+lab+">")
    if ((numb_b != -1)|(numb_e != -1)):
        h_str = lcont[(numb_b+len(lab)+1):(numb_e-1)] 
        numb_b = h_str.find(">")
        h_str_atr = h_str[:numb_b]
        h_str_val = h_str[(numb_b+1):]
        return h_str_atr, h_str_val
    return None, None

 
def save_matrix_to_str(np_2d, ax_1, ax_2):
    "Save matrix to string"
    lsout=[]
    line=["{:15}".format(ax_2.size)]
    line.extend(["{0:15.5f}".format(hh) for hh in ax_1])
    lsout.append("".join(line))
    for np_1d, hh in zip(np_2d.transpose(), ax_2):
        line=["{0:15.5f}".format(hh)]
        for val in np_1d:
            if numpy.isnan(val):
                line.append("           None")
            else:
                line.append("{0:15.5f}".format(val))
        lsout.append("".join(line))
    return lsout

class XML(object):
    """
    Class to describe container for XML data file
    """
    def __init__(self, manip = "unknown", wavelength=None, 
                 norm_monitor=1000000, flag_monitor=True, setup=None):
        super(XML,self).__init__()
        self._p_manip = None
        self._p_wavelength = None
        self._p_norm_monitor = None
        self._p_flag_monitor = None
        self._p_file_format = None
        self._p_setup=None
        self._list_frame = []
        self._refresh(manip, wavelength,
                      norm_monitor, flag_monitor, setup)
     
    def give_list_frame(self):
        return self._list_frame
     
    def __add__(self, var2):
        """
        output is XML
        """
        manip_1 = self._p_manip
        wavelength_1 = self._p_wavelength
        norm_monitor = self._p_norm_monitor
        flag_monitor = self._p_flag_monitor
        setup = self._p_setup
        l_frame_1 = self.give_list_frame()
        nl_1 = len(l_frame_1)

        manip_2 = var2.get_val("manip")
        wavelength_2 = var2.get_val("wavelength")
        l_frame_2 = var2.give_list_frame()
        nl_2 = len(l_frame_2)
        
        manip = None
        if ((manip_1 == manip_2)&(wavelength_1 == wavelength_2)&(nl_1 == nl_2)):
            manip = manip_1
            wavelength = wavelength_1
        else:
            print("Tho xml files can not be summarized")
            return
        xml = XML(manip=manip, wavelength=wavelength, 
                  norm_monitor=norm_monitor, flag_monitor=flag_monitor, 
                  setup=setup)
        for frame_1, frame_2 in zip(l_frame_1, l_frame_2):
            frame = frame_1+frame_2 
            xml.add_frame(frame)
        return xml
     
    def __radd__(self, var2):
        """
        output is XML
        """
        manip_1 = var2.get_val("manip")
        wavelength_1 = var2.get_val("wavelength")
        norm_monitor = var2.get_val("norm_monitor")
        flag_monitor = var2.get_val("flag_monitor")
        setup = var2.get_val("setup")
        l_frame_1 = self.give_list_frame()
        nl_1 = len(l_frame_1)

        manip_2 = self._p_manip
        wavelength_2 = self._p_wavelength
        l_frame_2 = self.give_list_frame()
        nl_2 = len(l_frame_2)
        
        manip = None
        if ((manip_1 == manip_2)&(wavelength_1 == wavelength_2)&(nl_1 == nl_2)):
            manip = manip_1
            wavelength = wavelength_1
        else:
            print("Tho xml files can not be summarized")
            return
        
        xml = XML(manip=manip, wavelength=wavelength, 
                  norm_monitor=norm_monitor, flag_monitor=flag_monitor, 
                  setup=setup)
        for frame_1, frame_2 in zip(l_frame_1, l_frame_2):
            frame = frame_1+frame_2 
            xml.add_frame(frame)
        return xml
         
    def __repr__(self):
        lsout = """XML (file format '{:}'): \n wavelength: {:},
 manip: {:}, flag_monitor: {:}, norm_monitor: {:},
 number of frames: {:}""".format(self._p_file_format, self._p_wavelength, 
                 self._p_manip, self._p_flag_monitor, self._p_norm_monitor, 
                 len(self._list_frame))
        return lsout
     
    def set_val(self, manip=None, wavelength=None, norm_monitor=None, 
                flag_monitor=None, setup=None):
        if not(isinstance(wavelength, type(None))):
            wavelength = float(wavelength)
        self._refresh(manip, wavelength, 
                      norm_monitor, flag_monitor, setup)
     
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
 setup_name, 
 nx_min, nx_max, 
 ny_min, ny_max
        """
        print(lsout )
     
    def _refresh(self, manip, wavelength, 
                 norm_monitor, flag_monitor, setup):
        if not(isinstance(manip, type(None))):
            self._p_manip = manip
        if not(isinstance(wavelength, type(None))):
            self._p_wavelength = wavelength
        if not(isinstance(norm_monitor, type(None))):
            self._p_norm_monitor = norm_monitor
        if not(isinstance(flag_monitor, type(None))):
            self._p_flag_monitor = flag_monitor
        if not(isinstance(setup, type(None))):
            self._p_setup = setup


     
    def _define_format(self, lcont):
        """
        Define format of xml file
        values: '5C1_old', '5C1_new', '6T2_old', '6T2_new', 'unknown'
        """
        file_format = 'unknown'
        
        lab ="Frame" 
        h_str_atr, h_str_val = in_tag(lcont, lab)
        if "MagneticField" in h_str_atr:
            file_format = '6T2_old'
        elif "monitorUpCount" in h_str_val:
            file_format = '6T2_new'
        elif "flipper_temp" in h_str_val:
            file_format = '5C1_new'
        self._p_file_format = file_format
        return file_format
         
    def read_file(self, f_name):
        """
        Read information about object from file
        """
        fid = open(f_name, 'r')
        lcont = fid.read()
        fid.close()
        
        self._list_frame = []
        
        
        file_format = self._define_format(lcont)    
        
        if file_format == "5C1_new":
            self._read_5C1_6T2_new(lcont)
        elif file_format == "6T2_old":
            self._read_6T2_old(lcont)
        elif file_format == "6T2_new":
            self._read_5C1_6T2_new(lcont)
        else:
            print("File format is unknown")
            return
        manip = self.get_val("manip")
        setup = Setup(manip)
        self._p_setup = setup
            
         
    def _read_6T2_old(self, lcont):
        self.set_val(manip="6T2")

        flag_frame = True
        lab ="Frame" 
        numb = 0
        while flag_frame: 
            h_str_atr, h_str_val = in_tag(lcont[numb:], lab)
            if h_str_val is None:
                flag_frame = False
            else:
                val = lcont[numb:].find(lab+">")
                numb += val+len(lab)+1
                frame = Frame()
                frame.read_str_6T2_old(h_str_atr, h_str_val)
                self._list_frame.append(frame)
    
     
    def _read_5C1_6T2_new(self, lcont):
        
        lab ="manip" 
        h_str_atr, s_manip = in_tag(lcont, lab)
        self.set_val(manip=s_manip)

        lab ="wavelenght" 
        h_str_atr, h_str_val = in_tag(lcont, lab)
        self.set_val(wavelength=h_str_val)

        lab ="wavelength" 
        h_str_atr, h_str_val = in_tag(lcont, lab)
        self.set_val(wavelength=h_str_val)

            
        flag_frame = True
        lab ="Frame" 
        numb = 0
        while flag_frame: 
            h_str_atr, h_str_val = in_tag(lcont[numb:], lab)
            if h_str_val is None:
                flag_frame = False
            else:
                val = lcont[numb:].find(lab+">")
                numb += val+len(lab)+1
                frame = Frame()
                frame.read_str_5C1_6T2_new(h_str_atr, h_str_val, s_manip)
                self._list_frame.append(frame)
         
    def add_frame(self, frame):
        self._list_frame.append(frame)
         
    def print_frame_info(self):
        ls_out = ["Number of frames is {:}".format(len(self._list_frame))]
        line = "Frame  temp.  field      phi  omega  gamma    time monitor     i_sum    i_diff   rel.%"
        ls_out.append(line)
        for iframe, frame in enumerate(self._list_frame):
            t = frame.get_val("temperature")
            f = frame.get_val("field")
            p = frame.get_val("phi")
            o = frame.get_val("omega")
            g = frame.get_val("gamma")
            time = frame.get_val("time_count")
            m = frame.get_val("monitor")
            s_i_u = frame.get_val("int_u").sum()
            s_i_d = frame.get_val("int_d").sum()
            s_i_sum = s_i_u + s_i_d
            s_i_dif = s_i_u - s_i_d

            line=" {:4} {:6} {:6}   {:6} {:6} {:6} {:7} {:7} {:9} {:9} {:7.1f}".format(
                    iframe+1, t, f, p, o, g, time, m, s_i_sum, s_i_dif, 100.*float(s_i_dif)*1./float(s_i_sum))
            ls_out.append(line)
        print("\n".join(ls_out))
        
         
    def give_t_p_g_n_iu_siu_id_sid(self):
        """
        it is only for powder (there is no phi and omega)
        """
        setup = self._p_setup
        if setup is None:
            print("Setup is not defined")
            return
        mask = setup.get_val("mask")

        l_frame = self._list_frame
        np_g = numpy.array([], dtype = float)
        np_n = numpy.array([], dtype = float)
        np_t = numpy.array([], dtype = float)
        np_p = numpy.array([], dtype = float)
        np_i_u = numpy.array([], dtype = float)
        np_si_u = numpy.array([], dtype = float)        
        np_i_d = numpy.array([], dtype = float)
        np_si_d = numpy.array([], dtype = float)    

        
        for frame in l_frame:
            gamma = frame.get_val("gamma")*numpy.pi/180.
            t_2d, p_2d, g_2d, n_2d = setup.calc_tpgn(gamma)
            i_u_2d = frame.get_val("int_u")
            i_d_2d = frame.get_val("int_d")
            #g_1d = numpy.ravel(g_2d, order="C")
            g_1d = g_2d[mask]
            n_1d = n_2d[mask]
            t_1d = t_2d[mask]
            p_1d = p_2d[mask]
            i_u_1d = i_u_2d[mask]
            i_d_1d = i_d_2d[mask]
            if self._p_flag_monitor:
                coeff = frame.get_val("monitor")*1./self._p_norm_monitor
            else:
                coeff = 1.
            i_u_fl = coeff * i_u_1d
            i_d_fl = coeff * i_d_1d
            si_u_fl = coeff * numpy.sqrt(i_u_1d)
            si_d_fl = coeff * numpy.sqrt(i_d_1d)

            np_t = numpy.hstack([np_t, t_1d])
            np_p = numpy.hstack([np_p, p_1d])
            np_g = numpy.hstack([np_g, g_1d])
            np_n = numpy.hstack([np_n, n_1d])
            np_i_u = numpy.hstack([np_i_u, i_u_fl])
            np_si_u = numpy.hstack([np_si_u, si_u_fl])
            np_i_d = numpy.hstack([np_i_d, i_d_fl])
            np_si_d = numpy.hstack([np_si_d, si_d_fl])
        return  np_t, np_p, np_g, np_n, np_i_u, np_si_u, np_i_d, np_si_d


     
    def calc_in_grid_gamma_nu(self, g_min=None, g_max=None, g_step=None, n_min=None, n_max=None, n_step=None):
        """
        parameters are given in degrees
        """
        setup = self._p_setup
        name = setup._p_name
        if (name in ["5C1", "VIP"]):
            g_min_r = 4.*numpy.pi*1./180.
            g_max_r = 80.*numpy.pi*1./180. 
            g_step_r = 0.1*numpy.pi*1./180.
            n_min_r = -4.*numpy.pi*1./180.
            n_max_r = 20.*numpy.pi*1./180. 
            n_step_r = 0.1072*numpy.pi*1./180.
        elif (name in ["6T2"]):
            g_min_r = 4.*numpy.pi*1./180.
            g_max_r = 80.*numpy.pi*1./180. 
            g_step_r = 0.2*numpy.pi*1./180.
            n_min_r = -11.*numpy.pi*1./180.
            n_max_r = 15.*numpy.pi*1./180. 
            n_step_r = 0.3*numpy.pi*1./180.
        else:
            print("Setup '{:}' is unknown".format(name))
            return
            
        if (g_min is not None): g_min_r = g_min*numpy.pi*1./180.
        if (g_max is not None): g_max_r = g_max*numpy.pi*1./180.
        if (g_step is not None): g_step_r = g_step*numpy.pi*1./180.

        if (n_min is not None): n_min_r = n_min*numpy.pi*1./180.
        if (n_max is not None): n_max_r = n_max*numpy.pi*1./180.
        if (n_step is not None): n_step_r = n_step*numpy.pi*1./180.
            
        
        grid_gamma_nu = GridGammaNu(g_min=g_min_r, g_max=g_max_r, g_step=g_step_r, 
                              n_min=n_min_r, n_max=n_max_r, n_step=n_step_r)

        np_g, np_n, np_i_u, np_si_u, np_i_d, np_si_d = self.give_t_p_g_n_iu_siu_id_sid()[2:]
        grid_gamma_nu.load_to_grid(np_g, np_n, np_i_u, np_si_u, np_i_d, np_si_d) 
        
        return grid_gamma_nu
     
    def calc_in_grid_ttheta_phi(self, tth_min=None, tth_max=None, tth_step=None, 
                             phi_min=None, phi_max=None, phi_step=None):
        """
        parameters are given in degrees
        
        it is for test
        """
        setup = self._p_setup
        name = setup._p_name
        if (name in ["5C1", "VIP"]):
            tth_min_r = 4.*numpy.pi*1./180.
            tth_max_r = 80.*numpy.pi*1./180. 
            tth_step_r = 0.1*numpy.pi*1./180.
            phi_min_r = -10.*numpy.pi*1./180.
            phi_max_r = 70.*numpy.pi*1./180. 
            phi_step_r = 0.5*numpy.pi*1./180.
        elif (name in ["6T2"]):
            tth_min_r = 4.*numpy.pi*1./180.
            tth_max_r = 80.*numpy.pi*1./180. 
            tth_step_r = 0.2*numpy.pi*1./180.
            phi_min_r = -20.*numpy.pi*1./180.
            phi_max_r = 40.*numpy.pi*1./180. 
            phi_step_r = 0.5*numpy.pi*1./180.
        else:
            print("Setup '{:}' is unknown".format(name))
            return
            
        if (tth_min is not None): tth_min_r = tth_min*numpy.pi*1./180.
        if (tth_max is not None): tth_max_r = tth_max*numpy.pi*1./180.
        if (tth_step is not None): tth_step_r = tth_step*numpy.pi*1./180.

        if (phi_min is not None): phi_min_r = phi_min*numpy.pi*1./180.
        if (phi_max is not None): phi_max_r = phi_max*numpy.pi*1./180.
        if (phi_step is not None): phi_step_r = phi_step*numpy.pi*1./180.
            
        
        grid_ttheta_phi = GridTthetaPhi(tth_min=tth_min_r, tth_max=tth_max_r, 
                                        tth_step=tth_step_r, phi_min=phi_min_r, 
                                        phi_max=phi_max_r, phi_step=phi_step_r)

        np_t, np_p, np_g, np_n, np_i_u, np_si_u, np_i_d, np_si_d = self.give_t_p_g_n_iu_siu_id_sid()
        grid_ttheta_phi.load_to_grid(np_t, np_p, np_i_u, np_si_u, np_i_d, np_si_d) 
        return grid_ttheta_phi 


class Frame(object):
    """
    Class to describe information from position sensitive detector
    """
    def __init__(self, nx=1, ny=1, temperature=0, field=0, phi=0, omega=0, 
                 gamma=0, time_count=1, monitor=1, int_u=0., int_d=0.):
        super(Frame,self).__init__()
        self._p_nx = None
        self._p_ny = None
        self._p_temperature = None
        self._p_field = None
        self._p_phi = None
        self._p_omega = None
        self._p_gamma = None
        self._p_time_count = None
        self._p_monitor = None
        self._p_int_u = None
        self._p_int_d = None
        self._refresh(nx, ny, temperature, field, phi, omega, gamma, 
                      time_count, monitor, int_u, int_d)
         
    def __repr__(self):
        lsout = """Frame: \n nx: {:}, ny: {:}\n temperature: {:}, field: {:}
 phi: {:}, omega: {:}, gamma: {:}
 time_count: {:}, monitor: {:}""".format(self._p_nx, self._p_ny, 
                 self._p_temperature, self._p_field, self._p_phi,self._p_omega, 
                 self._p_gamma, self._p_time_count, self._p_monitor)
        return lsout
     
    def get_all_param(self):
        x_1 = self._p_nx
        y_1 = self._p_ny 
        t_1 = self._p_temperature
        f_1 = self._p_field
        ph_1 = self._p_phi
        om_1 = self._p_omega
        ga_1 = self._p_gamma
        time_1 = self._p_time_count
        mon_1 = self._p_monitor
        i_u_1 = self._p_int_u 
        i_d_1 = self._p_int_d 
        return x_1, y_1, t_1, f_1, ph_1, om_1, ga_1, time_1, mon_1, i_u_1, i_d_1
        
     
    def __add__(self, var2):
        """
        output is Frame
        """
        x_1, y_1, t_1, f_1, ph_1, om_1, ga_1, time_1, mon_1, i_u_1, i_d_1 = self.get_all_param()
        x_2, y_2, t_2, f_2, ph_2, om_2, ga_2, time_2, mon_2, i_u_2, i_d_2 = var2.get_all_param()

    
        nx, ny, temperature, field = None, None, None, None
        gamma, phi, omega = None, None, None,
        
        if ((x_1 != x_2)|(y_1 != y_2)):
            print("The size of frames is different. ")
            return None
        else:
            int_u = i_u_1 + i_u_2
            int_d = i_d_1 + i_d_2
            
        if x_1 == x_2: nx = x_1
        if y_1 == y_2: ny = y_1
        if abs(t_1-t_2)<1.: temperature = t_1
        if abs(f_1-f_2)<0.1: field = f_1
        if abs(ga_1-ga_2)<0.1: gamma = ga_1
        if abs(ph_1-ph_2)<0.1: phi = ph_1
        if abs(om_1-om_2)<0.1: omega = om_1
            
        monitor = mon_1 + mon_2
        time_count = time_1 + time_2

        frame = Frame(nx=nx, ny=ny, temperature=temperature, field=field, 
              phi=phi, omega=omega, gamma=gamma, time_count=time_count, 
              monitor=monitor, int_u=int_u, int_d=int_d)
        return frame
     
    def __radd__(self, var2):
        """
        output is float
        """
        x_2, y_2, t_2, f_2, ph_2, om_2, ga_2, time_2, mon_2, i_u_2, i_d_2 = self.get_all_param()
        x_1, y_1, t_1, f_1, ph_1, om_1, ga_1, time_1, mon_1, i_u_1, i_d_1 = var2.get_all_param()
    
        nx, ny, temperature, field = None, None, None, None
        gamma, phi, omega = None, None, None,
        
        if ((x_1 != x_2)|(y_1 != y_2)):
            print("The size of frames is different. ")
            return None
        else:
            int_u = i_u_1 + i_u_2
            int_d = i_d_1 + i_d_2
            
        if x_1 == x_2: nx = x_1
        if y_1 == y_2: ny = y_1
        if abs(t_1-t_2)<1.: temperature = t_1
        if f_1 == f_2: field = f_1
        if ga_1 == ga_2: gamma = ga_1
        if ph_1 == ph_2: phi = ph_1
        if om_1 == om_2: omega = om_1
            
        monitor = mon_1 + mon_2
        time_count = time_1 + time_2

        frame = Frame(nx=nx, ny=ny, temperature=temperature, field=field, 
              phi=phi, omega=omega, gamma=gamma, time_count=time_count, 
              monitor=monitor, int_u=int_u, int_d=int_d)
        return frame 
         
    def _refresh(self, nx, ny, temperature, field, phi, omega, gamma, 
                      time_count, monitor, int_u, int_d):
        if not(isinstance(nx, type(None))):
            self._p_nx = int(nx)
        if not(isinstance(ny, type(None))):
            self._p_ny = int(ny)
        if not(isinstance(temperature, type(None))):
            self._p_temperature = float(temperature)
        if not(isinstance(field, type(None))):
            self._p_field = float(field)
        if not(isinstance(phi, type(None))):
            self._p_phi = float(phi)
        if not(isinstance(omega, type(None))):
            self._p_omega = float(omega)
        if not(isinstance(gamma, type(None))):
            self._p_gamma = float(gamma)
        if not(isinstance(time_count, type(None))):
            self._p_time_count = float(time_count)
        if not(isinstance(monitor, type(None))):
            self._p_monitor = float(monitor)
        if not(isinstance(int_u, type(None))):
            self._p_int_u = int_u
        if not(isinstance(int_d, type(None))):
            self._p_int_d = int_d
     
    def set_val(self, nx=None, ny=None, temperature=None, field=None, phi=None, 
                omega=None, gamma=None, time_count=None, monitor=None, 
                int_u=None, int_d=None):
        self._refresh(nx, ny, temperature, field, phi, omega, gamma, 
                      time_count, monitor, int_u, int_d)
     
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
 nx, ny is the size of frame
 temperature is a temperature on the sample
 field is a magnetic field on the sample
 phi,  omega,  gamma are angles to describe orientaion of the sample and detector
 time_count is time of measurements
 monitor is a number of neutrons on the sample
        """
        print(lsout)
     
    def read_str_6T2_old(self, latr, lcont):
        """
        read information about frame from string, format 6T2_old
        """
        stype = "6T2"
        order = "C"
        if stype == "5C1":
            order = "C"
        if stype == "6T2":
            order = "F"    
            
        lhelp_1 = latr.split()
        s_phi, s_omega, s_gamma, s_temperature = None, None, None, None
        s_field, s_time, s_monitor = None, None, None

        l_phi, l_omega, l_gamma = "Phi=", "Omega=", "DeuxTheta="  
        l_temperature, l_field = "Temperature=","MagneticField="
        l_time, l_monitor = "ResultatTimeTotal=", "ResultatMoniteurUp=" 
        
        for hh in lhelp_1:
            if hh.startswith(l_phi):
                s_phi = hh[(len(l_phi)+1):-1]
                s_phi = 0.01*float(s_phi)
            elif hh.startswith(l_omega):
                s_omega = hh[(len(l_omega)+1):-1]
                s_omega = 0.01*float(s_omega)
            elif hh.startswith(l_gamma):
                s_gamma = hh[(len(l_gamma)+1):-1]
                s_gamma = 0.01*float(s_gamma)
            elif hh.startswith(l_temperature):
                s_temperature = hh[(len(l_temperature)+1):-1]
            elif hh.startswith(l_field):
                s_field = hh[(len(l_field)+1):-1]
            elif hh.startswith(l_time):
                s_time = hh[(len(l_time)+1):-1]
            elif hh.startswith(l_monitor):
                s_monitor = hh[(len(l_monitor)+1):-1]
                


        lab ="Data" 
        s_unt_u, s_unt_d = None, None
        h_str_atr, h_str_val = in_tag(lcont, lab)
        if "Data0" in h_str_atr:
            s_unt_u = h_str_val 
        elif "Data1" in h_str_atr:
            s_unt_d = h_str_val 
            
        lhelp_1 = h_str_atr.split()
        s_nx, s_ny = None, None 
        for hh in lhelp_1:
            if hh.startswith("x="):
                #it looks like there is mistake in definition of x and y in xml file
                s_ny = hh[3:-1]
            if hh.startswith("y="):
                #it looks like there is mistake in definition of x and y in xml file
                s_nx = hh[3:-1]
 

        val = lcont.find(lab+">")
        numb = val+len(lab)+1
       
        lab ="Data" 
        h_str_atr, h_str_val = in_tag(lcont[numb:], lab)
        if "Data0" in h_str_atr:
            s_unt_u = h_str_val 
        elif "Data1" in h_str_atr:
            s_unt_d = h_str_val 
        
        int_u, int_d = None, None
        if s_unt_u is not None:
            np_int_u_1d = numpy.array([int(hh) for hh in s_unt_u.split(";")], 
                                       dtype=int)
            np_int_d_1d = numpy.array([int(hh) for hh in s_unt_d.split(";")], 
                                       dtype=int)
            int_u = numpy.reshape(np_int_u_1d, (int(s_nx), int(s_ny)), 
                                  order=order)
            int_d = numpy.reshape(np_int_d_1d, (int(s_nx), int(s_ny)), 
                                  order=order)
            if order=="F":
                h1 = numpy.flip(int_u, axis=0)
                int_u = numpy.flip(h1, axis=1)
                
                h1 = numpy.flip(int_d, axis=0)
                int_d = numpy.flip(h1, axis=1)
        
        self.set_val(nx=s_nx, ny=s_ny, temperature=s_temperature, field=s_field, 
                     phi=s_phi, 
                omega=s_omega, gamma=s_gamma, time_count=s_time, monitor=s_monitor, 
                int_u=int_u, int_d=int_d)
        
     
    def read_str_5C1_6T2_new(self, latr, lcont, stype):
        """
        read information about frame from string, format 5C1_new
        """
        order = "C"
        if stype == "VIP":
            order = "C"
        if stype == "6T2":
            order = "C"
        
        lab ="Phi" 
        h_str_atr, s_phi = in_tag(lcont , lab)

        lab ="Omega" 
        h_str_atr, s_omega = in_tag(lcont , lab)

        lab ="Gamma" 
        h_str_atr, s_gamma = in_tag(lcont , lab)

        lab ="temperature" 
        h_str_atr, s_temperature = in_tag(lcont , lab)
        
        lab ="magneticField" 
        h_str_atr, s_field = in_tag(lcont , lab)

        lab ="totaltimecount" 
        h_str_atr, s_time = in_tag(lcont , lab)

        lab ="totalmonitorcount" 
        h_str_atr, s_monitor = in_tag(lcont , lab)

        lab ="Data" 
        s_unt_u, s_unt_d = None, None
        h_str_atr, h_str_val = in_tag(lcont, lab)
        if "IntensityUp" in h_str_atr:
            s_unt_u = h_str_val 
        elif "IntensityDown" in h_str_atr:
            s_unt_d = h_str_val 
            
        lhelp_1 = h_str_atr.split()
        s_nx, s_ny = None, None 
        for hh in lhelp_1:
            if hh.startswith("x="):
                #it looks like there is mistake in definition of x and y in xml file
                s_ny = hh[3:-1]
            elif hh.startswith("y="):
                #it looks like there is mistake in definition of x and y in xml file
                s_nx = hh[3:-1]
 

        val = lcont.find(lab+">")
        numb = val+len(lab)+1
       
        lab ="Data" 
        h_str_atr, h_str_val = in_tag(lcont[numb:], lab)
        if "IntensityUp" in h_str_atr:
            s_unt_u = h_str_val 
        elif "IntensityDown" in h_str_atr:
            s_unt_d = h_str_val 
        
        int_u, int_d = None, None
        if s_unt_u is not None:
            np_int_u_1d = numpy.array([int(hh) for hh in s_unt_u.split(";")], 
                                       dtype=int)
            np_int_d_1d = numpy.array([int(hh) for hh in s_unt_d.split(";")], 
                                       dtype=int)
            
            int_u = numpy.reshape(np_int_u_1d, (int(s_nx), int(s_ny)), 
                                  order=order)
            int_d = numpy.reshape(np_int_d_1d, (int(s_nx), int(s_ny)), 
                                  order=order)
            if stype=="6T2":
                h1 = numpy.flip(int_u, axis=0)
                int_u = numpy.flip(h1, axis=1)
                
                h1 = numpy.flip(int_d, axis=0)
                int_d = numpy.flip(h1, axis=1)
        
        self.set_val(nx=s_nx, ny=s_ny, temperature=s_temperature, field=s_field, 
                     phi=s_phi, 
                omega=s_omega, gamma=s_gamma, time_count=s_time, monitor=s_monitor, 
                int_u=int_u, int_d=int_d)
     
    def plot_int_u(self):
        int_u = self._p_int_u
        if int_u is not None:
            matplotlib.pyplot.imshow(int_u.transpose(), aspect="auto", origin="lower")
     
    def plot_int_d(self):
        int_d = self._p_int_d
        if int_d is not None:
            matplotlib.pyplot.imshow(int_d.transpose(), aspect="auto", origin="lower")
     
    def plot_int_sum(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u+int_d).transpose(), aspect="auto", origin="lower")
     
    def plot_int_dif(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u-int_d).transpose(), aspect="auto", origin="lower")


class Setup(object):
    """
    Class to describe information from position sensitive detector
    """
    def __init__(self, name="unknown"):
        super(Setup, self).__init__()
        self._p_name = None
        self._p_mask = None 
        self._p_nx_psd = None 
        self._p_nx_min = None 
        self._p_nx_max = None 
        self._p_ny_psd = None 
        self._p_ny_min = None 
        self._p_ny_max = None 
        self._p_psd_gamma_0 = None
        self._p_x_step_r = None
        self._p_psd_equat = None
        self._p_psd_height = None
        self._p_psd_width = None
        self._p_dist_sd = None
        
        self._refresh(name)
         
    def __repr__(self):
        lsout = """Setup: \n name: {:} \n nx: {:3}   min: {:3}   max: {:3}
 ny: {:3}   min: {:3}   max: {:3} \n gamma_0: {:} rad.,  psd_equat: {:}
 psd_width: {:} cm.,  psd_height: {:} cm.\n distance sample-detector: {:} cm.
""".format(self._p_name, self._p_nx_psd, self._p_nx_min, self._p_nx_max,
           self._p_ny_psd, self._p_ny_min, self._p_ny_max,
           self._p_psd_gamma_0, self._p_psd_equat,
           self._p_psd_width, self._p_psd_height, self._p_dist_sd)
        return lsout
     
    def _recalc(self):
        if self._p_nx_psd is None:
            return
        val_x = numpy.linspace(0, self._p_nx_psd-1, num=self._p_nx_psd, dtype=int)
        mask_x = numpy.array(((val_x>=self._p_nx_min)&(val_x<self._p_nx_max)), 
                             dtype=bool)
        
        val_y = numpy.linspace(0, self._p_ny_psd-1, num=self._p_ny_psd, dtype=int)
        mask_y = numpy.array(((val_y>=self._p_ny_min)&(val_y<self._p_ny_max)), 
                             dtype=bool)
        
        mask_x_2d, mask_y_2d = numpy.meshgrid(mask_x, mask_y, indexing="ij")
        
        np_x, np_y = numpy.meshgrid(val_x, val_y, indexing="ij")
        self._p_np_x = np_x
        self._p_np_y = np_y
        self._p_mask = numpy.logical_and(mask_x_2d, mask_y_2d)
        
        
        if ((self._p_name == "5C1")|(self._p_name == "VIP")):
            np_g_0_5C1 = self._p_x_step_r * self._p_np_x - self._p_psd_gamma_0
            np_n_0_5C16T2 = numpy.arctan2((np_y-self._p_psd_equat)*
                      self._p_psd_height*(1./self._p_ny_psd), self._p_dist_sd)
            self._p_np_g_0_5C1 = np_g_0_5C1
            self._p_np_n_0_5C16T2 = np_n_0_5C16T2
            
        elif self._p_name == "6T2":
            np_g_0_6T2 =numpy.arctan2((np_x-0.5*(self._p_nx_psd-1))*self._p_psd_width*(
                1./self._p_nx_psd), self._p_dist_sd)
            np_n_0_5C16T2 = numpy.arctan2((np_y-self._p_psd_equat)*
                      self._p_psd_height*(1./self._p_ny_psd), self._p_dist_sd)
            self._p_np_g_0_6T2 = np_g_0_6T2
            self._p_np_n_0_5C16T2 = np_n_0_5C16T2

     
    def _refresh(self, name):
        if not(isinstance(name, type(None))):
            self._p_name = name
            if ((name == "5C1")|(name == "VIP")):
                self._set_5C1()
            elif name == "6T2":
                self._set_6T2()
        self._recalc()
     
    def set_val(self, name=None):
        self._refresh(name)
     
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
 name is the name of diffractometer '5C1' or 6T2
        """
        print(lsout)
     
    def _set_5C1(self):
        """
        set default parameters for diffractometer 5C1
        """
        self._p_nx_psd = 64
        self._p_ny_psd = 256

        self._p_dist_sd = 98.
        self._p_psd_height = 0.703*256./4.
        self._p_psd_equat = 56.
        self._p_x_step_r = 0.0225

        self._p_psd_gamma_0 = 0.

        self._p_nx_min = 2
        self._p_nx_max = 62

        self._p_ny_min = 4
        self._p_ny_max = 252
     
    def _set_6T2(self):
        """
        set default parameters for diffractometer 5C1
        """
        self._p_nx_psd = 128
        self._p_ny_psd = 128
        
        self._p_dist_sd = 57.
        self._p_psd_height = 25.6
        self._p_psd_width = 25.6
        self._p_psd_equat = 56.
        self._p_psd_gamma_0 = 0.

        self._p_nx_min = 0
        self._p_nx_max = 128

        self._p_ny_min = 0
        self._p_ny_max = 128
     
    def calc_gn(self, gamma):
        """
        transform coordinates of detector to gamma, nu
        """
            
        if ((self._p_name == "5C1")|(self._p_name == "VIP")):
            np_g, np_n = self._calc_gn_5C1(gamma)
        elif self._p_name == "6T2":
            np_g, np_n = self._calc_gn_6T2(gamma)
        return np_g, np_n
     
    def _calc_gn_5C1(self, gamma):
        """
        transform coordinates of detector to gamma, nu for 5C1 definition
        """
        
        np_g = self._p_np_g_0_5C1 + gamma
        np_n = self._p_np_n_0_5C16T2
        return np_g, np_n
     
    def _calc_gn_6T2(self, gamma):
        """
        transform coordinates of detector to gamma, nu for 6T2 definition
        """
        np_g = self._p_np_g_0_6T2 + gamma

        np_n = self._p_np_n_0_5C16T2
        return np_g, np_n
     
    def calc_tpgn(self, gamma, mode=True):
        """
        ttheta is diffraction angle
        phi is polar angle
        """
        np_g, np_n = self.calc_gn(gamma)
        
        np_tth = numpy.arccos(numpy.cos(np_g)*numpy.cos(np_n))
        np_phi = numpy.arctan2(numpy.tan(np_n), numpy.sin(np_g))
        return np_tth, np_phi, np_g, np_n 



class GridGammaNu(object):
    """
    Class to recalculate data in gamma nu grid
    """
    def __init__(self, g_min=4.*numpy.pi/180., g_max=80.*numpy.pi/180., 
                 g_step=0.1*numpy.pi/180., 
                 n_min=-20.*numpy.pi/180., n_max=20.*numpy.pi/180., 
                 n_step=0.1*numpy.pi/180.):
        
        super(GridGammaNu, self).__init__()
        self._p_g_min = None
        self._p_g_max = None
        self._p_g_step = None
        self._p_n_min = None
        self._p_n_max = None
        self._p_n_step = None
        self._p_grid_g = None
        self._p_grid_n = None
        self._p_n_g = None
        self._p_n_n = None
        
        self._p_int_u = None
        self._p_sint_u = None
        self._p_int_d = None
        self._p_sint_d = None
        
        self._refresh(g_min, g_max, g_step, n_min, n_max, n_step)
        
    def __repr__(self):
        lsout = """Frame: \n g_min: {:}, g_max: {:}, g_step: {:}
 n_min: {:}, n_max: {:}, n_step: {:}""".format(
 self._p_g_min*180./numpy.pi, self._p_g_max*180./numpy.pi, 
 self._p_g_step*180./numpy.pi, 
 self._p_n_min*180./numpy.pi, self._p_n_max*180./numpy.pi, 
 self._p_n_step*180./numpy.pi)
        return lsout
     
    def _calc_grid(self):
        g_min, g_max = self._p_g_min, self._p_g_max
        g_step = self._p_g_step
        cond_g = ((g_min is not None)&(g_max is not None)&(g_step is not None))

        n_min, n_max = self._p_n_min, self._p_n_max
        n_step = self._p_n_step
        cond_n = ((n_min is not None)&(n_max is not None)&(n_step is not None))
        
        if (cond_g & cond_n):
            n_g = int(round((g_max-g_min)*1./g_step)+1)
            val_g = numpy.linspace(g_min, g_max, n_g)
            self._p_g_step = (g_max-g_min)/float(n_g-1)
            self._p_n_g = n_g

            n_n = int(round((n_max-n_min)*1./n_step)+1)
            val_n = numpy.linspace(n_min, n_max, n_n)
            grid_g, grid_n = numpy.meshgrid(val_g, val_n, indexing="ij")
            self._p_n_step = (n_max-n_min)/float(n_n-1)
            self._p_n_n = n_n
        
            self._p_grid_g = grid_g
            self._p_grid_n = grid_n
            
            self._p_int_u  = numpy.nan * grid_g
            self._p_sint_u  = numpy.nan * grid_g
            self._p_int_d  = numpy.nan * grid_g
            self._p_sint_d  = numpy.nan * grid_g

            
     
    def _refresh(self, g_min, g_max, g_step, n_min, n_max, n_step):
        if not(isinstance(g_min, type(None))):
            self._p_g_min = g_min
        if not(isinstance(g_max, type(None))):
            self._p_g_max = g_max
        if not(isinstance(g_step, type(None))):
            self._p_g_step = g_step
        if not(isinstance(n_min, type(None))):
            self._p_n_min = n_min
        if not(isinstance(n_max, type(None))):
            self._p_n_max = n_max
        if not(isinstance(n_step, type(None))):
            self._p_n_step = n_step        
            
        self._calc_grid()
     
    def set_val(self, g_min=None, g_max=None, g_step=None, n_min=None, 
                n_max=None, n_step=None):
        self._refresh(g_min, g_max, g_step, n_min, n_max, n_step)
     
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
g_min, g_max, g_step are minimal, maximal and step of gamma (in radians)
n_min, n_max, n_step are minimal, maximal and step of nu (in radians)       
"""
        print(lsout)
        
     
    def load_to_grid(self, np_g, np_n, i_u, si_u, i_d, si_d):

        g_min = self._p_g_min
        g_step = self._p_g_step
        
        n_min = self._p_n_min
        n_step = self._p_n_step

        int_u  = numpy.nan * self._p_int_u
        sint_u  = numpy.nan * self._p_sint_u
        int_d  = numpy.nan * self._p_int_d
        sint_d  = numpy.nan * self._p_sint_d
        
        n_g = self._p_n_g 
        n_n = self._p_n_n 
        
        
        ind_g = (numpy.rint((np_g-g_min)*1./g_step)).astype(numpy.int)
        ind_n = (numpy.rint((np_n-n_min)*1./n_step)).astype(numpy.int)
        
        flag_g = numpy.logical_and(ind_g >= 0, ind_g < n_g)
        flag_n = numpy.logical_and(ind_n >= 0, ind_n < n_n)
        flag = numpy.logical_and(flag_g, flag_n)
        
        ind_g_1 = ind_g[flag]
        ind_n_1 = ind_n[flag]

        i_u_1 = i_u[flag]
        si_u_1 = si_u[flag] 
        i_d_1 = i_d[flag]
        si_d_1 = si_d[flag]
        
        ing_gn_1 = numpy.vstack([ind_g_1, ind_n_1]).transpose()
        ing_gn_uniq = numpy.unique(ing_gn_1, axis=0)
        ing_g_uniq = ing_gn_uniq[:,0]
        ing_n_uniq = ing_gn_uniq[:,1]
        flag_perc = False
        nmax = ing_gn_uniq.shape[0]
        #flag_all = numpy.ones(ind_g_1.size).astype(numpy.int)
        ind = 0
        for ind_g, ind_n in zip(ing_g_uniq, ing_n_uniq):
            flag = (ing_gn_1 == (ind_g, ind_n)).all(axis=1)
            
            if flag.any():
                flag_u_2 = numpy.nonzero(si_u_1[flag])
                i_u_av = i_u_1[flag][flag_u_2]
                si_u_av_isq = 1./(si_u_1[flag][flag_u_2])**2
                
                flag_d_2 = numpy.nonzero(si_d_1[flag])
                i_d_av = i_d_1[flag][flag_d_2]
                si_d_av_isq = 1./(si_d_1[flag][flag_d_2])**2
            
                si_u_val_isq = (si_u_av_isq).sum()
                i_u_val = (1./si_u_val_isq)*(i_u_av*si_u_av_isq).sum()
                si_u_val = 1./(si_u_val_isq**0.5)
        
                si_d_val_isq = (si_d_av_isq).sum()
                i_d_val = (1./si_d_val_isq)*(i_d_av*si_d_av_isq).sum()
                si_d_val = 1./(si_d_val_isq**0.5)
            
                int_u[ind_g, ind_n] = i_u_val
                sint_u[ind_g, ind_n] = si_u_val
                int_d[ind_g, ind_n] = i_d_val
                sint_d[ind_g, ind_n] = si_d_val 
                
                """
                flag_all[flag_all] = numpy.logical_not(flag)
                
                ing_gn_1  = numpy.delete(ing_gn_1, flag, axis=0)
                i_u_1 = numpy.delete(i_u_1, flag)
                si_u_1 = numpy.delete(si_u_1, flag)
                i_d_1 = numpy.delete(i_d_1, flag)
                si_d_1 = numpy.delete(si_d_1, flag)
                """
            perc = int(round(100.*float(ind)/float(nmax)))
            if ((perc%2 == 1)&(flag_perc)):
                flag_perc = False
                print(" {:3} %".format(perc))
            elif ((perc%2 == 0)&(not(flag_perc))):
                flag_perc = True
                print(" {:3} %".format(perc))
            ind += 1
            
        self._p_int_u = int_u
        self._p_sint_u = sint_u
        self._p_int_d = int_d
        self._p_sint_d = sint_d
      
    def plot_int_u(self):
        int_u = self._p_int_u
        if int_u is not None:
            matplotlib.pyplot.imshow(int_u.transpose(), aspect="auto", origin="lower")
     
    def plot_int_d(self):
        int_d = self._p_int_d
        if int_d is not None:
            matplotlib.pyplot.imshow(int_d.transpose(), aspect="auto", origin="lower")
     
    def plot_int_sum(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u+int_d).transpose(), aspect="auto", origin="lower")
     
    def plot_int_dif(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u-int_d).transpose(), aspect="auto", origin="lower")       
          
    def save_u_d_to_file(self, f_name):
        """
        save up and down in file
        """
        int_u = self._p_int_u
        sint_u = self._p_sint_u
        int_d = self._p_int_d
        sint_d = self._p_sint_d
        
        grid_g = self._p_grid_g
        grid_n = self._p_grid_n
        
        ax_1 = grid_g[:,0]*180./numpy.pi
        ax_2 = grid_n[0,:]*180./numpy.pi
        lstr_i_u = save_matrix_to_str(int_u, ax_1, ax_2)
        lstr_si_u = save_matrix_to_str(sint_u, ax_1, ax_2)
        lstr_i_d = save_matrix_to_str(int_d, ax_1, ax_2)
        lstr_si_d = save_matrix_to_str(sint_d, ax_1, ax_2)
        
        lstr = ("\n".join(lstr_i_u) + 3*"\n"+
                "\n".join(lstr_si_u) + 3*"\n"+
                "\n".join(lstr_i_d) + 3*"\n"+
                "\n".join(lstr_si_d) + 3*"\n")
        
        fid = open(f_name, "w")
        fid.write(lstr)
        fid.close()
        
    
         
    def save_sum_diff_to_file(self, f_name):
        """
        save sum and difference in file
        """
        int_u = self._p_int_u
        sint_u = self._p_sint_u
        int_d = self._p_int_d
        sint_d = self._p_sint_d
        
        grid_g = self._p_grid_g
        grid_n = self._p_grid_n
        
        ax_1 = grid_g[:,0]*180./numpy.pi
        ax_2 = grid_n[0,:]*180./numpy.pi
        int_sum = int_u+int_d
        sint_sumu = (sint_u**2+sint_d**2)**0.5
        int_dif = int_u-int_d
        sint_dif = sint_sumu
        lstr_i_sum = save_matrix_to_str(int_sum, ax_1, ax_2)
        lstr_si_sum = save_matrix_to_str(sint_sumu, ax_1, ax_2)
        lstr_i_dif = save_matrix_to_str(int_dif, ax_1, ax_2)
        lstr_si_dif = save_matrix_to_str(sint_dif, ax_1, ax_2)
        
        lstr = ("\n".join(lstr_i_sum) + 3*"\n"+
                "\n".join(lstr_si_sum) + 3*"\n"+
                "\n".join(lstr_i_dif) + 3*"\n"+
                "\n".join(lstr_si_dif) + 3*"\n")
        
        fid = open(f_name, "w")
        fid.write(lstr)
        fid.close()
        

class GridTthetaPhi(object):
    """
    Class to recalculate data in gamma nu grid
    """
    def __init__(self, tth_min=4.*numpy.pi/180., tth_max=80.*numpy.pi/180., 
                 tth_step=0.1*numpy.pi/180., 
                 phi_min=-10.*numpy.pi/180., phi_max=40.*numpy.pi/180., 
                 phi_step=0.5*numpy.pi/180.):
        
        super(GridTthetaPhi, self).__init__()
        self._p_tth_min = None
        self._p_tth_max = None
        self._p_tth_step = None
        self._p_phi_min = None
        self._p_phi_max = None
        self._p_phi_step = None
        self._p_grid_tth = None
        self._p_grid_phi = None
        self._p_n_tth = None
        self._p_n_phi = None
        
        self._p_int_u = None
        self._p_sint_u = None
        self._p_int_d = None
        self._p_sint_d = None
        
        self._p_1d_tth = None
        self._p_1d_int_u = None
        self._p_1d_sint_u = None
        self._p_1d_int_d = None
        self._p_1d_sint_d  = None 
        
        self._refresh(tth_min, tth_max, tth_step, phi_min, phi_max, phi_step)
        
    def __repr__(self):
        lsout = """Frame: \n tth_min: {:}, tth_max: {:}, tth_step: {:}
 phi_min: {:}, phi_max: {:}, phi_step: {:}""".format(
 self._p_tth_min*180./numpy.pi, self._p_tth_max*180./numpy.pi, 
 self._p_tth_step*180./numpy.pi, 
 self._p_phi_min*180./numpy.pi, self._p_phi_max*180./numpy.pi, 
 self._p_phi_step*180./numpy.pi)
        return lsout
     
    def _calc_grid(self):
        tth_min, tth_max = self._p_tth_min, self._p_tth_max
        tth_step = self._p_tth_step
        cond_tth = ((tth_min is not None)&(tth_max is not None)&(tth_step is not None))

        phi_min, phi_max = self._p_phi_min, self._p_phi_max
        phi_step = self._p_phi_step
        cond_phi = ((phi_min is not None)&(phi_max is not None)&(phi_step is not None))
        
        if (cond_tth & cond_phi):
            n_tth = int(round((tth_max-tth_min)*1./tth_step)+1)
            val_tth = numpy.linspace(tth_min, tth_max, n_tth)
            self._p_tth_step = (tth_max-tth_min)/float(n_tth-1)
            self._p_n_tth = n_tth

            n_phi = int(round((phi_max-phi_min)*1./phi_step)+1)
            val_phi = numpy.linspace(phi_min, phi_max, n_phi)
            grid_tth, grid_phi = numpy.meshgrid(val_tth, val_phi, indexing="ij")
            self._p_phi_step = (phi_max-phi_min)/float(n_phi-1)
            self._p_n_phi = n_phi
        
            self._p_grid_tth = grid_tth
            self._p_grid_phi = grid_phi
            
            self._p_int_u  = numpy.nan * grid_tth
            self._p_sint_u  = numpy.nan * grid_tth
            self._p_int_d  = numpy.nan * grid_tth
            self._p_sint_d  = numpy.nan * grid_tth

            
     
    def _refresh(self, tth_min, tth_max, tth_step, phi_min, phi_max, phi_step):
        if not(isinstance(tth_min, type(None))):
            self._p_tth_min = tth_min
        if not(isinstance(tth_max, type(None))):
            self._p_tth_max = tth_max
        if not(isinstance(tth_step, type(None))):
            self._p_tth_step = tth_step
        if not(isinstance(phi_min, type(None))):
            self._p_phi_min = phi_min
        if not(isinstance(phi_max, type(None))):
            self._p_phi_max = phi_max
        if not(isinstance(phi_step, type(None))):
            self._p_phi_step = phi_step        
            
        self._calc_grid()
     
    def set_val(self, tth_min=None, tth_max=None, tth_step=None, phi_min=None, 
                phi_max=None, phi_step=None):
        self._refresh(tth_min, tth_max, tth_step, phi_min, phi_max, phi_step)
     
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
tth_min, tth_max, tth_step are minimal, maximal and step of gamma (in radians)
phi_min, phi_max, phi_step are minimal, maximal and step of nu (in radians)       
"""
        print(lsout)
    
     
    def load_to_grid(self, np_tth, np_phi, i_u, si_u, i_d, si_d):

        tth_min = self._p_tth_min
        tth_step = self._p_tth_step
        
        phi_min = self._p_phi_min
        phi_step = self._p_phi_step

        int_u  = numpy.nan * self._p_int_u
        sint_u  = numpy.nan * self._p_sint_u
        int_d  = numpy.nan * self._p_int_d
        sint_d  = numpy.nan * self._p_sint_d
        
        n_tth = self._p_n_tth
        n_phi = self._p_n_phi
        
        
        ind_tth = (numpy.rint((np_tth-tth_min)*1./tth_step)).astype(numpy.int)
        ind_phi = (numpy.rint((np_phi-phi_min)*1./phi_step)).astype(numpy.int)
        
        flag_tth = numpy.logical_and(ind_tth >= 0, ind_tth < n_tth)
        flag_phi = numpy.logical_and(ind_phi >= 0, ind_phi < n_phi)
        flag = numpy.logical_and(flag_tth, flag_phi)
        
        ind_tth_1 = ind_tth[flag]
        ind_phi_1 = ind_phi[flag]

        i_u_1 = i_u[flag]
        si_u_1 = si_u[flag] 
        i_d_1 = i_d[flag]
        si_d_1 = si_d[flag]
        
        ing_tp_1 = numpy.vstack([ind_tth_1, ind_phi_1]).transpose()
        ing_tp_uniq = numpy.unique(ing_tp_1, axis=0)
        flag_perc = False
        nmax = ing_tp_uniq.shape[0]
        #flag_all = numpy.ones(ind_g_1.size).astype(numpy.int)
        for ind in range(ing_tp_uniq.shape[0]):
            ind_tth = ing_tp_uniq[ind, 0]
            ind_phi = ing_tp_uniq[ind, 1]
            
            flag = (ing_tp_1 == (ind_tth, ind_phi)).all(axis=1)
            
            if flag.any():
                flag_u_2 = numpy.nonzero(si_u_1[flag])
                i_u_av = i_u_1[flag][flag_u_2]
                si_u_av_isq = 1./(si_u_1[flag][flag_u_2])**2
                
                flag_d_2 = numpy.nonzero(si_d_1[flag])
                i_d_av = i_d_1[flag][flag_d_2]
                si_d_av_isq = 1./(si_d_1[flag][flag_d_2])**2
            
                si_u_val_isq = (si_u_av_isq).sum()
                i_u_val = (1./si_u_val_isq)*(i_u_av*si_u_av_isq).sum()
                si_u_val = 1./(si_u_val_isq**0.5)
        
                si_d_val_isq = (si_d_av_isq).sum()
                i_d_val = (1./si_d_val_isq)*(i_d_av*si_d_av_isq).sum()
                si_d_val = 1./(si_d_val_isq**0.5)
            
                int_u[ind_tth, ind_phi] = i_u_val
                sint_u[ind_tth, ind_phi] = si_u_val
                int_d[ind_tth, ind_phi] = i_d_val
                sint_d[ind_tth, ind_phi] = si_d_val 
                
                """
                flag_all[flag_all] = numpy.logical_not(flag)
                
                ing_gn_1  = numpy.delete(ing_gn_1, flag, axis=0)
                i_u_1 = numpy.delete(i_u_1, flag)
                si_u_1 = numpy.delete(si_u_1, flag)
                i_d_1 = numpy.delete(i_d_1, flag)
                si_d_1 = numpy.delete(si_d_1, flag)
                """
            
            perc = int(round(100.*float(ind)/float(nmax)))
            if ((perc%2 == 1)&(flag_perc)):
                flag_perc = False
                print(" {:3} %".format(perc))
            elif ((perc%2 == 0)&(not(flag_perc))):
                flag_perc = True
                print(" {:3} %".format(perc))
            
        self._p_int_u = int_u
        self._p_sint_u = sint_u
        self._p_int_d = int_d
        self._p_sint_d = sint_d
     
    def plot_int_u(self):
        int_u = self._p_int_u
        if int_u is not None:
            matplotlib.pyplot.imshow(int_u.transpose(), aspect="auto", origin="lower")
     
    def plot_int_d(self):
        int_d = self._p_int_d
        if int_d is not None:
            matplotlib.pyplot.imshow(int_d.transpose(), aspect="auto", origin="lower")
     
    def plot_int_sum(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u+int_d).transpose(), aspect="auto", origin="lower")
     
    def plot_int_dif(self):
        int_u = self._p_int_u
        int_d = self._p_int_d
        if (int_u is not None)&(int_d is not None):
            matplotlib.pyplot.imshow((int_u-int_d).transpose(), aspect="auto", origin="lower")       
            
    def save_u_d_to_file(self, f_name):
        """
        save up and down in file
        """
        int_u = self._p_int_u
        sint_u = self._p_sint_u
        int_d = self._p_int_d
        sint_d = self._p_sint_d
        
        grid_tth = self._p_grid_tth
        grid_phi = self._p_grid_phi
        
        ax_1 = grid_tth[:,0]*180./numpy.pi
        ax_2 = grid_phi[0,:]*180./numpy.pi
        lstr_i_u = save_matrix_to_str(int_u, ax_1, ax_2)
        lstr_si_u = save_matrix_to_str(sint_u, ax_1, ax_2)
        lstr_i_d = save_matrix_to_str(int_d, ax_1, ax_2)
        lstr_si_d = save_matrix_to_str(sint_d, ax_1, ax_2)
        
        lstr = ("\n".join(lstr_i_u) + 3*"\n"+
                "\n".join(lstr_si_u) + 3*"\n"+
                "\n".join(lstr_i_d) + 3*"\n"+
                "\n".join(lstr_si_d) + 3*"\n")
        
        fid = open(f_name, "w")
        fid.write(lstr)
        fid.close()
        
    
         
    def save_sum_diff_to_file(self, f_name):
        """
        save sum and difference in file
        """
        int_u = self._p_int_u
        sint_u = self._p_sint_u
        int_d = self._p_int_d
        sint_d = self._p_sint_d
        
        grid_tth = self._p_grid_tth
        grid_phi = self._p_grid_phi
        
        ax_1 = grid_tth[:,0]*180./numpy.pi
        ax_2 = grid_phi[0,:]*180./numpy.pi
        int_sum = int_u+int_d
        sint_sumu = (sint_u**2+sint_d**2)**0.5
        int_dif = int_u-int_d
        sint_dif = sint_sumu
        lstr_i_sum = save_matrix_to_str(int_sum, ax_1, ax_2)
        lstr_si_sum = save_matrix_to_str(sint_sumu, ax_1, ax_2)
        lstr_i_dif = save_matrix_to_str(int_dif, ax_1, ax_2)
        lstr_si_dif = save_matrix_to_str(sint_dif, ax_1, ax_2)
        
        lstr = ("\n".join(lstr_i_sum) + 3*"\n"+
                "\n".join(lstr_si_sum) + 3*"\n"+
                "\n".join(lstr_i_dif) + 3*"\n"+
                "\n".join(lstr_si_dif) + 3*"\n")
        
        fid = open(f_name, "w")
        fid.write(lstr)
        fid.close()
     
    def calc_1D_profile(self, phi_min=-2., phi_max=10.):
        """
        Calculate  1D diffraction profile
        """
        phi_min_r = phi_min*numpy.pi*1./180.
        phi_max_r = phi_max*numpy.pi*1./180.
        
        int_u = self._p_int_u
        sint_u = self._p_sint_u
        int_d = self._p_int_d
        sint_d = self._p_sint_d

        grid_tth = self._p_grid_tth
        grid_phi = self._p_grid_phi
        
        ltth, li_u, lsi_u, li_d, lsi_d = [], [], [], [], []
        for tth_1d, phi_1d, i_u_1d, si_u_1d, i_d_1d, si_d_1d in zip(
                grid_tth, grid_phi, int_u, sint_u, int_d, sint_d):
            flag = numpy.logical_and(phi_1d >= phi_min_r, phi_1d <= phi_max_r)
            tth = tth_1d[0]
            si_u_isd = 1./(si_u_1d[flag])**2
            si_d_isd = 1./(si_d_1d[flag])**2
            
            si_u_sd_val = 1./si_u_isd.sum()
            i_u_val = si_u_sd_val*(i_u_1d[flag]*si_u_isd).sum()
            si_u_val = si_u_sd_val**0.5
            
            si_d_sd_val = 1./si_d_isd.sum()
            i_d_val = si_d_sd_val*(i_d_1d[flag]*si_d_isd).sum()
            si_d_val = si_d_sd_val**0.5

            cond_u= not(numpy.isnan(i_u_val)|numpy.isnan(si_u_val))
            cond_d= not(numpy.isnan(i_d_val)|numpy.isnan(si_d_val))
            
            if (cond_u&cond_d):
                ltth.append(float(tth))
                li_u.append(float(i_u_val))
                lsi_u.append(float(si_u_val)) 
                li_d.append(float(i_d_val))
                lsi_d.append(float(si_d_val)) 
            
        self._p_1d_tth = numpy.array(ltth, dtype=float)
        self._p_1d_int_u = numpy.array(li_u, dtype=float)
        self._p_1d_sint_u = numpy.array(lsi_u, dtype=float)
        self._p_1d_int_d = numpy.array(li_d, dtype=float)
        self._p_1d_sint_d  = numpy.array(lsi_d, dtype=float)
        
    def plot_1d_u_d(self):
        """
        plot up and down
        """
        col_1 = "#FF0000"
        col_2 = "#0000FF"
        tth = self._p_1d_tth
        int_u = self._p_1d_int_u 
        sint_u = self._p_1d_sint_u 
        int_d = self._p_1d_int_d 
        sint_d = self._p_1d_sint_d 
        
        if int_u is not None:
            matplotlib.pyplot.plot(tth*180./numpy.pi, int_u, "k-", 
                                   linewidth=1.0)
            
            matplotlib.pyplot.errorbar(tth*180./numpy.pi, int_u, yerr = sint_u, 
                                       ecolor = col_1,  fmt='.', color=col_1, 
                                       linewidth = 0.5)

        if int_d is not None:
            matplotlib.pyplot.plot(tth*180./numpy.pi, int_d, "k-", 
                                   linewidth=1.0)
            
            matplotlib.pyplot.errorbar(tth*180./numpy.pi, int_d, yerr = sint_d, 
                                       ecolor = col_2,  fmt='.', color=col_2, 
                                       linewidth = 0.5)
     
    def plot_1d_sum_diff(self):
        """
        plot sum and difference
        """
        col_1 = "#FF0000"
        col_2 = "#0000FF"
        tth = self._p_1d_tth
        int_sum = self._p_1d_int_u + self._p_1d_int_d 
        sint_sum = (self._p_1d_sint_u**2 + self._p_1d_sint_d**2)**0.5
        int_dif = self._p_1d_int_u - self._p_1d_int_d 
        
        if int_sum is not None:
            matplotlib.pyplot.plot(tth*180./numpy.pi, int_sum, "k-", 
                                   linewidth=1.0)

            matplotlib.pyplot.errorbar(tth*180./numpy.pi, int_sum, yerr = sint_sum, 
                                       ecolor = col_1,  fmt='.', color=col_1, 
                                       linewidth = 0.5)

        if int_dif is not None:
            matplotlib.pyplot.plot(tth*180./numpy.pi, int_dif, "k-", 
                                   linewidth=1.0)
            matplotlib.pyplot.errorbar(tth*180./numpy.pi, int_dif, yerr = sint_sum, 
                                       ecolor = col_2,  fmt='.', color=col_2, 
                                       linewidth = 0.5)

     
    def save_1D_u_d(self, f_name):
        """
        save 1d profile for up and down
        """
        tth = self._p_1d_tth 
        int_u = self._p_1d_int_u 
        sint_u = self._p_1d_sint_u
        int_d = self._p_1d_int_d 
        sint_d = self._p_1d_sint_d
        if ((tth is None)|(int_u is None)|(sint_u is None)|(int_d is None)|
                (sint_d is None)):
            self.calc_1D_profile()
            tth = self._p_1d_tth 
            int_u = self._p_1d_int_u 
            sint_u = self._p_1d_sint_u
            int_d = self._p_1d_int_d 
            sint_d = self._p_1d_sint_d
        lsout = ["#    ttheta       IntUP      sIntUP     IntDOWN    sIntDOWN"]
        for h1, h2, h3, h4, h5 in zip(tth, int_u, sint_u, int_d, sint_d) :
            lsout.append(" {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}".format(
                    h1* 180./numpy.pi, h2, h3, h4, h5))
            
        fid = open(f_name, "w")
        fid.write("\n".join(lsout))
        fid.close()
     
    def save_1D_sum_diff(self, f_name):
        """
        save 1d profile for sum and difference
        """
        tth = self._p_1d_tth * 180./numpy.pi
        int_u = self._p_1d_int_u 
        sint_u = self._p_1d_sint_u
        int_d = self._p_1d_int_d 
        sint_d = self._p_1d_sint_d
        if ((tth is None)|(int_u is None)|(sint_u is None)|(int_d is None)|
                (sint_d is None)):
            self.calc_1D_profile()
            tth = self._p_1d_tth * 180./numpy.pi
            int_u = self._p_1d_int_u 
            sint_u = self._p_1d_sint_u
            int_d = self._p_1d_int_d 
            sint_d = self._p_1d_sint_d

        int_sum = int_u+int_d
        sint_sum = (sint_u**2 + sint_d**2)**0.5
        int_dif = int_u-int_d
        sint_dif = sint_sum
        lsout = ["#    ttheta      IntSUM     sIntSUM     IntDIFF    sIntDIFF"]
        for h1, h2, h3, h4, h5 in zip(tth, int_sum, sint_sum, int_dif, 
                                      sint_dif) :
            lsout.append(" {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}  {:10.3f}".format(
                    h1, h2, h3, h4, h5))
            
        fid = open(f_name, "w")
        fid.write("\n".join(lsout))
        fid.close()
 
def to_do_job(xml_1, ljob, f_name):
    """
    Description of typical scenarios
    """
    f_base = ".".join(f_name.split(".")[:-1])
    if (("1" in ljob)|("2" in ljob)):
        print("\nConverting data to gamma-nu coordinates")
        grid_gamma_nu = xml_1.calc_in_grid_gamma_nu()
        
    if ("1" in ljob):
        f_out = f_base+"_u_d_gn.dat"
        grid_gamma_nu.save_u_d_to_file(f_out)
        print("2D profile for 'up' and 'down' in gamma-nu coordinates is saved in the file '{:}'".format(f_out))
    if ("2" in ljob):
        f_out = f_base+"_sum_diff_gn.dat"
        grid_gamma_nu.save_sum_diff_to_file(f_out)
        print("2D profile for 'up'+'down' and 'up'-'down' in gamma-nu coordinates is saved in the file '{:}'".format(f_out))


        
    if (("3" in ljob)|("4" in ljob)|("5" in ljob)|("6" in ljob)):
        print("\nConverting data to ttheta-phi coordinates")
        grid_ttheta_phi = xml_1.calc_in_grid_ttheta_phi()
    
    if ("3" in ljob):
        f_out = f_base+"_u_d_tp.dat"
        grid_ttheta_phi.save_u_d_to_file(f_out)
        print("2D profile for 'up' and 'down' in ttheta_phi coordinates is saved in the file '{:}'".format(f_out))
    if ("4" in ljob):
        f_out = f_base+"_sum_diff_tp.dat"
        grid_ttheta_phi.save_sum_diff_to_file(f_out)
        print("2D profile for 'up'+'down' and 'up'-'down' in ttheta_phi coordinates is saved in the file '{:}'".format(f_out))
    if ("5" in ljob):
        f_out = f_base+"_u_d_1d.dat"
        grid_ttheta_phi.save_1D_u_d(f_out)
        print("1D profile for 'up' and 'down' is saved in the file '{:}'".format(f_out))
    if ("6" in ljob):
        f_out = f_base+"_sum_diff_1d.dat"
        grid_ttheta_phi.save_1D_sum_diff(f_out)
        print("1D profile for 'up'+'down' and 'up'-'down' is saved in the file '{:}'".format(f_out))
    return 
    

def mainos(larg):
    fdir = os.getcwd()
    lfiles = os.listdir(fdir)
    lxml = [hh for hh in lfiles if (hh.split(".")[-1] == "xml")]
    if len(lxml) == 0:
        print("Can not find any '.xml' file in the folder:\n{:}\n\nProgram is terminated".format(fdir))
        return
    print("The following '.xml' files have been found in the folder:\n{:}\n".format(fdir))
    for ihh, hh in enumerate(lxml):
        print(" ({:}) {:}".format(ihh+1, hh))

    if len(lxml) > 1:        
        print("""\nWith which files you would like to work?
(choose one, several or all of them: 1-3, 5 or all)""")
        reponse = raw_input("........  ")
        reponse=reponse.strip()
        if (reponse.startswith("quit")|reponse.startswith("exit")):
            return
        l_xml = []
        if reponse=="all":
            l_xml = [hh for hh in lxml]
        else:
            lind = []
            lhelp = reponse.split(",")
            for hh in lhelp:
                if hh.strip().isdigit():
                    lind.append(int(hh)-1)
                else:
                    lhh_2 = hh.split("-")
                    num_1 = int(lhh_2[0])-1
                    num_2 = int(lhh_2[1])
                    lind.extend(range(num_1,num_2))
            l_xml = [lxml[ind] for ind in lind]
        print("\nYou choose: \n"+", ".join(l_xml))
        
        flag_xml_sum = False
        if len(l_xml) > 1:
            print("Do you want to sum choosen files?\n(yes/no)")
            reponse = raw_input("........  ")
            reponse=reponse.strip()
            if (reponse.startswith("quit")|reponse.startswith("exit")):
                return
            if reponse.startswith("yes"):
                flag_xml_sum = True
        
        
    else:
        l_xml = [hh for hh in lxml]
        flag_xml_sum = False
    
    print("""\nWhat do you want to do?
    
(1) calculate 2D profile 'up' and 'down' in gamma-nu coordinates
(2)               'up'+'down' and 'up'-'down' 

(3) calculate 2D profile 'up' and 'down' in ttheta-phi coordinates
(4)               'up'+'down' and 'up'-'down' 

(5) calculate 1D profile 'up' and 'down' in ttheta-phi coordinates
(6)               'up'+'down' and 'up'-'down' 

(several answers can be given: 1,3,5,6 or all)
    """)
    reponse = raw_input("........  ")
    reponse=reponse.strip()
    if (reponse.startswith("quit")|reponse.startswith("exit")):
        return
    ljob = []    
    if reponse=="all":
            ljob.extend(["1","2","3","4","5","6"])
    else:
        lind = []
        lhelp = reponse.split(",")
        for hh in lhelp:
            if hh.strip().isdigit():
                ljob.append(hh.strip())
            else:
                lhh_2 = hh.split("-")
                num_1 = int(lhh_2[0])
                num_2 = int(lhh_2[1])+1
                lind = range(num_1,num_2)
                ljob.extend(["{:}".format(hh) for hh in lind])
    
    if (flag_xml_sum):
        print("\nSumming profiles")
        f_name = l_xml[0]
        print("\nFile: {:}".format(f_name))
        xml_1 = XML()
        xml_1.read_file(f_name)
        xml_1.print_frame_info()
        for f_name in l_xml[1:]:
            print("\nFile: {:}".format(f_name))
            xml_2 = XML()
            xml_2.read_file(f_name)
            xml_2.print_frame_info()
            xml_1 = xml_1 + xml_2 
        print("\nInformation about summed frames")
        xml_1.print_frame_info()
        to_do_job(xml_1, ljob, l_xml[0])
    else:
        for f_name in l_xml:
            print("File: {:}".format(f_name))
            xml_1 = XML()
            xml_1.read_file(f_name)
            xml_1.print_frame_info()
            to_do_job(xml_1, ljob, f_name)
            
if (__name__ == "__main__"):
    print(80*"*")
    print("Program {} merges experimental data to 1D or 2D matrix".format('intPSD.py'))
    print(80*"*")
    mainos(sys.argv)
    print(80*"*")
    print("It is done.")
    print(80*"*")
    
