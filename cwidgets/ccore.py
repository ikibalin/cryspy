"""
define main class cagent

subclasses:
    ag_ph, ag_epx, ag_ph_at, ag_exp_ph, ag_ref

it is contained the information which are presented in the parameter xml file.
It is a bridge between GUI and core parameters

"""

import os
import sys
import traceback

import math
import cmath

import numpy
import numpy.linalg

import scipy
import scipy.optimize

#it is not correct in general case when start from any folder !!! SHOULD BE CORRECTED
fdirprog = os.getcwd()

import xml
import xml.etree.ElementTree


import cfunc
import xml2var

def set_val_by_link(obj, slink, val):
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
    if (sname_atr_obj in obj.__dict__.keys()):
        if (len(lhelp)==2):
            sval_el_obj = lhelp[1]
            for el_obj in obj.__dict__[sname_atr_obj]:
                if (sname_el_obj in el_obj.__dict__.keys()):
                    if el_obj.__dict__[sname_el_obj] == sval_el_obj:
                        val_new, lmessage_new = set_val_by_link(el_obj, slink[(num2+1):], val)
                        lmessage.extend(lmessage_new)
                        break
                else:
                    lmessage.append("there is no parameter '{:}' in the object".format(sname_el_obj))
            if (val_new == None)&(lmessage_new == []):
                lmessage.append("the name '{:}' for '{:}' not found".format(sval_el_obj, sname_atr_obj))
        elif (len(lhelp)==1):
            val_old = obj.__dict__[sname_atr_obj]
            if (val == None):
                val_new = val_old
            elif (not(isinstance(val, float))):
                val_new = [hh for hh in val]
                obj.__dict__[sname_atr_obj] = val_new
            else:
                val_new = [val, val_old[1], val_old[2]]
                obj.__dict__[sname_atr_obj] = val_new
        else:
            lmessage.append("unknown link {:} for the object".format(slink))
    else:
        lmessage.append("there is no parameter '{:}' in the object".format(sname_atr_obj))
    return val_new, lmessage


class ccore(object):
    """
    have whole information for input parameter file and some output parameters
    """
    def __init__(self, cagent = None):
        self.exp= []
        self.ph = []
        self.ref = [ccore_ref()]
        self.fdirxml = None
        self.fdirprog = None
        self.fxml = None
        self.valchi2 = None
        self.n_total = None
        if cagent != None:
            self.take_from_agent(cagent)

    def set_fdirxml(self, fdirxml):
        self.fdirxml = fdirxml
        for exp in self.exp:
            exp.set_fdirxml(fdirxml)
        for ref in self.ref:
            ref.set_fdirxml(fdirxml)
            
    def set_fxml(self, fxml):
        self.fxml = fxml
        for exp in self.exp:
            exp.set_fxml(fxml)
            
    def set_fdirprog(self, fdirprog):
        self.fdirprog = fdirprog
        for exp in self.exp:
            exp.set_fdirprog(fdirprog)
        for ph in self.ph:
            ph.set_fdirprog(fdirprog)
            
    def take_from_agent(self, cagent):
        """
        take information from cagent
        """
        self.exp = [ccore_exp(cexp) for cexp in cagent.exp]
        self.ph = [ccore_ph(cph) for cph in cagent.ph]
        self.ref = [ccore_ref(cref) for cref in cagent.ref]
        self.fdirxml = cagent.fdirxml
        self.fdirprog = cagent.fdirprog
        self.fxml = cagent.fxml

    def run_refinement(self):
        """
        refinement over parameters
        """
        lmessage = []
        for cph in self.ph:
            lmessage = cph.load_exp_data()
            if lmessage != []:
                return lmessage
        for cexp in self.exp:
            lmessage = cexp.load_exp_data(self.ph)
            if lmessage != []:
                return lmessage
        lpar_const = self.take_param_const()
        lpar_ref = self.take_param_ref()
        lmessage = self.reload_const(lpar_const)
        if lmessage != []:
            return lmessage

        def tempfunc(lparam):
            self.load_param_to_model(lpar_ref, lparam, lpar_const)
            self.calc_chi2()
            return self.valchi2

        lpar_ref_val = []
        for hh in lpar_ref:
            val, lmessage = self.set_val_by_link(hh, None)
            lpar_ref_val.append(val[0])
            if lmessage!= []:
                return lmessage

        if self.ref[0].refin:
            res = scipy.optimize.minimize(tempfunc, lpar_ref_val, method='Nelder-Mead')
            print "refinement is finished"
            lparam = [float(hh) for hh in res.x]
            lmessage = self.load_param_to_model(lpar_ref, lparam, lpar_const)
            if lmessage != []:
                return lmessage
            lhelp = [list(hh) for hh in res.final_simplex[0]]
            lhelp = cfunc.transposeM(lhelp)
            lpar_ref_val = [hh for hh in lparam]
            #ldeltaparam = [max(hh)-min(hh) for hh in lhelp]
            #self.guiLstatus.setText("Minimized (chi2/N)**0.5 for all experiments is {:.3f}".format((res.fun*1./self.n_total)**0.5))
        else:
            tempfunc(lpar_ref_val)
            #ldeltaparam = [hh*0.01 if (hh != 0) else 0.01 for hh in lpar_ref_val]
            #self.guiLstatus.setText("Minimized (chi2/N)**0.5 for all experiments is {:.3f} (no refinement)".format((self.valchi2*1./self.n_total)**0.5))
        for cexp in self.exp:
            cexp.save_mod_data(self.ph)

        cref = self.ref[0]
        cref.paramreferrors = None
        cref.hessian = None
        cref.ihessian = None
        cref.corrM = None
        cref.paramref = lpar_ref
        cref.paramrefval = lpar_ref_val
        if cref.sigmas:
            self.run_errors()
        self.out_to_string()

        return lmessage

    def run_errors(self):
        """
        numerical calculation of hessian and errorbars
        """
        for cexp in self.exp:
            cexp.load_exp_data(self.ph)


        lpar_ref = self.take_param_ref()
        lpar_ref_val = []
        for hh in lpar_ref:
            val, lmessage = self.set_val_by_link(hh, None)
            lpar_ref_val.append(val[0])
            if lmessage != []:
                return lmessage

        ldpar = [abs(hh)*0.1 if (hh != 0) else 0.01 for hh in lpar_ref_val]
        def tempfunc(lpar):
            self.load_param_to_model(self.take_param_ref(), lpar, self.take_param_const())
            self.calc_chi2()
            return self.valchi2

        #lHessian = [[0. for hh2 in lRefVal] for hh1 in lRefVal]
        lhes = []
        chi = tempfunc(lpar_ref_val)
        n_total = 0
        for cexp in self.exp:
            n_total += cexp.n_total
        self.n_total = n_total
        for ip1 in range(len(ldpar)):
            lhes.append([])
            for ip2 in range(len(ldpar)):
                lpar = [hh for hh in lpar_ref_val]
                chi_help = 0.
                if ip1 != ip2:
                    lpar[ip1] += ldpar[ip1]
                    lpar[ip2] += ldpar[ip2]
                    chi_help += tempfunc(lpar)

                    lpar[ip2] -= 2*ldpar[ip2]
                    chi_help -= tempfunc(lpar)

                    lpar[ip1] -= 2*ldpar[ip1]
                    chi_help += tempfunc(lpar)

                    lpar[ip2] += 2*ldpar[ip2]
                    chi_help -= tempfunc(lpar)

                    sder = chi_help*1./(4*ldpar[ip1]*ldpar[ip2])
                    #lHessian[ipar1][ipar2] = sder
                else:
                    lpar[ip1] += ldpar[ip1]
                    chi_help += tempfunc(lpar)

                    lpar[ip1] -= 2*ldpar[ip1]
                    chi_help += tempfunc(lpar)

                    sder = (chi_help-2.*chi)*1./(ldpar[ip1]**2)
                    #lHessian[ipar1][ipar1] = sder
                #i do not know what should be there, actually something like sder, but for this case the errorbars is to low
                lhes[-1].append(sder*n_total*1./chi)
        chi = tempfunc(lpar_ref_val)
        print "Hessian matrix is"
        for hh in lhes:
            print "".join(["{:9.1f} ".format(hh2) for hh2 in hh])
        try:
            hes_inv = numpy.linalg.inv(lhes)
            lhes_inv = [[hh2 for hh2 in hh1] for hh1 in hes_inv]
            print "\nInverse Hessian matrix is"
            for hh in lhes_inv:
                print "".join(["{:9.5f} ".format(hh2) for hh2 in hh])
            lcorr = cfunc.calcCorrM(lhes_inv)
            print "\nCorrelation matrix is"
            for hh in lcorr:
                print "".join(["{:6.3f} ".format(hh2) for hh2 in hh])
            lerrors = [(abs(hh[ip1]))**0.5 for ip1, hh in enumerate(lhes_inv)]
            for cref in self.ref:
                cref.paramref = lpar_ref
                lpar_ref_val = []
                for hh in lpar_ref:
                    val, lmessage = self.set_val_by_link(hh, None)
                    lpar_ref_val.append(val[0])
                    if lmessage!= []:
                        return lmessage
                cref.paramconstr = []
                cref.paramreferrors = lerrors
                cref.hessian = lhes
                cref.ihessian = lhes_inv
                cref.corrM = lcorr
        except:
            print "Mistake in 'run_errors' of 'ccore.py'"
        return []

    def calc_chi2(self):
        """
        method to calc chi2
        """
        #self.SaveXML()
        lchi2,ln_total = [],[]
        for cexp in self.exp:
            cexp.calc_chi2(self.ph)
            lchi2.append(cexp.valchi2)
            ln_total.append(cexp.n_total)
        self.valchi2 = sum(lchi2)
        self.n_total = sum(ln_total)
        #print sum(lchi2)

    def out_to_string(self):
        """
        form the output file
        """
        lsinfo = []
        lsinfo.append(70*"*")
        lsinfo.append("Information about phases")
        for cph in self.ph:
            lsinfo.extend(cph.out_to_string())
        lsinfo.append(70*"*")
        lsinfo.append("Information about experiments")
        for cexp in self.exp:
            lsinfo.extend(cexp.out_to_string())
        for cref in self.ref:
            lsinfo.extend(cref.out_to_string())
        fout = self.ref[0].output
        fid = open(os.path.join(self.fdirxml, fout), 'w')
        fid.write("\n".join(lsinfo))
        fid.close()

    def set_val_by_link(self, slink, val):
        """
        slink = (experiment,PND),(exp_ph,Tb),(scale)
        """
        val_new, lmessages = set_val_by_link(self, slink, val)
        return val_new, lmessages

    def take_param_ref(self):
        """
        take parameters which are refined
        """
        lparam_ref = []
        for cexp in self.exp:
            lparam_ref.extend(cexp.take_param_ref())
        for cph in self.ph:
            lparam_ref.extend(cph.take_param_ref())
        return lparam_ref

    def take_param_const(self):
        """
        take parameters which are constrained
        """
        lparam_const = []
        for cexp in self.exp:
            lparam_const.extend(cexp.take_param_const())
        for cph in self.ph:
            lparam_const.extend(cph.take_param_const())
        return lparam_const

    def reload_const(self, lpar_const):
        """
        reload constrainted values from other parameters
        """
        lmessage_new = []
        for slink in lpar_const:
            lsconst, lmessage = self.set_val_by_link(slink, None)
            if lmessage != []:
                return lmessage
            sconst = lsconst[2]
            sfunc = sconst.split("[")[0]
            lparams_help = sconst.split("[")[1:]
            lparams = [(hh1.split("]")[0]).strip() for hh1 in lparams_help]
            for iparam, slink2 in enumerate(lparams):
                val, lmessage = self.set_val_by_link(slink2, None)
                lmessage_new.extend(lmessage)
                if lmessage != []:
                    return lmessage
                exec("x{} = {}".format(iparam+1, val[0]))
            valconst = eval(sfunc)
            val_new, lmessage = self.set_val_by_link(slink, valconst)
            lmessage_new.extend(lmessage)
        return lmessage_new

    def load_param_to_model(self, lpar_ref, lpar_ref_val, lpar_const):
        """
        load parameters to model
        """
        for slink, val in zip(lpar_ref, lpar_ref_val):
            val_new, lmessage = self.set_val_by_link(slink, val)
            if lmessage != []:
                return lmessage
        lmessage_new = self.reload_const(lpar_const)
        return lmessage_new




class ccore_exp(object):
    """
    have information of input parameter for experiment
    """
    def __init__(self, cag_exp = None):
        self.name = ""
        self.input = ""
        self.output = ""
        self.bkgr = ""
        self.tth_min = "0.1"
        self.tth_max = "179.9"
        self.phi_min = "-89.9"
        self.phi_max = "89.9"
        self.powder = True
        self.mode_2dpd = False 
        self.pup = [0., False, ""]
        self.pdown = [0., False, ""]
        self.U = [0., False, ""]
        self.V = [0., False, ""]
        self.W = [0., False, ""]
        self.x = [0., False, ""]
        self.y = [0., False, ""]
        self.zeroshift = [0.,False,""]
        self.zshift_a  = [0.,False,""]
        self.p_asym_1  = [0.,False,""]
        self.p_asym_2  = [0.,False,""]
        self.p_asym_3  = [0.,False,""]
        self.p_asym_4  = [0.,False,""]
        self.excl2theta = "(0.0,1.0)"
        self.excl2theta_flags = []
        self.exp_ph = []
        self.modechi2_up = True
        self.modechi2_down = True
        self.modechi2_diff = True
        self.valchi2 = 0.
        self.valchi2up = 0.
        self.valchi2down = 0.
        self.valchi2diff = 0
        self.n_total = 0
        self.fdirprog = None
        self.fdirxml = None
        self.fxml = None
        if cag_exp != None:
            self.take_from_agent(cag_exp)
            
    def set_fdirxml(self, fdirxml):
        self.fdirxml = fdirxml
            
    def set_fxml(self, fxml):
        self.fxml = fxml
            
    def set_fdirprog(self, fdirprog):
        self.fdirprog = fdirprog
            
    def take_from_agent(self, cag_exp):
        self.name = cag_exp.name
        self.input = cag_exp.input
        self.output = cag_exp.output
        self.bkgr = cag_exp.bkgr
        self.powder = cag_exp.powder
        self.pup = [hh for hh in cag_exp.pup]
        self.pdown = [hh for hh in cag_exp.pdown]
        self.U = [hh for hh in cag_exp.U]
        self.V = [hh for hh in cag_exp.V]
        self.W = [hh for hh in cag_exp.W]
        self.x = [hh for hh in cag_exp.x]
        self.y = [hh for hh in cag_exp.y]
        self.zeroshift = [hh for hh in cag_exp.zeroshift]
        self.zshift_a = [hh for hh in cag_exp.zshift_a]
        self.p_asym_1 = [hh for hh in cag_exp.p_asym_1]
        self.p_asym_2 = [hh for hh in cag_exp.p_asym_2]
        self.p_asym_3 = [hh for hh in cag_exp.p_asym_3]
        self.p_asym_4 = [hh for hh in cag_exp.p_asym_4]
        self.excl2theta = cag_exp.excl2theta
        self.exp_ph = [ccore_exp_ph(hh) for hh in cag_exp.exp_ph]
        self.modechi2 = [cag_exp.modechi2_up, cag_exp.modechi2_down, cag_exp.modechi2_diff]
        self.fdirprog = cag_exp.fdirprog
        self.fdirxml = cag_exp.fdirxml
        self.fxml = cag_exp.fxml
        self.tth_min = cag_exp.tth_min
        self.tth_max = cag_exp.tth_max
        self.phi_min = cag_exp.phi_min
        self.phi_max = cag_exp.phi_max
        self.mode_2dpd = cag_exp.mode_2dpd 

    def load_exp_data(self,lcphase):
        """
        make prelimiary calculations and load experimental points before refienement
        """
        lmessage = []
        try:
            if ((self.powder)&(self.mode_2dpd)):
                #2d dimensional data
                ll_int, l_tth, l_phi, dinput = cfunc.read_2dmatrices_np(os.path.join(self.fdirxml,self.input))
            else:
                #1d dimensional data
                dinput=cfunc.readexpdata(os.path.join(self.fdirxml,self.input))
            wavelength = dinput['wavelength']
            self.field=dinput['field']
        except:
            lmessage = ["Can not correctly read experimental data in file '{:}'".format(os.path.join(self.fdirxml,self.input))]
            return lmessage
        self.wavelength=wavelength

        if self.powder:
            
            try:
                dbkgr=cfunc.readexpdata(os.path.join(self.fdirxml,self.bkgr))
                ltthetaBKGRin  = dbkgr["ttheta"]
                lIntBKGRin     = dbkgr["IntBKGR"]
            except:
                lmessage=["Can not correctly read background file'{:}'".format(os.path.join(self.fdirxml,self.bkgr))]
                return lmessage

            if self.mode_2dpd:
                np_int_e_u = numpy.array(ll_int[0], float)
                np_int_e_su = numpy.array(ll_int[1], float)
                np_int_e_d = numpy.array(ll_int[2], float)
                np_int_e_sd = numpy.array(ll_int[3], float)
                np_int_e_diff = numpy.array(ll_int[4], float)
                np_int_e_sdiff = numpy.array(ll_int[5], float)
                    
                    
                phi_min, phi_max = float(self.phi_min), float(self.phi_max)
                tth_min, tth_max = float(self.tth_min), float(self.tth_max) 
                    
                    
                np_phi = numpy.array(l_phi[0], float)
                np_phi_flag = ((np_phi > phi_min) & (np_phi < phi_max))
                    
                np_tth = numpy.array(l_tth[0], float)
                np_tth_flag = ((np_tth > tth_min) & (np_tth < tth_max))
                   
                    
                np_int_e_u = np_int_e_u[np_phi_flag][:, np_tth_flag]
                np_int_e_su = np_int_e_su[np_phi_flag][:, np_tth_flag]
                np_int_e_d = np_int_e_d[np_phi_flag][:, np_tth_flag]
                np_int_e_sd = np_int_e_sd[np_phi_flag][:, np_tth_flag]
                np_int_e_diff = np_int_e_diff[np_phi_flag][:, np_tth_flag]
                np_int_e_sdiff = np_int_e_sdiff[np_phi_flag][:, np_tth_flag]
                np_phi = np_phi[np_phi_flag]
                np_tth = np_tth[np_tth_flag]

                np_bkg_1d = numpy.interp(np_tth, ltthetaBKGRin, lIntBKGRin)
                np_bkg_2d = numpy.meshgrid(np_bkg_1d, np_phi)[0]
                
                self.np_int_e_u = np_int_e_u
                self.np_int_e_su = np_int_e_su
                self.np_int_e_d = np_int_e_d
                self.np_int_e_sd = np_int_e_sd
                self.np_int_e_diff = np_int_e_diff
                self.np_int_e_sdiff = np_int_e_sdiff
                self.np_phi = np_phi
                self.np_tth = np_tth                
                self.np_bkg_2d = np_bkg_2d 

                
                sthovlMIN = math.sin(math.radians(0.5*tth_min))/wavelength
                sthovlMAX = math.sin(math.radians(0.5*tth_max))/wavelength
                
            else:

                ltthetaUDin    = dinput["ttheta"]
                lIntExpUPin    = dinput["IntUP"]
                lsIntExpUPin   = dinput["sIntUP"]
                lIntExpDOWNin  = dinput["IntDOWN"]
                lsIntExpDOWNin = dinput["sIntDOWN"]

                shelp=self.excl2theta
                lexcltthRanges=[]
                if (shelp.strip()!=""):
                    for hh1 in shelp.split(")"):
                        if (hh1.strip()!=""):
                            tthmin=float(hh1.split(",")[0].split("(")[1])
                            tthmax=float(hh1.split(",")[1])
                            lexcltthRanges.append((tthmin,tthmax))
                self.valexcl2theta=lexcltthRanges
                lflags=[True for hh in ltthetaUDin]
                for iind,excltth in enumerate(lexcltthRanges):
                    lflagsh=[((hh1<excltth[0])|(hh1>excltth[1])) for hh1 in ltthetaUDin]
                    lflags=[(hh1&hh2) for hh1,hh2 in zip(lflags,lflagsh)]
                self.excl2theta_flags = lflags
                lIntBKGRUD     = [scipy.interp(hh,ltthetaBKGRin,lIntBKGRin) for hh in ltthetaUDin]

                self.tthetaUDin = ltthetaUDin
                self.IntExpUPin = lIntExpUPin
                self.sIntExpUPin = lsIntExpUPin
                self.IntExpDOWNin = lIntExpDOWNin
                self.sIntExpDOWNin = lsIntExpDOWNin
                self.IntBKGRUD = lIntBKGRUD
                
                sthovlMIN = math.sin(math.radians(0.5*min(ltthetaUDin)))/wavelength
                sthovlMAX = math.sin(math.radians(0.5*max(ltthetaUDin)))/wavelength

            for cexp_ph in self.exp_ph:
                lmessage = cexp_ph.load_exp_data(lcphase,sthovlMinMax=[sthovlMIN,sthovlMAX])
                if lmessage!=[]:
                    return lmessage
        else:
            lhh1=dinput['orientation']
            self.mU=[lhh1[:3],lhh1[3:6],lhh1[6:9]]
            lhkl=[(hh1,hh2,hh3) for hh1,hh2,hh3 in zip(dinput['h'],dinput['k'],dinput['l'])]
            self.hkl=lhkl
            self.FRexp=[hh1 for hh1 in dinput['FR']]
            self.sFRexp=[hh1 for hh1 in dinput['sFR']]
            for cphExp in self.exp_ph:
                lmessage = cphExp.load_exp_data(lcphase,listHKL=lhkl)
                if lmessage!=[]:
                    return lmessage
                cphExp.mUB=cfunc.multMAT(self.mU,cphExp.mB)
        return lmessage

    def calc_chi2(self,lcph):
        """
        calculation of chi2 for one experiment
        """
        wavelength = self.wavelength
        field = self.field
        normH = (sum([hh**2 for hh in field]))**0.5
        e_up = [hh*1./normH for hh in field]
        p_up = self.pup[0]
        p_down = self.pdown[0]

        if (self.powder):
            if self.mode_2dpd:
                np_int_e_u = self.np_int_e_u
                np_int_e_su = self.np_int_e_su
                np_int_e_d = self.np_int_e_d
                np_int_e_sd = self.np_int_e_sd
                np_int_e_diff = self.np_int_e_diff
                np_int_e_sdiff = self.np_int_e_sdiff
                np_bkg_2d = self.np_bkg_2d
                np_phi = self.np_phi
                np_tth = self.np_tth   
                
                np_phi_rad = np_phi * math.pi/180.
                np_tth_rad = np_tth * math.pi/180.
                    
                np_tth_rad_2d, np_phi_rad_2d = numpy.meshgrid(np_tth_rad, np_phi_rad)
                    
                np_c2_alpha_2d = numpy.cos(0.5*np_tth_rad_2d)**2*numpy.sin(np_phi_rad_2d)**2
                np_s2_alpha_2d = 1.0 - np_c2_alpha_2d 
                
                #vertical contribution to the Lorentz factor
                np_lor_vert = (1.0-(numpy.sin(np_tth_rad_2d)*numpy.sin(np_phi_rad_2d))**2)**0.5
                
                np_int_m_u = numpy.zeros(np_s2_alpha_2d.shape, float)
                np_int_m_d = numpy.zeros(np_s2_alpha_2d.shape, float)      
                
            else:
                ltthetaUDin = self.tthetaUDin
                lIntExpUPin = self.IntExpUPin
                lsIntExpUPin = self.sIntExpUPin
                lIntExpDOWNin = self.IntExpDOWNin
                lsIntExpDOWNin = self.sIntExpDOWNin
                lIntBKGRUD = self.IntBKGRUD
            zeroshift = self.zeroshift[0]
            zshift_a = self.zshift_a[0]
            p_asym = [self.p_asym_1[0], self.p_asym_2[0], self.p_asym_3[0], self.p_asym_4[0]]
            lIntUPModExp, lIntDOWNModExp = [], []
        else:
            mU = self.mU
            lhkl = self.hkl
            lFRexp = self.FRexp
            lsFRexp = self.sFRexp

        for cexp_ph in self.exp_ph:
            cphmod = [cph for iphase,cph in enumerate(lcph) if cph.name==cexp_ph.name][0]
            ucp = [cphmod.cella[0],cphmod.cellb[0],cphmod.cellc[0],cphmod.cellalpha[0],cphmod.cellbeta[0],cphmod.cellgamma[0]]
            scale = cexp_ph.scale[0]
            lhkl = cexp_ph.hkl
            lkloc = cexp_ph.kloc
            ext = [cexp_ph.ext11[0], cexp_ph.ext22[0], cexp_ph.ext33[0], cexp_ph.ext12[0], cexp_ph.ext13[0], cexp_ph.ext23[0]]
            #for hkl,kloc in zip(lhkl,lkloc):
            #    print "hkl: {0[0]:3} {0[1]:3} {0[2]:3}   klocl: {1[0]:6.4f} {1[1]:6.4f} {1[2]:6.4f}".format(hkl,kloc)
            #    quit()
            lsthovlhkl = [cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]
            lFNhkl = cphmod.calcFNhkl(lhkl)

            lSFThkl = cphmod.calcSFThkl(lhkl)
            lthetahkl = [math.asin((wavelength*hh1)) if wavelength*hh1<1 else 0.5*math.pi for hh1 in lsthovlhkl]#sin2(theta_hkl)
            ltthetahkl = [2.*math.degrees(hh1) for hh1 in lthetahkl]

            if (self.powder):
                lmulthkl = cexp_ph.multhkl
                UVWIgxy = (self.U[0], self.V[0], self.W[0], cexp_ph.igsize[0], self.x[0], self.y[0])
                lmT = [cfunc.calcmT(kloc) for kloc in lkloc]
                lmTT = [cfunc.transposeM(mT) for mT in lmT]
                lTHETA = [cfunc.multMAT(cfunc.multMAT(mTT,SFT),mT) for mT,mTT,SFT in zip(lmT,lmTT,lSFThkl)]
                lTHETA11 = [THETA[0][0] for THETA in lTHETA]
                lTHETA22 = [THETA[1][1] for THETA in lTHETA]
                lTHETA12 = [THETA[0][1] for THETA in lTHETA]
                lMSFpsq = [(normH**2)*abs(0.5*(THETA11*THETA11.conjugate()+THETA22*THETA22.conjugate())+THETA12*THETA12.conjugate()) for THETA11,THETA22,THETA12 in zip(lTHETA11,lTHETA22,lTHETA12)]
                lMSFp_field = [normH*0.5*(THETA11+THETA22) for THETA11,THETA22,THETA12 in zip(lTHETA11,lTHETA22,lTHETA12)]
                lcross = [2.*(hh1.real*hh2.real+hh1.imag*hh2.imag) for hh1,hh2 in zip(lFNhkl,lMSFp_field)]
                
                if (not(self.mode_2dpd)):
                    #lNUCL=[abs(hh1)**2 for hh1 in lFNhkl]
                    lIntModUPhkl, lIntModDOWNhkl = [], []

                    #for Ising model, not tested
                    #lFMhkl_sq_aver, lFMhkl_perp_aver = cphmod.calcFMhkl_ising_powder(lhkl,lmTT)
                    #cexp_ph.FMhkl_sq_aver = lFMhkl_sq_aver
                    #cexp_ph.FMhkl_perp_aver = lFMhkl_perp_aver

                    #for hkl, fn, fm_perp_eup, fm_perp_sq, multhkl, fm2_sq, fm2_perp in zip(lhkl, lFNhkl, lMSFp_field, lMSFpsq, lmulthkl, lFMhkl_sq_aver, lFMhkl_perp_aver):
                    for hkl, fn, fm_perp_eup, fm_perp_sq, multhkl in zip(lhkl, lFNhkl, lMSFp_field, lMSFpsq, lmulthkl):
                        #I_p, I_m = cfunc.calc_int_extinc_powder(hkl, fn, fm2_perp, float(fm2_sq), ext, p_up, p_down, ucp, wavelength)
                        I_p, I_m = cfunc.calc_int_extinc_powder(hkl,fn,fm_perp_eup,fm_perp_sq,ext,p_up,p_down,ucp,wavelength)
                        lIntModUPhkl.append(scale*multhkl*I_p)
                        lIntModDOWNhkl.append(scale*multhkl*I_m)
                        
                    cexp_ph.IntModUPhkl = lIntModUPhkl
                    cexp_ph.IntModDOWNhkl = lIntModDOWNhkl
                    lIntModUPinPhase = cfunc.calcprofile(ltthetahkl,lIntModUPhkl,ltthetaUDin,zeroshift,UVWIgxy,zshift_a,p_asym)
                    lIntModDOWNinPhase = cfunc.calcprofile(ltthetahkl,lIntModDOWNhkl,ltthetaUDin,zeroshift,UVWIgxy,zshift_a,p_asym)
                    if (lIntUPModExp != []):
                        lIntUPModExp = [hh1+hh2 for hh1,hh2 in zip(lIntUPModExp,lIntModUPinPhase)]
                    else:
                        lIntUPModExp = [hh for hh in lIntModUPinPhase]
                    if (lIntDOWNModExp != []):
                        lIntDOWNModExp = [hh1+hh2 for hh1,hh2 in zip(lIntDOWNModExp,lIntModDOWNinPhase)]
                    else:
                        lIntDOWNModExp = [hh for hh in lIntModDOWNinPhase]

                else:
                    lTHETA13 = [THETA[0][2] for THETA in lTHETA]
                    lTHETA23 = [THETA[1][2] for THETA in lTHETA]
                    

                    lnp_shape_hkl = []
                    for hkl, ttheta_hkl in zip(lhkl, ltthetahkl):    
                        np_shape_hkl = cfunc.calc_np_shape_hkl(np_tth, ttheta_hkl, zeroshift, UVWIgxy, zshift_a, p_asym)
                        lnp_shape_hkl.append(np_shape_hkl)
                    
                    

                    
                    for fn, fm_perp_e_up, fm_perp_sq, tth_13, tth_23, multhkl, np_shape_hkl in zip(lFNhkl, lMSFp_field, lMSFpsq, lTHETA13, lTHETA23, lmulthkl, lnp_shape_hkl):    
                        np_val_1 = abs(fn)**2 * numpy.ones(np_s2_alpha_2d.shape, float)
                        np_val_2 = ((fm_perp_e_up*fn.conjugate()+fm_perp_e_up.conjugate()*fn).real) * np_s2_alpha_2d
                        np_val_3 = fm_perp_sq * np_s2_alpha_2d  + 0.5*normH**2 * abs(tth_13*tth_13.conjugate() + tth_23*tth_23.conjugate()) * np_c2_alpha_2d
                        #just test
                        #np_val_3 = fm_perp_sq * np_s2_alpha_2d**2  + 0.5*normH**2 * abs(tth_13*tth_13.conjugate() + tth_23*tth_23.conjugate()) * np_c2_alpha_2d
                        #without magnetic field
                        #np_val_2 = numpy.zeros(np_s2_alpha_2d.shape, float)
                        #np_val_3 = fm_perp_sq * numpy.ones(np_s2_alpha_2d.shape, float)
                        

                        np_shape_hkl_2d = numpy.meshgrid(np_shape_hkl, np_phi_rad)[0]
                        #no extinction are taken into account
                        np_int_m_u += scale*multhkl*(np_val_1+p_up*np_val_2+np_val_3)*np_shape_hkl_2d
                        np_int_m_d += scale*multhkl*(np_val_1-p_down*np_val_2+np_val_3)*np_shape_hkl_2d

                #correction on the vertical component of the Lorentz factor
                np_int_m_u = np_int_m_u*np_lor_vert
                np_int_m_d = np_int_m_d*np_lor_vert
                
                    
                cexp_ph.mT = lmT
                cexp_ph.THETA = lTHETA
                cexp_ph.MSFpsq = lMSFpsq
                cexp_ph.MSFp = lMSFp_field
                cexp_ph.cross = lcross
            else:

                phiD, chiD, omegaD = math.radians(0.), math.radians(0.), math.radians(0.)
                PHId = [[ math.cos(phiD),math.sin(phiD),0.],
                        [-math.sin(phiD),math.cos(phiD),0.],
                        [             0.,            0.,1.]]
                OMEGAd = [[ math.cos(omegaD),math.sin(omegaD),0.],
                          [-math.sin(omegaD),math.cos(omegaD),0.],
                          [               0.,              0.,1.]]
                CHId = [[ math.cos(chiD),0.,math.sin(chiD)],
                        [             0.,1.,            0.],
                        [-math.sin(chiD),0.,math.cos(chiD)]]
                mDU = cfunc.multMAT(OMEGAd, cfunc.multMAT(CHId, cfunc.multMAT(PHId,mU)))

                Hloc = cfunc.multMv(cfunc.transposeM(mDU),field)  #magnetic field in local coordinate system
                e_up_loc = cfunc.multMv(cfunc.transposeM(mDU),e_up)
                Puloc = [p_up*hh1*1./normH for hh1 in Hloc]
                Pdloc = [-1.*p_down*hh1*1./normH for hh1 in Hloc]

                lMSFloc = [cfunc.multMv(SFT,Hloc) for SFT in lSFThkl]#magnetic moment in local coordinate system
                lMSFploc = [cfunc.vecprod(cfunc.vecprod(kloc,msfloc),kloc) for kloc,msfloc in zip(lkloc,lMSFloc)]#its perpendicular component

                lFMloc = cphmod.calcFMhkl_ising_mono(lhkl, Hloc)
                lFMploc = [cfunc.vecprod(cfunc.vecprod(kloc,msfloc),kloc) for kloc,msfloc in zip(lkloc,lFMloc)]#its perpendicular component

                #tensor aproach
                #lFRmod_ext = [cfunc.calc_fr_extinc(hkl,fn,fm_perp,e_up_loc,ext,p_up,p_down,ucp,wavelength) for hkl,fn,fm_perp in zip(lhkl,lFNhkl,lMSFploc)]
                #lFRmod = [cfunc.calcFR(FN,msfp,Puloc,Pdloc) for FN,msfp in zip(lFNhkl,lMSFploc)]

                #Ising model
                #lFRmod_ext = [cfunc.calc_fr_extinc(hkl,fn,fm_perp,e_up_loc,ext,p_up,p_down,ucp,wavelength) for hkl,fn,fm_perp in zip(lhkl,lFNhkl,lFMploc)]
                #lFRmod = [cfunc.calcFR(FN,msfp,Puloc,Pdloc) for FN,msfp in zip(lFNhkl,lFMploc)]


                #tensor + Ising (it should be redone)
                lFRmod_ext = [cfunc.calc_fr_extinc(hkl, fn, [hh1+hh2 for hh1, hh2 in zip(msfp, fm_perp)], e_up_loc, ext, p_up, p_down, ucp, wavelength) for hkl, fn, msfp, fm_perp in zip(lhkl, lFNhkl, lMSFploc, lFMploc)]
                lFRmod = [cfunc.calcFR(FN, [hh1+hh2 for hh1, hh2 in zip(msfp, fm_perp)], Puloc, Pdloc) for FN, msfp, fm_perp in zip(lFNhkl, lMSFploc, lFMploc)]

                lkglob = [cfunc.multMv(mDU,kloc) for kloc in lkloc]#unity vector k_s - k_i in global coordinate system
                lMSFhkl = [cfunc.multMv(mDU,msfloc) for msfloc in lMSFloc]#magnetic moment in global coordinate system
                lMSFperp = [cfunc.multMv(mDU,msfploc) for msfploc in lMSFploc]#its perpendicular component

                lFMhkl = [cfunc.multMv(mDU,msfloc) for msfloc in lFMloc]
                lFMphkl = [cfunc.multMv(mDU,msfloc) for msfloc in lFMploc]


                cexp_ph.FRmod = lFRmod
                cexp_ph.FRmod_ext = lFRmod_ext
                cexp_ph.MSFhkl = lMSFhkl
                cexp_ph.MSFperp = lMSFperp
                cexp_ph.kglob = lkglob

                cexp_ph.FMhkl = lFMhkl
                cexp_ph.FMphkl = lFMphkl


            cexp_ph.FNhkl = lFNhkl
            cexp_ph.SFThkl = lSFThkl
            cexp_ph.thetahkl = lthetahkl

        if (self.powder):
            if self.mode_2dpd:
                np_int_m_diff = np_int_m_u - np_int_m_d
                np_int_m_u += np_bkg_2d
                np_int_m_d += np_bkg_2d
                    
                np_chi2_u = ((np_int_e_u - np_int_m_u)/np_int_e_su)**2
                np_chi2_d = ((np_int_e_d - np_int_m_d)/np_int_e_sd)**2
                np_chi2_diff = ((np_int_e_diff - np_int_m_diff)/np_int_e_sdiff)**2
                    
                chi2_u = (np_chi2_u[numpy.logical_not(numpy.isnan(np_chi2_u))]).sum()
                n_u = numpy.logical_not(numpy.isnan(np_chi2_u)).sum()
                chi2_d = (np_chi2_d[numpy.logical_not(numpy.isnan(np_chi2_d))]).sum()
                n_d = numpy.logical_not(numpy.isnan(np_chi2_d)).sum()
                chi2_diff = (np_chi2_diff[numpy.logical_not(numpy.isnan(np_chi2_diff))]).sum()
                n_diff = numpy.logical_not(numpy.isnan(np_chi2_diff)).sum()
                
                self.valchi2Up = chi2_u
                self.valchi2Down = chi2_d
                self.valchi2Diff = chi2_diff
                self.np_int_m_u = np_int_m_u
                self.np_int_m_d = np_int_m_d
                self.np_int_m_diff = np_int_m_diff
                self.valchi2 = chi2_u*self.modechi2_up+chi2_d*self.modechi2_down+chi2_diff*self.modechi2_diff
                self.n_total = n_u*self.modechi2_up + n_d*self.modechi2_down + n_diff*self.modechi2_diff

            else:
                pass
                lflags = self.excl2theta_flags
                lIntUPModExp = [hh1+hh2 for hh1,hh2 in zip(lIntUPModExp,lIntBKGRUD)]
                lIntDOWNModExp = [hh1+hh2 for hh1,hh2 in zip(lIntDOWNModExp,lIntBKGRUD)]
                lIntDIFFModExp = [hh1-hh2 for hh1,hh2 in zip(lIntUPModExp,lIntDOWNModExp)]
                lIntExpDIFFin = [hh1-hh2 for hh1,hh2 in zip(lIntExpUPin,lIntExpDOWNin)]
                chi2UP = sum([((hh1-hh2)*1./hh3)**2 for hh1,hh2,hh3,cond in zip(lIntExpUPin,lIntUPModExp,lsIntExpUPin,lflags) if (cond)])
                chi2DOWN = sum([((hh1-hh2)*1./hh3)**2 for hh1,hh2,hh3,cond in zip(lIntExpDOWNin,lIntDOWNModExp,lsIntExpDOWNin,lflags) if (cond)])
                chi2DIFF = sum([(((hh1-hh2)**2)*1./(hh3**2+hh4**2)) for hh1,hh2,hh3,hh4 in zip(lIntExpDIFFin,lIntDIFFModExp,lsIntExpUPin,lsIntExpDOWNin)])
                #chi2UP *= 1./len(lIntExpUPin)
                #chi2DOWN *= 1./len(lIntExpDOWNin)
                #chi2DIFF *= 1./len(lIntExpUPin)
                self.valchi2Up = chi2UP
                self.valchi2Down = chi2DOWN
                self.valchi2Diff = chi2DIFF
                self.IntUPModExp = lIntUPModExp
                self.IntDOWNModExp = lIntDOWNModExp
                self.IntDIFFModExp = lIntDIFFModExp
                self.valchi2=chi2UP*self.modechi2_up+chi2DOWN*self.modechi2_down+chi2DIFF*self.modechi2_diff
                self.n_total=sum(lflags)*self.modechi2_up+sum(lflags)*self.modechi2_down+len(lflags)*self.modechi2_diff
        else:
            self.valchi2 = sum([((exp-mod)*1./sigma)**2 for exp,mod,sigma in zip(lFRexp,lFRmod_ext,lsFRexp)])
            self.n_total=len(lFRexp)
        #print "chi2 {:.2f}".format(self.valchi2)

    def save_mod_data(self,lphase):
        """
        save the calculated points in the output file
        """
        fdirxml=self.fdirxml
        fout=self.output
        if (self.powder):
            if self.mode_2dpd:
                np_phi = self.np_phi
                np_tth = self.np_tth
                np_int_e_u = self.np_int_e_u
                np_int_m_u = self.np_int_m_u
                np_int_e_su = self.np_int_e_su
                np_int_e_d = self.np_int_e_d
                np_int_m_d = self.np_int_m_d
                np_int_e_sd = self.np_int_e_sd
                np_int_e_diff = self.np_int_e_diff
                np_int_m_diff = self.np_int_m_diff
                np_int_e_sdiff = self.np_int_e_sdiff

                    
                np_int_em_u = np_int_e_u - np_int_m_u
                np_int_em_d = np_int_e_d - np_int_m_d
                np_int_em_diff = np_int_e_diff - np_int_m_diff
                
                finp = os.path.join(fdirxml,"e_"+self.output)
                cfunc.save_2dmatrices_np(finp, np_phi, np_tth, np_int_e_u, np_int_e_su, np_int_e_d, np_int_e_sd, np_int_e_diff, np_int_e_sdiff)
                finp = os.path.join(fdirxml,"m_"+self.output)
                cfunc.save_2dmatrices_np(finp, np_phi, np_tth, np_int_m_u, np_int_e_su, np_int_m_d, np_int_e_sd, np_int_m_diff, np_int_e_sdiff)
                finp = os.path.join(fdirxml,"diff_"+self.output)
                cfunc.save_2dmatrices_np(finp, np_phi, np_tth, np_int_em_u, np_int_e_su, np_int_em_d, np_int_e_sd, np_int_em_diff, np_int_e_sdiff)
                llsgpl = [""]

            else:
                lflags = self.excl2theta_flags
                lIntUPModExp   = self.IntUPModExp
                lIntDOWNModExp = self.IntDOWNModExp
                #lIntDIFFModExp = self.IntDIFFModExp
                ltthetaUDin    = self.tthetaUDin
                lIntExpUPin    = self.IntExpUPin
                lsIntExpUPin   = self.sIntExpUPin
                lIntExpDOWNin  = self.IntExpDOWNin
                lsIntExpDOWNin = self.sIntExpDOWNin
                lIntBKGRUD     = self.IntBKGRUD
                wavelength = self.wavelength
                zeroshift = self.zeroshift[0]
                zshift_a = self.zshift_a[0]
                llsgpl=["    ttheta   IntUPexp  sIntUPexp   IntUPmod IntDOWNexp sIntDOWNexp IntDOWNmod   IntBKGR flag"]
                for ttheta,IntUExp,sIntUExp,IntUMod,IntDExp,sIntDExp,IntDMod,IntBMod,cond in zip(ltthetaUDin,lIntExpUPin,lsIntExpUPin,lIntUPModExp,lIntExpDOWNin,lsIntExpDOWNin,lIntDOWNModExp,lIntBKGRUD,lflags):
                    llsgpl.append("{:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:10.3f} {:2}".format(ttheta,IntUExp,sIntUExp,IntUMod,IntDExp,sIntDExp,IntDMod,IntBMod,1*cond))
                llsgpl.append("\n\n")
                for iphase,cphExp in enumerate(self.exp_ph):
                    llsgpl.append("#phase {}".format(iphase+1))
                    cPhMod=[cph for cph in lphase if cph.name==cphExp.name][0]
                    ucp=[cPhMod.cella[0],cPhMod.cellb[0],cPhMod.cellc[0],
                     cPhMod.cellalpha[0],cPhMod.cellbeta[0],cPhMod.cellgamma[0]]
                    lsthovlhkl=[cfunc.calcsthovl(hkl,ucp) for hkl in cphExp.hkl]
                    ltthetahkl=[math.degrees(2.*math.asin(sthovl*wavelength)) for sthovl in lsthovlhkl]
                    llsgpl.append("  h   k   l mult    sthovl    ttheta")
                    llsgpl.extend(["{0[0]:3} {0[1]:3} {0[2]:3} {1:4} {2:9.5f} {3:9.3f}".format(hkl,mult,sthovl,ttheta+zeroshift+zshift_a*math.cos(0.5*math.radians(ttheta))) for hkl,mult,sthovl,ttheta in zip(cphExp.hkl,cphExp.multhkl,lsthovlhkl,ltthetahkl)])
                    llsgpl.append("\n\n")
        else:
            for iphase,cphExp in enumerate(self.exp_ph):
                llsgpl=["     h    k    l     FRexp    sFRexp    FRmod_ext   FRmod"]
                for hkl,FRmod,FRmod_ext,FRexp,sFRexp in zip(cphExp.hkl,cphExp.FRmod,cphExp.FRmod_ext,self.FRexp,self.sFRexp):
                    llsgpl.append(" {0[0]:4} {0[1]:4} {0[2]:4} {1:9.4f} {2:9.4f}{3:9.4f}{4:9.4f}".format(hkl,FRexp,sFRexp,FRmod_ext,FRmod))
                llsgpl.append("\n\n")
        fid=open(os.path.join(fdirxml,fout),"w")
        fid.write("\n".join(llsgpl))
        fid.close()

    def set_val_by_link(self,slink,val):
        """slink = (U,)
           slink = (exp_ph,Tb),(scale,)
        """
        val_new, lmessages = set_val_by_link(self, slink, val)
        return val_new, lmessages

    def take_param_ref(self):
        lRefParam=[]
        if (self.pup[1]): lRefParam.append("(exp,{})(pup)".format(self.name))
        if (self.pdown[1]): lRefParam.append("(exp,{})(pdown)".format(self.name))
        if (self.powder):
            if (self.U[1]): lRefParam.append("(exp,{})(U)".format(self.name))
            if (self.V[1]): lRefParam.append("(exp,{})(V)".format(self.name))
            if (self.W[1]): lRefParam.append("(exp,{})(W)".format(self.name))
            if (self.x[1]): lRefParam.append("(exp,{})(x)".format(self.name))
            if (self.y[1]): lRefParam.append("(exp,{})(y)".format(self.name))
            if (self.zeroshift[1]): lRefParam.append("(exp,{})(zeroshift)".format(self.name))
            if (self.zshift_a[1]): lRefParam.append("(exp,{})(zshift_a)".format(self.name))
            if (self.p_asym_1[1]): lRefParam.append("(exp,{})(p_asym_1)".format(self.name))
            if (self.p_asym_2[1]): lRefParam.append("(exp,{})(p_asym_2)".format(self.name))
            if (self.p_asym_3[1]): lRefParam.append("(exp,{})(p_asym_3)".format(self.name))
            if (self.p_asym_4[1]): lRefParam.append("(exp,{})(p_asym_4)".format(self.name))
            for cPhExp in self.exp_ph:
                lRefParamPhExp=cPhExp.take_param_ref()
                if (lRefParamPhExp!=[]):
                    lRefParam.extend(["(exp,{}){}".format(self.name,hh) for hh in lRefParamPhExp])
        else:
            if self.exp_ph != []:
                #potentially mistake, correct it
                lRefParamPhExp = self.exp_ph[0].take_param_ref()
                if (lRefParamPhExp!=[]):
                    lRefParam.extend(["(exp,{}){}".format(self.name,hh) for hh in lRefParamPhExp])

        return lRefParam

    def take_param_const(self):
        lRefParam=[]
        if (self.pup[2]!=""): lRefParam.append("(exp,{})(pup)".format(self.name))
        if (self.pdown[2]!=""): lRefParam.append("(exp,{})(pdown)".format(self.name))
        if (self.powder):
            if (self.U[2]!=""): lRefParam.append("(exp,{})(U)".format(self.name))
            if (self.V[2]!=""): lRefParam.append("(exp,{})(V)".format(self.name))
            if (self.W[2]!=""): lRefParam.append("(exp,{})(W)".format(self.name))
            if (self.x[2]!=""): lRefParam.append("(exp,{})(x)".format(self.name))
            if (self.y[2]!=""): lRefParam.append("(exp,{})(y)".format(self.name))
            if (self.zeroshift[2]!=""): lRefParam.append("(exp,{})(zeroshift)".format(self.name))
            if (self.zshift_a[2]!=""): lRefParam.append("(exp,{})(zshift_a)".format(self.name))
            if (self.p_asym_1[2]!=""): lRefParam.append("(exp,{})(p_asym_1)".format(self.name))
            if (self.p_asym_2[2]!=""): lRefParam.append("(exp,{})(p_asym_2)".format(self.name))
            if (self.p_asym_3[2]!=""): lRefParam.append("(exp,{})(p_asym_3)".format(self.name))
            if (self.p_asym_4[2]!=""): lRefParam.append("(exp,{})(p_asym_4)".format(self.name))
            for cPhExp in self.exp_ph:
                lRefParamPhExp=cPhExp.take_param_const()
                if (lRefParamPhExp!=[]):
                    lRefParam.extend(["(exp,{}){}".format(self.name,hh) for hh in lRefParamPhExp])
        else:
            if self.exp_ph != []:
                #potentially mistake, correct it
                lRefParamPhExp=self.exp_ph[0].take_param_const()
                if (lRefParamPhExp!=[]):
                    lRefParam.extend(["(exp,{}){}".format(self.name,hh) for hh in lRefParamPhExp])
        return lRefParam

    def out_to_string(self):
        """
        output information about experiment
        """
        llsout = []
        llsout.append(70*"*")
        smode="powder" if (self.powder) else "mono"
        llsout.append("Experiment {} ({}):".format(self.name,smode))
        if (self.powder):
            #npoints=len(self.tthetaUDin)
            llsout.append("wavelength = {:10.5f}".format(self.wavelength))
            llsout.append("U = {:10.5f}, V = {:10.5f}, W = {:10.5f}".format(self.U[0],self.V[0],self.W[0]))
            llsout.append("x = {:10.5f}, y = {:10.5f}".format(self.x[0],self.y[0]))
            llsout.append("zeroshift = {:10.5f}; zshift_a = {:10.5f};".format(self.zeroshift[0],self.zshift_a[0]))
            llsout.append("asymmetry parameters: {:10.5f} {:10.5f} {:10.5f} {:10.5f}".format(self.p_asym_1[0],self.p_asym_2[0],self.p_asym_3[0],self.p_asym_4[0]))
            llsout.append("polarization of beam:  p_up = {:10.5f}, p_down = {:10.5f}".format(self.pup[0],self.pdown[0]))
            for cexp_ph in self.exp_ph:
                llsout.append(70*"*")
                llsoutPhExp = cexp_ph.out_to_string(self.powder,self.wavelength)
                llsout.extend(llsoutPhExp)
            llsout.append(70*"*")
            llsout.append("Experiment is calculated.\nrchi2 UP: {:.2f} DOWN: {:.2f} DIFF: {:.2f} per point".format((self.valchi2Up)**0.5,(self.valchi2Down)**0.5,(self.valchi2Diff)**0.5))
            llsout.append(70*"*")
        else:
            for cexp_ph in self.exp_ph:
                llsout.append(70*"*")
                llsoutPhExp = cexp_ph.out_to_string(self.powder,self.wavelength,lFRexp=self.FRexp,lsFRexp=self.sFRexp)
                llsout.extend(llsoutPhExp)
            llsout.append("Chi2 is {:.3f}:".format(self.valchi2))
        return llsout

class ccore_exp_ph(object):
    """
    have information of input parameter for phase in the experiment
    """
    def __init__(self, cag_exp_ph = None):
        self.name = " "
        self.scale = [1.00000, False, ""]
        self.igsize = [0.00000, False, ""]
        self.ext11 = [0.00000, False, ""]
        self.ext22 = [0.00000, False, ""]
        self.ext33 = [0.00000, False, ""]
        self.ext12 = [0.00000, False, ""]
        self.ext13 = [0.00000, False, ""]
        self.ext23 = [0.00000, False, ""]
        if cag_exp_ph != None:
            self.take_from_agent(cag_exp_ph)

    def take_from_agent(self, cag_exp_ph):
        """
        load information from cagent
        """
        self.name = cag_exp_ph.name
        self.scale = [hh for hh in cag_exp_ph.scale]
        self.igsize = [hh for hh in cag_exp_ph.igsize]
        self.ext11 = [hh for hh in cag_exp_ph.ext11]
        self.ext22 = [hh for hh in cag_exp_ph.ext22]
        self.ext33 = [hh for hh in cag_exp_ph.ext33]
        self.ext12 = [hh for hh in cag_exp_ph.ext12]
        self.ext13 = [hh for hh in cag_exp_ph.ext13]
        self.ext23 = [hh for hh in cag_exp_ph.ext23]

    def load_exp_data(self,lcphase,sthovlMinMax=None,listHKL=None):
        """
        preliminary calculations before parameters' refinement
        """
        lmessage = []
        cphmod=[cph for iphase,cph in enumerate(lcphase) if cph.name == self.name][0]
        ucp=[cphmod.cella[0],cphmod.cellb[0],cphmod.cellc[0],
             cphmod.cellalpha[0],cphmod.cellbeta[0],cphmod.cellgamma[0]]
        lmessage = cphmod.load_exp_data()
        if lmessage != []:
            return lmessage
        dsymm=cphmod.dsymm
        if (listHKL!=None):
            lhkl=listHKL
        if (sthovlMinMax!=None):
            [sthovlMIN,sthovlMAX]=sthovlMinMax
            lhkl,lmulthkl=cfunc.listhkl(sthovlMIN,sthovlMAX,ucp,dsymm)
            self.multhkl=lmulthkl
        self.hkl=lhkl
        mB=cphmod.mB
        lkloc=[cfunc.calck(hkl,mB) for hkl in lhkl]
        self.mB=mB#the information is doubled, original in cphmod.mB
        self.kloc=lkloc
        self.ucp=ucp#just for output
        return lmessage

    def set_val_by_link(self,slink,val):
        val_new, lmessages = set_val_by_link(self, slink, val)
        return val_new, lmessages

    def take_param_ref(self):
        lparam_ref = []
        if (self.scale[1]): lparam_ref.append("(exp_ph,{})(scale)".format(self.name))
        if (self.igsize[1]): lparam_ref.append("(exp_ph,{})(igsize)".format(self.name))
        if (self.ext11[1]): lparam_ref.append("(exp_ph,{})(ext11)".format(self.name))
        if (self.ext22[1]): lparam_ref.append("(exp_ph,{})(ext22)".format(self.name))
        if (self.ext33[1]): lparam_ref.append("(exp_ph,{})(ext33)".format(self.name))
        if (self.ext12[1]): lparam_ref.append("(exp_ph,{})(ext12)".format(self.name))
        if (self.ext13[1]): lparam_ref.append("(exp_ph,{})(ext13)".format(self.name))
        if (self.ext23[1]): lparam_ref.append("(exp_ph,{})(ext23)".format(self.name))
        return lparam_ref

    def take_param_const(self):
        lparam_const = []
        if (self.scale[2]!=""): lparam_const.append("(exp_ph,{})(scale)".format(self.name))
        if (self.igsize[2]!=""): lparam_const.append("(exp_ph,{})(igsize)".format(self.name))
        if (self.ext11[2]!=""): lparam_const.append("(exp_ph,{})(ext11)".format(self.name))
        if (self.ext22[2]!=""): lparam_const.append("(exp_ph,{})(ext22)".format(self.name))
        if (self.ext33[2]!=""): lparam_const.append("(exp_ph,{})(ext33)".format(self.name))
        if (self.ext12[2]!=""): lparam_const.append("(exp_ph,{})(ext12)".format(self.name))
        if (self.ext13[2]!=""): lparam_const.append("(exp_ph,{})(ext13)".format(self.name))
        if (self.ext23[2]!=""): lparam_const.append("(exp_ph,{})(ext23)".format(self.name))
        return lparam_const

    def out_to_string(self,powder,wavelength,lFRexp=None,lsFRexp=None):
        """
        output information about phase in the experiment
        """
        llsout=[]
        lhkl=self.hkl
        lsthovl=[cfunc.calcsthovl(hkl,self.ucp) for hkl in lhkl]

        if (powder):
            llsout.append("Integral intensities for phase {} (scale is {:10.7f}):\n(NUCL2, MSFp2 and cross are in 10-24cm)".format(self.name,self.scale[0]))
            llsout.append("Ig size is \n {:8.5f}\n ".format(self.igsize[0]))
            llsout.append("Extinction is \n {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f}\n ".format(self.ext11[0],self.ext22[0],self.ext33[0],self.ext12[0],self.ext13[0],self.ext23[0]))
            #for hkl,sthovl,mult,FN,MSFp2,cross,IntU,IntD,SFTensor,THETA,mT,fm2_sq,fm2_perp in zip(lhkl,lsthovl,self.multhkl,self.FNhkl,self.MSFpsq,self.cross,self.IntModUPhkl,self.IntModDOWNhkl,self.SFThkl,self.THETA,self.mT,self.FMhkl_sq_aver,self.FMhkl_perp_aver):
            for hkl,sthovl,mult,FN,MSFp2,cross,SFTensor,THETA,mT in zip(lhkl,lsthovl,self.multhkl,self.FNhkl,self.MSFpsq,self.cross,self.SFThkl,self.THETA,self.mT):
                NUCL=abs(FN)**2
                llsout.append("  h  k  l mult   ReFN    NUCL2    MSFp2     cross     sthovl")
                llsout.append("{0[0]:3}{0[1]:3}{0[2]:3}{1:4}{4:8.3f}{6:9.3f}{7:9.3f}{8:10.3f}{2:8.3f}".format(hkl,mult,sthovl,None,FN.real,FN.imag,NUCL,MSFp2,cross))
                #llsout.append("  averaged fm_sq and fm along H for Ising model  IT IS NOT USED, JUST TEST")
                #llsout.append(" {:9.4f}+ j * {:9.4f},  {:9.4f} + j * {:9.4f}".format(fm2_sq.real, fm2_sq.imag, fm2_perp.real, fm2_perp.imag))
                llsout.append("  structure factor tensor in local Cartesian system")
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[0][0].real,SFTensor[0][1].real,SFTensor[0][2].real,SFTensor[0][0].imag,SFTensor[0][1].imag,SFTensor[0][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[1][0].real,SFTensor[1][1].real,SFTensor[1][2].real,SFTensor[1][0].imag,SFTensor[1][1].imag,SFTensor[1][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[2][0].real,SFTensor[2][1].real,SFTensor[2][2].real,SFTensor[2][0].imag,SFTensor[2][1].imag,SFTensor[2][2].imag))
                llsout.append("  structure factor tensor after T transformation (invT * SFT * T)")
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(THETA[0][0].real,THETA[0][1].real,THETA[0][2].real,THETA[0][0].imag,THETA[0][1].imag,THETA[0][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(THETA[1][0].real,THETA[1][1].real,THETA[1][2].real,THETA[1][0].imag,THETA[1][1].imag,THETA[1][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(THETA[2][0].real,THETA[2][1].real,THETA[2][2].real,THETA[2][0].imag,THETA[2][1].imag,THETA[2][2].imag))
                llsout.append("  matrix T")
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}".format(mT[0][0].real,mT[0][1].real,mT[0][2].real))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}".format(mT[1][0].real,mT[1][1].real,mT[1][2].real))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}".format(mT[2][0].real,mT[2][1].real,mT[2][2].real))
                llsout.append(" ")
        else:
            llsout.append("Phase {}:".format(self.name))
            llsout.append("UB matrix is \n {0[0][0]:8.5f} {0[0][1]:8.5f} {0[0][2]:8.5f}\n {0[1][0]:8.5f} {0[1][1]:8.5f} {0[1][2]:8.5f}\n {0[2][0]:8.5f} {0[2][1]:8.5f} {0[2][2]:8.5f}\n ".format(self.mUB))
            llsout.append("Extinction is \n {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f}\n ".format(self.ext11[0],self.ext22[0],self.ext33[0],self.ext12[0],self.ext13[0],self.ext23[0]))
            #llsout.append("U  matrix is \n {0[0][0]:8.5f} {0[0][1]:8.5f} {0[0][2]:8.5f}\n {0[1][0]:8.5f} {0[1][1]:8.5f} {0[1][2]:8.5f}\n {0[2][0]:8.5f} {0[2][1]:8.5f} {0[2][2]:8.5f}\n ".format(self.mU))
            for hkl, sthovl, FRmod_ext, FRmod, FRexp, sFRexp, msfperp, msf, FN, SFTensor, k in zip(self.hkl,lsthovl,self.FRmod_ext,self.FRmod,lFRexp,lsFRexp,self.MSFperp,self.MSFhkl,self.FNhkl,self.SFThkl,self.kglob):
                llsout.append("   h  k  l     FRexp    sFRexp FRmod_ext     ReFN      ImFN     |FMp| 2|FMpZ|  sthovl   FRmod")
                llsout.append(" {0[0]:3}{0[1]:3}{0[2]:3} {1:9.4f} {2:9.4f}{3:9.4f} {4:9.4f} {5:9.4f} {6:9.4f}{7:8.3f}{8:8.3f}{9:8.3f}".format(hkl,FRexp,sFRexp,FRmod_ext,FN.real,FN.imag,sum([abs(hh1)**2 for hh1 in msfperp])**0.5,2*abs(msfperp[2]),sthovl,FRmod))
                llsout.append("  MSF / scat.vec. k / MSF perp. to k in global Cartesian system")
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(msf[0].real,msf[1].real,msf[2].real,msf[0].imag,msf[1].imag,msf[2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}".format(k[0].real,k[1].real,k[2].real))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(msfperp[0].real,msfperp[1].real,msfperp[2].real,msfperp[0].imag,msfperp[1].imag,msfperp[2].imag))
                llsout.append("  structure factor tensor in local Cartesian system")
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[0][0].real,SFTensor[0][1].real,SFTensor[0][2].real,SFTensor[0][0].imag,SFTensor[0][1].imag,SFTensor[0][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[1][0].real,SFTensor[1][1].real,SFTensor[1][2].real,SFTensor[1][0].imag,SFTensor[1][1].imag,SFTensor[1][2].imag))
                llsout.append(" {:9.4f},{:9.4f},{:9.4f}    + j * {:9.4f},{:9.4f},{:9.4f}".format(SFTensor[2][0].real,SFTensor[2][1].real,SFTensor[2][2].real,SFTensor[2][0].imag,SFTensor[2][1].imag,SFTensor[2][2].imag))
                llsout.append(" ")
        return llsout

class ccore_ph(object):
    """
    have information of input parameter for phase
    """
    def __init__(self, cag_ph = None):
        self.name = "Empty"
        self.spgr = "P1"
        self.spgr_code = "1"
        self.cella = [0., False, ""]
        self.cellb = [0., False, ""]
        self.cellc = [0., False, ""]
        self.cellalpha = [0.,False, ""]
        self.cellbeta = [0.,False, ""]
        self.cellgamma = [0.,False, ""]
        self.atom = []
        self.fdirprog = None
        if cag_ph != None:
            self.take_from_agent(cag_ph)
            
    def set_fdirprog(self, fdirprog):
        self.fdirprog = fdirprog

    def take_from_agent(self, cag_ph):
        self.name = cag_ph.name
        self.spgr = cag_ph.spgr
        self.spgr_code = cag_ph.spgr_code
        self.fdirprog = cag_ph.fdirprog
        self.cella = [hh for hh in cag_ph.cella]
        self.cellb = [hh for hh in cag_ph.cellb]
        self.cellc = [hh for hh in cag_ph.cellc]
        self.cellalpha = [hh for hh in cag_ph.cellalpha]
        self.cellbeta = [hh for hh in cag_ph.cellbeta]
        self.cellgamma = [hh for hh in cag_ph.cellgamma]
        latom = []
        for hh in cag_ph.atom:
            c_at = ccore_ph_at()
            c_at.take_from_agent(hh)
            latom.append(c_at)
        self.atom = latom

    def load_exp_data(self):
        """
        preiminary calculations before refinement
        """
        lmessage = []
        self.getsymm()
        self.calcmB()
        dsymm = self.dsymm
        if dsymm == None:
            lmessage = ["Can not find the elements of symmetry for phase '{:}'".format(self.name)]
            return lmessage

        for cat in self.atom:
            coord = [cat.coordx[0], cat.coordy[0], cat.coordz[0]]
            lelat, leluniqat = cfunc.els4pos(dsymm, coord)
            if (cat.type) != "":
                cat.get_j0j2(self.fdirprog)
            cat.elat = lelat
            cat.eluniqat = leluniqat
            ntot = len(dsymm['elsymm'])*len(dsymm['orig'])
            if (dsymm['centr']): ntot *= 2
            cat.mult = int(ntot*1./len(lelat))
            cat.symmconstr(lelat)
        return lmessage

    def calcmB(self):
        """B matrix:
        crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
        This definition corresponds to the paper Busing, Levy Acta Cryst. (1967) 22, 457
        """
        a = self.cella[0]
        b = self.cellb[0]
        c = self.cellc[0]
        alpha = self.cellalpha[0]
        beta  = self.cellbeta[0]
        gamma = self.cellgamma[0]
        ucp=[a,b,c,alpha,beta,gamma]
        rucp=cfunc.calcrucp(ucp)

        [ra,rb,rc,ralpha,rbeta,rgamma]=rucp
        alphar,betar,gammar=math.radians(alpha),math.radians(beta),math.radians(gamma)
        ralphar,rbetar,rgammar=math.radians(ralpha),math.radians(rbeta),math.radians(rgamma)
        ca,cb,cg=math.cos(alphar),math.cos(betar),math.cos(gammar)
        rca,rcb,rcg=math.cos(ralphar),math.cos(rbetar),math.cos(rgammar)
        sa,sb,sg=math.sin(alphar),math.sin(betar),math.sin(gammar)
        rsa,rsb,rsg=math.sin(ralphar),math.sin(rbetar),math.sin(rgammar)
        self.mB=[[ra,  rb*rcg,  rc*rcb],
            [0.,  rb*rsg, -rc*rsb*ca],
            [0.,      0.,  1./c]]
        return

    def get_constr_ucp(self):
        lconstr = []
        ph_name = self.name
        pa = "(ph,{:})(cella)".format(ph_name)
        pb = "(ph,{:})(cellb)".format(ph_name)
        pc = "(ph,{:})(cellc)".format(ph_name)
        palpha = "(ph,{:})(cellalpha)".format(ph_name)
        pbeta = "(ph,{:})(cellbeta)".format(ph_name)
        pgamma = "(ph,{:})(cellgamma)".format(ph_name)

        self.getsymm()
        dsymm = self.dsymm
        no = dsymm["number"]
        n_choise = dsymm["choise"]
        if ((no >= 16)&(no <= 74)):
            #orthorombic
            lconstr.append((palpha, 90.))
            lconstr.append((pbeta, 90.))
            lconstr.append((pgamma, 90.))
        elif ((no >= 75)&(no <= 142)):
            #tetragonal
            lconstr.append((pb, 1., pa))
            lconstr.append((palpha, 90.))
            lconstr.append((pbeta, 90.))
            lconstr.append((pgamma, 90.))
        elif ((no >= 143)&(no <= 167)):
            if n_choise == 1:
                #hexagonal
                lconstr.append((pb, 1., pa))
                lconstr.append((palpha, 90.))
                lconstr.append((pbeta, 90.))
                lconstr.append((pgamma, 120.))
            else:
                #trigonal
                lconstr.append((pb, 1., pa))
                lconstr.append((pc, 1., pa))
                lconstr.append((pbeta, 1., palpha))
                lconstr.append((pgamma, 1., palpha))
        elif ((no >= 168)&(no <= 194)):
            #hexagonal
            lconstr.append((pb, 1., pa))
            lconstr.append((palpha, 90.))
            lconstr.append((pbeta, 90.))
            lconstr.append((pgamma, 120.))
        elif ((no >= 195)&(no <= 230)):
            #cubic
            lconstr.append((pb, 1., pa))
            lconstr.append((pc, 1., pa))
            lconstr.append((palpha, 90.))
            lconstr.append((pbeta, 90.))
            lconstr.append((pgamma, 90.))
        return lconstr

    def get_constr_adp(self):
        lconstr = []
        ph_name = self.name
        print "\nphase: {:}".format(ph_name)
        self.getsymm()
        dsymm = self.dsymm
        lelsymm = dsymm["elsymm"]
        for atom in self.atom:
            atom.symmconstr(lelsymm)
            matrix = atom.condm
            at_name = atom.name
            pbeta11 = "(ph,{:})(atom,{:})(beta11)".format(ph_name, at_name)
            pbeta22 = "(ph,{:})(atom,{:})(beta22)".format(ph_name, at_name)
            pbeta33 = "(ph,{:})(atom,{:})(beta33)".format(ph_name, at_name)
            pbeta12 = "(ph,{:})(atom,{:})(beta12)".format(ph_name, at_name)
            pbeta13 = "(ph,{:})(atom,{:})(beta13)".format(ph_name, at_name)
            pbeta23 = "(ph,{:})(atom,{:})(beta23)".format(ph_name, at_name)

            pchi11 = "(ph,{:})(atom,{:})(chi11)".format(ph_name, at_name)
            pchi22 = "(ph,{:})(atom,{:})(chi22)".format(ph_name, at_name)
            pchi33 = "(ph,{:})(atom,{:})(chi33)".format(ph_name, at_name)
            pchi12 = "(ph,{:})(atom,{:})(chi12)".format(ph_name, at_name)
            pchi13 = "(ph,{:})(atom,{:})(chi13)".format(ph_name, at_name)
            pchi23 = "(ph,{:})(atom,{:})(chi23)".format(ph_name, at_name)
            #it is totaly bad solution but enough simple !!!!!!!!!!!!!!!!!!!!!!
            if (1, -1, 0, 0, 0, 0) in matrix:
                lconstr.append((pbeta22, 1., pbeta11))
                lconstr.append((pchi22, 1., pchi11))
            if (1, 0, -1, 0, 0, 0) in matrix:
                lconstr.append((pbeta33, 1., pbeta11))
                lconstr.append((pchi33, 1., pchi11))

            if ( 1, 0, 0, 2, 0, 0) in matrix:
                lconstr.append((pbeta12, 0.5, pbeta11))
                lconstr.append((pchi12, 0.5, pchi11))
            if ( 1, 0, 0, 0, 2, 0) in matrix:
                lconstr.append((pbeta13, 0.5, pbeta11))
                lconstr.append((pchi13, 0.5, pchi11))
            if ( 1, 0, 0, 0, 0, 2) in matrix:
                lconstr.append((pbeta23, 0.5, pbeta11))
                lconstr.append((pchi23, 0.5, pchi11))

            if ( 0, 0, 0, 1, 0, 0) in matrix:
                lconstr.append((pbeta12, 0.))
                lconstr.append((pchi12, 0.))
            else:
                if ( 0, 0, 0, 1, -1, 0) in matrix:
                    lconstr.append((pbeta13, 1., pbeta12))
                    lconstr.append((pchi13, 1., pchi12))
                if ( 0, 0, 0, 1, 0, -1) in matrix:
                    lconstr.append((pbeta23, 1., pbeta12))
                    lconstr.append((pchi23, 1., pchi12))

            if ( 0, 0, 0, 0, 1, 0) in matrix:
                lconstr.append((pbeta13, 0.))
                lconstr.append((pchi13, 0.))

            if ( 0, 0, 0, 0, 0, 1) in matrix:
                lconstr.append((pbeta23, 0.))
                lconstr.append((pchi23, 0.))
        return lconstr

    def getsymm(self):
        """get symmetry from space group
        fdir is the folder which contents file 'symmetry.xml'

        spgr = "227 2" (numbder 2227 choice 2) or "227" (numer 227 choice 1) or "Fd-3m 2" or "Fd-3m"
        """
        fdir = self.fdirprog
        number_or_name_space_origin = "".join(self.spgr.split())
        n_choise = "{:}".format(self.spgr_code)
        if isinstance(self.spgr_code, float):
            n_choise = "{:}".format(int(self.spgr_code))
        

        if number_or_name_space_origin.isdigit():
            spgr_n = number_or_name_space_origin
            spgr_name = ""
        else:
            spgr_n = ""
            spgr_name = number_or_name_space_origin

        fitables=os.path.join(fdir,"itables.txt")
        ldcard = cfunc.read_el_cards(fitables)

        flag = False
        for dcard in ldcard:
            if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choice"][0] == n_choise)):
                flag = True
                break
        if (not flag):
            self.dsymm=None
            return
        lelsymm = []
        for ssymm in dcard["symmetry"]:
            lelsymm.append(cfunc.tr2elsymm(ssymm))
        centr = dcard["centr"][0]=="true"
        pcentr = [float(hh) for hh in dcard["pcentr"][0].split(",")]
        fletter = dcard["name"][0]
        spgr = dcard["name"]
        number = int(dcard["number"])
        if ((fletter == "P")|(fletter == "R")):
            lorig = [(0, 0, 0)]
        elif fletter == "C":
            lorig = [(0, 0, 0), (0.5, 0.5, 0)]
        elif fletter == "I":
            lorig = [(0, 0, 0), (0.5, 0.5, 0.5)]
        elif fletter == "F":
            lorig = [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)]
        else:
            print "Undefined syngony"
        dsymm={"centr":centr, "elsymm":lelsymm, "orig":lorig, "pcentr":pcentr,
               "spgr":spgr, "number":number, "choise": int(n_choise)}

        self.dsymm=dsymm

    def calcFNhkl(self,lhkl):
        ucp=[self.cella[0],self.cellb[0],self.cellc[0],
             self.cellalpha[0],self.cellbeta[0],self.cellgamma[0],]
        lsthovl=[cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]
        dsymm=self.dsymm
        lelsymm,lorig,centr=dsymm["elsymm"],dsymm["orig"],dsymm["centr"]
        lFhelp=[0]
        lFhkl = []
        for hkl,sthovl in zip(lhkl,lsthovl):
            for cat in self.atom:
                occ,mult=cat.occ[0],cat.mult
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                bscat=cat.bscat[0]
                AA=mult*occ*bscat
                if (AA!=0.):
                    for isymm,symm in enumerate(lelsymm):
                        hkls=[hkl[0]*symm[1]+hkl[1]*symm[5]+hkl[2]*symm[9], hkl[0]*symm[2]+hkl[1]*symm[6]+hkl[2]*symm[10], hkl[0]*symm[3]+hkl[1]*symm[7]+hkl[2]*symm[11]]
                        #if (cat.modedw=='iso'):
                        biso = cat.biso[0]
                        DWF_i = math.exp(-biso*sthovl**2)
                        #elif (cat.modedw=='aniso'):
                        #MISTAKE: introduce constraints on the betas
                        beta11,beta22=cat.beta11[0],cat.beta22[0]
                        beta33,beta12=cat.beta33[0],cat.beta12[0]
                        beta13,beta23=cat.beta13[0],cat.beta23[0]
                        #check it, not sure about 2 in front of beta12, beta13 and beta 23
                        #it is not correct, introduce equivalent reflections
                        DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2.*beta12*hkls[0]*hkls[1]+2.*beta13*hkls[0]*hkls[2]+2.*beta23*hkls[1]*hkls[2]))
                        #else:
                            #print 'Unknown thermal parameters for atom {0[name]:}'.format(datom)
                        #    print 'Program is stopped'
                        #    quit()
                        chelp = hkls[0]*coord[0]+hkl[0]*symm[0]+hkls[1]*coord[1]+hkl[1]*symm[4]+hkls[2]*coord[2]+hkl[2]*symm[8]
                        lFhelp.append(AA * DWF_i * DWF_a * cmath.exp(2*math.pi*1j*chelp))
            Fhklas=sum(lFhelp)*1./len(lelsymm)
            lFhelp=[]
            for orig in lorig:
                lFhelp.append(Fhklas*cmath.exp( 2*math.pi*1j* (hkl[0]*orig[0]+hkl[1]*orig[1]+hkl[2]*orig[2])))
            Fhklas=sum(lFhelp)*1./len(lFhelp)
            lFhelp=[]
            if (centr):
                orig=dsymm["pcentr"]
                Fhkl=(Fhklas+Fhklas.conjugate()*cmath.exp(2*2*math.pi*1j* (hkl[0]*orig[0]+hkl[1]*orig[1]+hkl[2]*orig[2])))*0.5
            else:
                Fhkl=Fhklas
            lFhkl.append(Fhkl)
        return lFhkl


    def calcFMhkl_ising_powder(self, lhkl, lmTT):
        """
        calculate the <M_{||, H}**2> and <abs(M_{||, H}*h)> along magnetic field for each atom for powder in case of Ising model
        the orientation of magnetic moment along the Ising axis is given by vector components mag_ia, mag_ib, mag_ic
        lmTT -- transposed (inversed) T matrix
        """
        def calc_av_acosbsin(a,b):
            """
            calculate integral(abs(a*cos(phi)+b*sin(spi)), (phi, 0, 2*pi)) /2*pi
            """
            phi_1 = math.atan2(-a, b)
            return ((2./math.pi)*abs(-a*math.sin(phi_1)+b*math.cos(phi_1)))

        def calc_aibiajbj(ai, bi, aj, bj, phi_1, phi_2):
            """
            calculate 1/2pi*integral((ai*cos(phi)+bi*sin(phi))*(aj*cos(phi)+bj*sin(phi)),(phi, phi_1, phi_2))
            """
            cval = 0.5*(math.sin(2.*phi_2)-math.sin(2.*phi_1))
            cos_phi_sq = ((phi_2-phi_1)+cval)/(4.*math.pi)
            sin_phi_sq = ((phi_2-phi_1)-cval)/(4.*math.pi)
            cos_sin_phi = ((math.sin(phi_2))**2-(math.sin(phi_1))**2)/(4.*math.pi)
            res = ai*aj*cos_phi_sq + (ai*bj+aj*bi)*cos_sin_phi + bi*bj*sin_phi_sq
            return res

        def calc_ab_aibiajbj(ai, bi, aj, bj):
            """
            calculate 1/2pi*integral(abs((ai*cos(phi)+bi*sin(phi))*(aj*cos(phi)+bj*sin(phi))),(phi, 0, 2*pi))
            """
            phi_1 = math.atan2(-ai, bi)%(math.pi)
            phi_2 = math.atan2(-aj, bj)%(math.pi)
            if phi_1 <= phi_2:
                phi_A, phi_B = phi_1, phi_2
            else:
                phi_A, phi_B = phi_2, phi_1
            res_1 = abs(calc_aibiajbj(ai, bi, aj, bj, phi_A, phi_B))
            res_2 = abs(calc_aibiajbj(ai, bi, aj, bj, phi_B, phi_A+math.pi))
            return 2.*(res_1+res_2)

        lFMhkl_sq_aver = []
        lFMhkl_perp_aver =  [0. for hkl in lhkl]
        ucp = [self.cella[0], self.cellb[0], self.cellc[0],
               self.cellalpha[0], self.cellbeta[0], self.cellgamma[0]]
        lsthovl = [cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]

        llhelp_param =[[] for hkl in lhkl]
        for cat in self.atom:
            if (cat.modemagn):
                lmff = [cat.calcMFFSpherMod(sthovl,kappa=cat.kappa[0]) for sthovl in lsthovl]
                occ = cat.occ[0]
                #chi in 10-12 cm; chim in muB (it is why here 0.2695)
                #it is not correct
                #it should be some over elements symmetry which transfer atom from one site to another
                leluniqat = cat.eluniqat

                #lchiLOCuniqat=cat.chiLOCuniqat
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                lcoordeq=[((elseq[0]+elseq[1]*coord[0]+elseq[ 2]*coord[1]+elseq[ 3]*coord[2])%1.,
                       (elseq[4]+elseq[5]*coord[0]+elseq[ 6]*coord[1]+elseq[ 7]*coord[2])%1.,
                       (elseq[8]+elseq[9]*coord[0]+elseq[10]*coord[1]+elseq[11]*coord[2])%1.) for elseq in leluniqat]
                mag = [cat.mag_ia[0],cat.mag_ib[0],cat.mag_ic[0]]
                #condM=cat.condm
                #for hh in condM:
                #    hh3=abs(sum([hh2*chim for hh2,chim in zip(hh,chi)]))
                #    if (hh3>0.0001):
                #        print "Mistake with constraints"
                #chic = chi
                lmagr=[cfunc.magrot(mag,elseq) for elseq in leluniqat]#rotation of chi by element of symmetry
                for elseq,coordeq,magr in zip(leluniqat,lcoordeq,lmagr):
                    #magrm in 10-12 cm; magr in muB (it is why here 0.2695)
                    magrm=[0.2695*hh2 for hh2 in magr]
                    mag_loc_1 = cfunc.magLOC(magrm,ucp)#representation of chi in Crystallographical coordinate system defined as x||a*, z||c, y= [z x] (right handed)
                    lmag_loc = [cfunc.multMv(mTT, mag_loc_1) for mTT in lmTT] #representation of chi in new Crystallographical coordinate system defined by T matrix
                    lansw_perp = []
                    for hkl, sthovl, mff, lhelp_param, FMhkl_perp_aver, mag_loc in zip(lhkl, lsthovl, lmff, llhelp_param, lFMhkl_perp_aver, lmag_loc):
                        phase=sum([2*math.pi*1j*hh1*hh2 for hh1,hh2 in zip(hkl,coordeq)])
                        #if (cat.modedw=='iso'):
                        biso = cat.biso[0]
                        DWF_i = math.exp(-biso*sthovl**2)
                        #elif (cat.modedw=='aniso'):
                        #    #not sure about hkls
                        #    #not sure about 2 in front of beta12, beta13, beta23
                        beta11,beta22=cat.beta11[0],cat.beta22[0]
                        beta33,beta12=cat.beta33[0],cat.beta12[0]
                        beta13,beta23=cat.beta13[0],cat.beta23[0]
                        hkls=[hkl[0]*elseq[1]+hkl[1]*elseq[5]+hkl[2]*elseq[9], hkl[0]*elseq[2]+hkl[1]*elseq[6]+hkl[2]*elseq[10], hkl[0]*elseq[3]+hkl[1]*elseq[7]+hkl[2]*elseq[11]]
                        DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2*beta12*hkls[0]*hkls[1]+2*beta13*hkls[0]*hkls[2]+2*beta23*hkls[1]*hkls[2]))
                        valhelp = occ*DWF_i*DWF_a*mff*cmath.exp(phase)
                        param_to_save = (mag_loc[0], mag_loc[1], valhelp)
                        val_perp = calc_av_acosbsin(mag_loc[0], mag_loc[1])*valhelp
                        lhelp_param.append(param_to_save)
                        lansw_perp.append(FMhkl_perp_aver+val_perp)
                    lFMhkl_perp_aver = lansw_perp# can be complex
        for lhelp_param in llhelp_param:
            lval_sq = []
            for param_1 in lhelp_param:
                for param_2 in lhelp_param:
                    val_sq = calc_ab_aibiajbj(param_1[0], param_1[1], param_2[0], param_2[1])*param_1[2]*param_2[2].conjugate()
                    lval_sq.append(val_sq)
            lFMhkl_sq_aver.append(sum(lval_sq))# should be real
        return lFMhkl_sq_aver, lFMhkl_perp_aver


    def calcFMhkl_ising_mono(self, lhkl, h_loc):
        """
        calculate magnetic structure factor for given magnetic moment for each atom (magneti momen is not function of applied magnetic field) in frame of the ising model

        the orientation of magnetic moments are along the applied magnetic fields
        if the magnetic moment is perpendicular to the applied magnetic field then it is absent
        """
        res = [0.,0.,0.]
        lres = [res for hkl in lhkl]
        ucp = [self.cella[0],self.cellb[0],self.cellc[0],
               self.cellalpha[0],self.cellbeta[0],self.cellgamma[0],]
        lsthovl = [cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]
        for cat in self.atom:
            if (cat.modemagn):
                lmff = [cat.calcMFFSpherMod(sthovl,kappa=cat.kappa[0]) for sthovl in lsthovl]
                occ = cat.occ[0]
                #chi in 10-12 cm; chim in muB (it is why here 0.2695)
                #it is not correct
                #it should be some over elements symmetry which transfer atom from one site to another
                leluniqat = cat.eluniqat

                #lchiLOCuniqat=cat.chiLOCuniqat
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                lcoordeq=[((elseq[0]+elseq[1]*coord[0]+elseq[ 2]*coord[1]+elseq[ 3]*coord[2])%1.,
                       (elseq[4]+elseq[5]*coord[0]+elseq[ 6]*coord[1]+elseq[ 7]*coord[2])%1.,
                       (elseq[8]+elseq[9]*coord[0]+elseq[10]*coord[1]+elseq[11]*coord[2])%1.) for elseq in leluniqat]
                mag = [cat.mag_ia[0], cat.mag_ib[0], cat.mag_ic[0]]
                tol  = cat.tol_ising[0]
                #condM=cat.condm
                #for hh in condM:
                #    hh3=abs(sum([hh2*chim for hh2,chim in zip(hh,chi)]))
                #    if (hh3>0.0001):
                #        print "Mistake with constraints"
                #chic = chi
                lmagr=[cfunc.magrot(mag, elseq) for elseq in leluniqat]#rotation of chi by element of symmetry
                for elseq,coordeq,magr in zip(leluniqat,lcoordeq,lmagr):
                    #magrm in 10-12 cm; magr in muB (it is why here 0.2695)
                    magrm = [0.2695*hh2 for hh2 in magr]
                    mag_loc = cfunc.magLOC(magrm, ucp)#representation of chi in Crystallographical coordinate system defined as x||a*, z||c, y= [z x] (right handed)
                    mult_m_h = sum([hh1*hh2 for hh1, hh2 in zip(h_loc, mag_loc)])
                    if abs(mult_m_h) > tol:
                        sign_mult_m_h = float(2*(mult_m_h > 0.)-1)
                        lansw=[]
                        for hkl,sthovl,mff,res in zip(lhkl,lsthovl,lmff,lres):
                            phase=sum([2*math.pi*1j*hh1*hh2 for hh1,hh2 in zip(hkl,coordeq)])
                            #if (cat.modedw=='iso'):
                            biso = cat.biso[0]
                            DWF_i = math.exp(-biso*sthovl**2)
                            #elif (cat.modedw=='aniso'):
                            #    #not sure about hkls
                            #    #not sure about 2 in front of beta12, beta13, beta23
                            beta11,beta22=cat.beta11[0],cat.beta22[0]
                            beta33,beta12=cat.beta33[0],cat.beta12[0]
                            beta13,beta23=cat.beta13[0],cat.beta23[0]
                            hkls=[hkl[0]*elseq[1]+hkl[1]*elseq[5]+hkl[2]*elseq[9], hkl[0]*elseq[2]+hkl[1]*elseq[6]+hkl[2]*elseq[10], hkl[0]*elseq[3]+hkl[1]*elseq[7]+hkl[2]*elseq[11]]
                            DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2*beta12*hkls[0]*hkls[1]+2*beta13*hkls[0]*hkls[2]+2*beta23*hkls[1]*hkls[2]))
                            valhelp=occ*DWF_i*DWF_a*mff*cmath.exp(phase)
                            val=[sign_mult_m_h*valhelp*hh2  for hh2 in mag_loc]
                            lansw.append([hh3+hh4 for hh3,hh4 in zip(res,val)])
                        lres = lansw
        return lres


    def calcFMhkl(self, lhkl):
        """
        calculate magnetic structure factor for given magnetic moment for each atom (magneti momen is not function of applied magnetic field)
        """
        res = [0.,0.,0.]
        lres = [res for hkl in lhkl]
        ucp = [self.cella[0],self.cellb[0],self.cellc[0],
               self.cellalpha[0],self.cellbeta[0],self.cellgamma[0],]
        lsthovl = [cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]
        for cat in self.atom:
            if (cat.modemagn):
                lmff = [cat.calcMFFSpherMod(sthovl,kappa=cat.kappa[0]) for sthovl in lsthovl]
                occ = cat.occ[0]
                #chi in 10-12 cm; chim in muB (it is why here 0.2695)
                #it is not correct
                #it should be some over elements symmetry which transfer atom from one site to another
                leluniqat = cat.eluniqat

                #lchiLOCuniqat=cat.chiLOCuniqat
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                lcoordeq=[((elseq[0]+elseq[1]*coord[0]+elseq[ 2]*coord[1]+elseq[ 3]*coord[2])%1.,
                       (elseq[4]+elseq[5]*coord[0]+elseq[ 6]*coord[1]+elseq[ 7]*coord[2])%1.,
                       (elseq[8]+elseq[9]*coord[0]+elseq[10]*coord[1]+elseq[11]*coord[2])%1.) for elseq in leluniqat]
                mag = [cat.mag_ia[0],cat.mag_ib[0],cat.mag_ic[0]]
                #condM=cat.condm
                #for hh in condM:
                #    hh3=abs(sum([hh2*chim for hh2,chim in zip(hh,chi)]))
                #    if (hh3>0.0001):
                #        print "Mistake with constraints"
                #chic = chi
                lmagr=[cfunc.magrot(mag,elseq) for elseq in leluniqat]#rotation of chi by element of symmetry
                for elseq,coordeq,magr in zip(leluniqat,lcoordeq,lmagr):
                    #magrm in 10-12 cm; magr in muB (it is why here 0.2695)
                    magrm=[0.2695*hh2 for hh2 in magr]
                    mag_loc=cfunc.magLOC(magrm,ucp)#representation of chi in Crystallographical coordinate system defined as x||a*, z||c, y= [z x] (right handed)
                    lansw=[]
                    for hkl,sthovl,mff,res in zip(lhkl,lsthovl,lmff,lres):
                        phase=sum([2*math.pi*1j*hh1*hh2 for hh1,hh2 in zip(hkl,coordeq)])
                        #if (cat.modedw=='iso'):
                        biso = cat.biso[0]
                        DWF_i = math.exp(-biso*sthovl**2)
                        #elif (cat.modedw=='aniso'):
                        #    #not sure about hkls
                        #    #not sure about 2 in front of beta12, beta13, beta23
                        beta11,beta22=cat.beta11[0],cat.beta22[0]
                        beta33,beta12=cat.beta33[0],cat.beta12[0]
                        beta13,beta23=cat.beta13[0],cat.beta23[0]
                        hkls=[hkl[0]*elseq[1]+hkl[1]*elseq[5]+hkl[2]*elseq[9], hkl[0]*elseq[2]+hkl[1]*elseq[6]+hkl[2]*elseq[10], hkl[0]*elseq[3]+hkl[1]*elseq[7]+hkl[2]*elseq[11]]
                        DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2*beta12*hkls[0]*hkls[1]+2*beta13*hkls[0]*hkls[2]+2*beta23*hkls[1]*hkls[2]))
                        valhelp=occ*DWF_i*DWF_a*mff*cmath.exp(phase)
                        val=[valhelp*hh2  for hh2 in mag_loc]
                        lansw.append([hh3+hh4 for hh3,hh4 in zip(res,val)])
                    lres = lansw
        return lres


    def calcSFThkl(self, lhkl):
        res = [[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]
        lres = [res for hkl in lhkl]
        ucp = [self.cella[0],self.cellb[0],self.cellc[0],
               self.cellalpha[0],self.cellbeta[0],self.cellgamma[0],]
        lsthovl = [cfunc.calcsthovl(hkl,ucp) for hkl in lhkl]
        for cat in self.atom:
            if (cat.modemagn):
                lmff = [cat.calcMFFSpherMod(sthovl,kappa=cat.kappa[0]) for sthovl in lsthovl]
                occ = cat.occ[0]
                #chi in 10-12 cm; chim in muB (it is why here 0.2695)
                #it is not correct
                #it should be some over elements symmetry which transfer atom from one site to another
                leluniqat = cat.eluniqat
                #lchiLOCuniqat=cat.chiLOCuniqat
                coord=[cat.coordx[0],cat.coordy[0],cat.coordz[0]]
                lcoordeq=[((elseq[0]+elseq[1]*coord[0]+elseq[ 2]*coord[1]+elseq[ 3]*coord[2])%1.,
                       (elseq[4]+elseq[5]*coord[0]+elseq[ 6]*coord[1]+elseq[ 7]*coord[2])%1.,
                       (elseq[8]+elseq[9]*coord[0]+elseq[10]*coord[1]+elseq[11]*coord[2])%1.) for elseq in leluniqat]
                chi=[cat.chi11[0],cat.chi22[0],cat.chi33[0],cat.chi12[0],cat.chi13[0],cat.chi23[0]]
                condM=cat.condm
                for hh in condM:
                    hh3=abs(sum([hh2*chim for hh2,chim in zip(hh,chi)]))
                    if (hh3>0.0001):
                        print "Mistake with constraints"
                chic = chi
                lchir=[cfunc.chirot(chic,elseq) for elseq in leluniqat]#rotation of chi by element of symmetry
                for elseq,coordeq,chir in zip(leluniqat,lcoordeq,lchir):
                    #chi in 10-12 cm; chim in muB (it is why here 0.2695)
                    chirm=[[0.2695*hh2 for hh2 in hh1] for hh1 in chir]
                    chi_loc=cfunc.chiLOC(chirm,ucp)#representation of chi in Crystallographical coordinate system defined as x||a*, z||c, y= [z x] (right handed)
                    lansw=[]
                    for hkl,sthovl,mff,res in zip(lhkl,lsthovl,lmff,lres):
                        phase=sum([2*math.pi*1j*hh1*hh2 for hh1,hh2 in zip(hkl,coordeq)])
                        #if (cat.modedw=='iso'):
                        biso=cat.biso[0]
                        DWF_i = math.exp(-biso*sthovl**2)
                        #elif (cat.modedw=='aniso'):
                        #not sure about hkls
                        #not sure about 2 in front of beta12, beta13, beta23
                        beta11,beta22 = cat.beta11[0],cat.beta22[0]
                        beta33,beta12 = cat.beta33[0],cat.beta12[0]
                        beta13,beta23 = cat.beta13[0],cat.beta23[0]
                        hkls = [hkl[0]*elseq[1]+hkl[1]*elseq[5]+hkl[2]*elseq[9], hkl[0]*elseq[2]+hkl[1]*elseq[6]+hkl[2]*elseq[10], hkl[0]*elseq[3]+hkl[1]*elseq[7]+hkl[2]*elseq[11]]
                        DWF_a = math.exp(-1*(beta11*hkls[0]**2+beta22*hkls[1]**2+beta33*hkls[2]**2+
                                             2*beta12*hkls[0]*hkls[1]+2*beta13*hkls[0]*hkls[2]+2*beta23*hkls[1]*hkls[2]))
                        valhelp=occ*DWF_i*DWF_a*mff*cmath.exp(phase)
                        val=[[valhelp*hh2  for hh2 in hh1] for hh1 in chi_loc]
                        lansw.append([[hh3+hh4 for hh3,hh4 in zip(hh1,hh2)] for hh1,hh2 in zip(res,val)])
                    lres = lansw
        return lres

    def set_val_by_link(self,slink,val):
        """slink = (U)
           slink = (exp_ph,Tb)(scale)
        """
        val_new, lmessages = set_val_by_link(self, slink, val)
        return val_new, lmessages

    def take_param_ref(self):
        lparam_ref = []
        if (self.cella[1]): lparam_ref.append("(ph,{})(cella)".format(self.name))
        if (self.cellb[1]): lparam_ref.append("(ph,{})(cellb)".format(self.name))
        if (self.cellc[1]): lparam_ref.append("(ph,{})(cellc)".format(self.name))
        if (self.cellalpha[1]): lparam_ref.append("(ph,{})(cellalpha)".format(self.name))
        if (self.cellbeta[1]): lparam_ref.append("(ph,{})(cellbeta)".format(self.name))
        if (self.cellgamma[1]): lparam_ref.append("(ph,{})(cellgamma)".format(self.name))
        for cat in self.atom:
            lparam_ref_at=cat.take_param_ref()
            if (lparam_ref_at!=[]):
                lparam_ref.extend(["(ph,{}){}".format(self.name,hh) for hh in lparam_ref_at])
        return lparam_ref

    def take_param_const(self):
        lparam_const = []
        if (self.cella[2]!=""): lparam_const.append("(ph,{})(cella)".format(self.name))
        if (self.cellb[2]!=""): lparam_const.append("(ph,{})(cellb)".format(self.name))
        if (self.cellc[2]!=""): lparam_const.append("(ph,{})(cellc)".format(self.name))
        if (self.cellalpha[2]!=""): lparam_const.append("(ph,{})(cellalpha)".format(self.name))
        if (self.cellbeta[2]!=""): lparam_const.append("(ph,{})(cellbeta)".format(self.name))
        if (self.cellgamma[2]!=""): lparam_const.append("(ph,{})(cellgamma)".format(self.name))
        for cAt in self.atom:
            lparam_const_at=cAt.take_param_const()
            if (lparam_const_at!=[]):
                lparam_const.extend(["(ph,{}){}".format(self.name,hh) for hh in lparam_const_at])
        return lparam_const

    def out_to_string(self):
        #print output information
        llsout=[]
        llsout.append(70*"*")
        llsout.append("\nPhase {}:".format(self.name))
        ucp=[self.cella[0],self.cellb[0],self.cellc[0],self.cellalpha[0],self.cellbeta[0],self.cellgamma[0]]
        llsout.append(" unit cell: \n  {0[0]:10.5f} {0[1]:10.5f} {0[2]:10.5f} \n     {0[3]:7.2f}    {0[4]:7.2f}    {0[5]:7.2f}".format(ucp))
        llsout.append(" B matrix is \n {0[0][0]:8.5f} {0[0][1]:8.5f} {0[0][2]:8.5f}\n {0[1][0]:8.5f} {0[1][1]:8.5f} {0[1][2]:8.5f}\n {0[2][0]:8.5f} {0[2][1]:8.5f} {0[2][2]:8.5f}\n ".format(self.mB))
        llsout.append("symmetry: ")
        dsymm=self.dsymm
        if (dsymm['centr']):
            llsout.append(" inversion center in position {0[0]:5.3f} {0[1]:5.3f} {0[2]:5.3f}".format(dsymm['pcentr']))
        else:
            llsout.append(" no center of inversion")
        llsout.append(" elements of symmetry:")
        lhelp=[[] for hh1 in range(4)]
        for ihh1,elsymm in enumerate(dsymm['elsymm']):
            lhelp[0].append("      B    Rot.")
            lhelp[1].append("  {0[0]:5.3f} {0[1]:2}{0[2]:2}{0[3]:2};".format(elsymm))
            lhelp[2].append("  {0[4]:5.3f} {0[5]:2}{0[6]:2}{0[7]:2};".format(elsymm))
            lhelp[3].append("  {0[8]:5.3f} {0[9]:2}{0[10]:2}{0[11]:2};".format(elsymm))
            if (((ihh1+1)%5)==0):
                llsout.append("".join(lhelp[0]))
                llsout.append("".join(lhelp[1]))
                llsout.append("".join(lhelp[2]))
                llsout.append("".join(lhelp[3]))
                lhelp=[[] for hh1 in range(4)]
        if (lhelp[0]!=[]):
            llsout.append("".join(lhelp[0]))
            llsout.append("".join(lhelp[1]))
            llsout.append("".join(lhelp[2]))
            llsout.append("".join(lhelp[3]))
        llsout.append("\n shifting:")
        llsout.append("".join([" +{0[0]:5.3f}  ".format(orig) for orig in dsymm['orig']]))
        llsout.append("".join([" +{0[1]:5.3f}  ".format(orig) for orig in dsymm['orig']]))
        llsout.append("".join([" +{0[2]:5.3f}  ".format(orig) for orig in dsymm['orig']]))

        for cat in self.atom:
            llsout.append(70*"*")
            llsout.extend(cat.out_to_string(ucp))
        llsout.append("\n")
        return llsout


class ccore_ph_at(object):
    """
    have information of input parameter for atom in the phase
    """
    def __init__(self):
        self.name = 'Empty'
        self.type = 'O'
        self.type_n = 'O'
        self.modedw = 'iso'
        self.modechi = 'aniso'
        self.modemagn = False
        self.coordx = [0.0000, False, ""]
        self.coordy = [0.0000, False, ""]
        self.coordz = [0.0000, False, ""]
        self.bscat = [0.0000, False, ""]
        self.occ = [1.0000, False, ""]
        self.lfactor = [2.0000, False, ""]
        self.coeff0 = [0., 0., 0., 0., 0., 0., 0.]
        self.coeff2 = [0., 0., 0., 0., 0., 0., 0.]
        self.biso = [0.0000, False, ""]
        self.beta11, self.beta22, self.beta33 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.beta12, self.beta13, self.beta23 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.kappa = [1.0000, False, ""]
        self.chiiso = [0.0000, False, ""]
        self.chi11, self.chi22, self.chi33 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.chi12, self.chi13, self.chi23 = [0., False, ""], [0., False, ""], [0., False, ""]
        cval = 0.0/(3**0.5)
        self.mag_ia, self.mag_ib, self.mag_ic = [cval, False, ""], [cval, False, ""], [cval, False, ""]
        self.al_field, self.tol_ising = [1., False, ""], [0.2, False, ""]
        self.eluniqat = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.coordeq = [0., 0., 0.]
        self.chir = [0., 0., 0., 0., 0., 0.]
        self.condm = []
    def take_from_agent(self, cag_ph):
        self.name = cag_ph.name
        self.type = cag_ph.type
        self.type_n = cag_ph.type_n
        self.modedw = cag_ph.modedw
        self.modechi = cag_ph.modechi
        self.modemagn = cag_ph.modemagn
        self.coordx = [hh for hh in cag_ph.coordx]
        self.coordy = [hh for hh in cag_ph.coordy]
        self.coordz = [hh for hh in cag_ph.coordz]
        self.bscat = [hh for hh in cag_ph.bscat]
        self.occ = [hh for hh in cag_ph.occ]
        self.lfactor = [hh for hh in cag_ph.lfactor]
        self.biso = [hh for hh in cag_ph.biso]
        self.beta11, self.beta22, self.beta33 = [hh for hh in cag_ph.beta11], [hh for hh in cag_ph.beta22], [hh for hh in cag_ph.beta33]
        self.beta12, self.beta13, self.beta23 = [hh for hh in cag_ph.beta12], [hh for hh in cag_ph.beta13], [hh for hh in cag_ph.beta23]
        self.kappa = [hh for hh in cag_ph.kappa]
        self.chiiso = [hh for hh in cag_ph.chiiso]
        self.chi11, self.chi22, self.chi33 = [hh for hh in cag_ph.chi11], [hh for hh in cag_ph.chi22], [hh for hh in cag_ph.chi33]
        self.chi12, self.chi13, self.chi23 = [hh for hh in cag_ph.chi12], [hh for hh in cag_ph.chi13], [hh for hh in cag_ph.chi23]
        self.mag_ia, self.mag_ib, self.mag_ic = [hh for hh in cag_ph.mag_ia], [hh for hh in cag_ph.mag_ib], [hh for hh in cag_ph.mag_ic]
        self.al_field, self.tol_ising = [hh for hh in cag_ph.al_field], [hh for hh in cag_ph.tol_ising]


    def calcMFFSpherMod(self, sthovl,kappa = 1.):
        """Calculate magnetic form factor in frame of Spherical model (Int.Tabl.C.p.592)\n
        LFactor is Lande factor\n
        coeff0 is a list [A,a,B,b,C,c,D] at n=0\n
        coeff2 is a list [A,a,B,b,C,c,D] at n=2\n
        lsthovl is list sin(theta)/lambda in Angstrems**-1

        Calculation of magnetic form factor <j0>,<j2>,<j4>,<j6>\n
        coeff is a list [A,a,B,b,C,c,D] at n=0,2,4,6
        For help see International Table Vol.C p.460
        """
        #not sure about kappa, it is here just for test, by default it is 1.0
        sthovl2 = (sthovl*1./kappa)**2
        [A, a, B, b, C, c, D] = self.coeff0
        res0 = (A*math.exp(-a*sthovl2)+B*math.exp(-b*sthovl2)+C*math.exp(-c*sthovl2)+D)
        [A, a, B, b, C, c, D] = self.coeff2
        res2 = (A*math.exp(-a*sthovl2)+B*math.exp(-b*sthovl2)+C*math.exp(-c*sthovl2)+D)*sthovl2
        res = res0+(1.0-2.0/self.lfactor[0])*res2
        return res

    def set_val_by_link(self, slink, val):
        """
        set value by link
        """
        val_new, lmessages = set_val_by_link(self, slink, val)
        return val_new, lmessages

    def take_param_ref(self):
        """
        take refinement parameters from model
        """
        lparam_ref = []
        if self.coordx[1]:
            lparam_ref.append("(atom,{})(coordx)".format(self.name))
        if self.coordy[1]:
            lparam_ref.append("(atom,{})(coordy)".format(self.name))
        if self.coordz[1]:
            lparam_ref.append("(atom,{})(coordz)".format(self.name))
        if self.bscat[1]:
            lparam_ref.append("(atom,{})(bscat)".format(self.name))
        if self.occ[1]:
            lparam_ref.append("(atom,{})(occ)".format(self.name))
        #if self.modedw == "aniso":
        if self.beta11[1]:
            lparam_ref.append("(atom,{})(beta11)".format(self.name))
        if self.beta22[1]:
            lparam_ref.append("(atom,{})(beta22)".format(self.name))
        if self.beta33[1]:
            lparam_ref.append("(atom,{})(beta33)".format(self.name))
        if self.beta12[1]:
            lparam_ref.append("(atom,{})(beta12)".format(self.name))
        if self.beta13[1]:
            lparam_ref.append("(atom,{})(beta13)".format(self.name))
        if self.beta23[1]:
            lparam_ref.append("(atom,{})(beta23)".format(self.name))
        #else:
        if self.biso[1]:
            lparam_ref.append("(atom,{})(biso)".format(self.name))
        if self.modemagn:
            if self.lfactor[1]:
                lparam_ref.append("(atom,{})(lfactor)".format(self.name))
            if self.kappa[1]:
                lparam_ref.append("(atom,{})(kappa)".format(self.name))
            if self.modechi == "aniso":
                if self.chi11[1]:
                    lparam_ref.append("(atom,{})(chi11)".format(self.name))
                if self.chi22[1]:
                    lparam_ref.append("(atom,{})(chi22)".format(self.name))
                if self.chi33[1]:
                    lparam_ref.append("(atom,{})(chi33)".format(self.name))
                if self.chi12[1]:
                    lparam_ref.append("(atom,{})(chi12)".format(self.name))
                if self.chi13[1]:
                    lparam_ref.append("(atom,{})(chi13)".format(self.name))
                if self.chi23[1]:
                    lparam_ref.append("(atom,{})(chi23)".format(self.name))
            else:
                if self.chiiso[1]:
                    lparam_ref.append("(atom,{})(chiiso)".format(self.name))
            if self.mag_ia[1]:
                lparam_ref.append("(atom,{})(mag_ia)".format(self.name))
            if self.mag_ib[1]:
                lparam_ref.append("(atom,{})(mag_ib)".format(self.name))
            if self.mag_ic[1]:
                lparam_ref.append("(atom,{})(mag_ic)".format(self.name))
        return lparam_ref

    def take_param_const(self):
        """
        take constrained parameters from model
        """
        lparam_const = []
        if self.coordx[2] != "":
            lparam_const.append("(atom,{})(coordx)".format(self.name))
        if self.coordy[2] != "":
            lparam_const.append("(atom,{})(coordy)".format(self.name))
        if self.coordz[2] != "":
            lparam_const.append("(atom,{})(coordz)".format(self.name))
        if self.bscat[2] != "":
            lparam_const.append("(atom,{})(bscat)".format(self.name))
        if self.occ[2] != "":
            lparam_const.append("(atom,{})(occ)".format(self.name))
        #if self.modedw == "aniso":
        if self.beta11[2] != "":
            lparam_const.append("(atom,{})(beta11)".format(self.name))
        if self.beta22[2] != "":
            lparam_const.append("(atom,{})(beta22)".format(self.name))
        if self.beta33[2] != "":
            lparam_const.append("(atom,{})(beta33)".format(self.name))
        if self.beta12[2] != "":
            lparam_const.append("(atom,{})(beta12)".format(self.name))
        if self.beta13[2] != "":
            lparam_const.append("(atom,{})(beta13)".format(self.name))
        if self.beta23[2] != "":
            lparam_const.append("(atom,{})(beta23)".format(self.name))
        #else:
        if self.biso[2] != "":
            lparam_const.append("(atom,{})(biso)".format(self.name))
        if self.modemagn:
            if self.lfactor[2] != "":
                lparam_const.append("(atom,{})(lfactor)".format(self.name))
            if self.kappa[2] != "":
                lparam_const.append("(atom,{})(kappa)".format(self.name))
            if self.modechi == "aniso":
                if self.chi11[2] != "":
                    lparam_const.append("(atom,{})(chi11)".format(self.name))
                if self.chi22[2] != "":
                    lparam_const.append("(atom,{})(chi22)".format(self.name))
                if self.chi33[2] != "":
                    lparam_const.append("(atom,{})(chi33)".format(self.name))
                if self.chi12[2] != "":
                    lparam_const.append("(atom,{})(chi12)".format(self.name))
                if self.chi13[2] != "":
                    lparam_const.append("(atom,{})(chi13)".format(self.name))
                if self.chi23[2] != "":
                    lparam_const.append("(atom,{})(chi23)".format(self.name))
            else:
                if self.chiiso[2] != "":
                    lparam_const.append("(atom,{})(chiiso)".format(self.name))

            if self.mag_ia[2] != "":
                lparam_const.append("(atom,{})(mag_ia)".format(self.name))
            if self.mag_ib[2] != "":
                lparam_const.append("(atom,{})(mag_ib)".format(self.name))
            if self.mag_ic[2] != "":
                lparam_const.append("(atom,{})(mag_ic)".format(self.name))

        return lparam_const

    def out_to_string(self, ucp):
        """
        output to string
        """
        llsout = []
        llsout.append("")
        llsout.append(70*"*")
        llsout.append("atom {}:".format(self.name))
        llsout.append("type {}".format(self.type))
        llsout.append("\n        x       y       z   occ. mult.  bscat.(cm-12)")
        coord = [self.coordx[0], self.coordy[0], self.coordz[0]]
        llsout.append("  {:.5f} {:.5f} {:.5f}  {:.3f} {:5}  {:.5f}".format(
            coord[0], coord[1], coord[2], self.occ[0], self.mult, self.bscat[0]))
        #if self.modedw == 'aniso':
        llsout.append("\n beta11 beta22 beta33 beta12 beta13 beta23")
        llsout.append(" {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n".format(
                self.beta11[0], self.beta22[0], self.beta33[0], self.beta12[0], self.beta13[0],
                self.beta23[0]))
        #elif self.modedw == 'iso':
        llsout.append("\n    biso = {:6.3f} (Ang**2)".format(self.biso[0]))
        if self.modemagn:
            llsout.append("coeff j0 {0[0]:8.5f} {0[1]:8.5f} {0[2]:8.5f} {0[3]:8.5f} {0[4]:8.5f} {0[5]:8.5f} {0[6]:8.5f}".format(self.coeff0))
            llsout.append("coeff j2 {0[0]:8.5f} {0[1]:8.5f} {0[2]:8.5f} {0[3]:8.5f} {0[4]:8.5f} {0[5]:8.5f} {0[6]:8.5f}".format(self.coeff2))
            llsout.append("lfactor = {:.2f} kappa = {:.2f}".format(self.lfactor[0],self.kappa[0]))
            chian = []
            if self.modechi == 'aniso':
                chian = [self.chi11[0], self.chi22[0], self.chi33[0], self.chi12[0], self.chi13[0],
                         self.chi23[0]]
                llsout.append("\n  chi11  chi22  chi33  chi12  chi13  chi23 (susceptibility, muB/T)")
                llsout.append(" {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f}\n".format(
                    chian[0], chian[1], chian[2], chian[3], chian[4], chian[5]))
            elif self.modechi == 'iso':
                llsout.append("\n  chiiso = {:6.3f} (susceptibility, muB/T)\n".format(self.chiiso[0]))
            condm = self.condm
            if condm != []:
                llsout.append("Symmetry restriction on chi")
                llsout.append("     11  22  33  12  13  23  =  0")
                llsout.append("\n".join(
                    ["    "+"".join(["{:3} ".format(hh2) for hh2 in hh]) for hh in condm]))
            llsout.append("\n magnetic form factor from sthovl (from 0.00 with step 0.01 Ang.):")
            lsthovl = [0.01*hh1 for hh1 in range(100)]
            lmff = [self.calcMFFSpherMod(sthovl, kappa = self.kappa[0]) for sthovl in lsthovl]
            lhelp = [[]for ih1 in range(11)]
            for ih1 in range(10):
                lhelp[ih1+1].append("".join([" {:5.2f}".format(lmff[10*hh1+ih1]) for hh1 in range(10)]))
            for hh1 in lhelp:
                llsout.append("".join(hh1))

        if self.modemagn:
            if chian != []:
                llsout.append("\n equivalent position of atom in unit cell with chi (in muB/T) and     f_mag_0 (in muB), it is not used, just for test")
                llsout.append("       x      y      z    chi11  chi22  chi33  chi12  chi13  chi23    mag_ia mag_ib mag_ic")
                llshelp = [[""], [""], [""]]
                mag = [self.mag_ia[0],self.mag_ib[0],self.mag_ic[0]]
                ihelp = 0
                for els in self.eluniqat:
                    if ihelp >= 3:
                        llshelp.append([" "])
                        llshelp.extend([[""], [""], [""]])
                        ihelp = 0
                    coordeq = ((els[0]+els[1]*coord[0]+ els[2]*coord[1]+ els[3]*coord[2])%1.,
                               (els[4]+els[5]*coord[0]+ els[6]*coord[1]+ els[7]*coord[2])%1.,
                               (els[8]+els[9]*coord[0]+els[10]*coord[1]+els[11]*coord[2])%1.)
                    chir = cfunc.chirot(chian, els)
                    magr=cfunc.magrot(mag,els)
                    llsout.append(" {0[0]:7.3f}{0[1]:7.3f}{0[2]:7.3f}  {1[0][0]:7.3f}{1[1][1]:7.3f}{1[2][2]:7.3f}{1[0][1]:7.3f}{1[0][2]:7.3f}{1[1][2]:7.3f}   {2[0]:7.3f}{2[1]:7.3f}{2[2]:7.3f} ".format(coordeq, chir,magr))
                    chirm = [[0.2695*hh2 for hh2 in hh1] for hh1 in chir]
                    chi_loc = cfunc.chiLOC(chirm, ucp)
                    ihelp += 1
                    for hh, lshelp in zip(chi_loc, llshelp[-3:]):
                        lshelp.append("".join(["{:7.3f}".format(hh2) for hh2 in hh])+"  ")
                llsout.append("\n CHI in orthogonal coordinate system x||a* z||c (in 10**-12cm)")
                for lshelp in llshelp:
                    llsout.append("".join(lshelp))

        else:
            llsout.append("\n equivalent position of atom in unit cell")
            llsout.append("       x      y      z    ")
            for els in self.eluniqat:
                coordeq = ((els[0]+els[1]*coord[0]+ els[2]*coord[1]+ els[3]*coord[2])%1.,
                           (els[4]+els[5]*coord[0]+ els[6]*coord[1]+ els[7]*coord[2])%1.,
                           (els[8]+els[9]*coord[0]+els[10]*coord[1]+els[11]*coord[2])%1.)
                llsout.append(" {0[0]:7.3f}{0[1]:7.3f}{0[2]:7.3f} ".format(coordeq))

        return llsout

    def symmconstr(self, lelsymm):
        """
        symmetry constraint on the thermal vibrations and susceptibility
        """
        mcond, mcond1 = [], []
        m1 = (0, 1, 2, 3, 4, 5)
        m2 = m1
        kl1 = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))
        kl2 = kl1
        for els in lelsymm:
            [b1, r11, r12, r13, b2, r21, r22, r23, b3, r31, r32, r33] = els
            #elemets of rotation matrix is integers: -1, 0, 1 and nothing else.
            mR = [[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]]
            #miR=[[int(round(hh2)) for hh2 in hh1] for hh1 in numpy.linalg.inv(mR)]
            mA = [[0 for hh2 in m2] for hh1 in m1]
            for nm1, nkl1 in zip(m1, kl1):
                i1, j1 = nkl1
                for nm2, nkl2 in zip(m2, kl2):
                    k1, l1 = nkl2
                    mA[nm1][nm2] = mR[i1][k1]*mR[j1][l1]
                    if nm2 > 2:
                        mA[nm1][nm2] += mR[i1][l1]*mR[j1][k1]
                cond1 = sum([abs(hh) if ihh != nm1 else abs(hh-1) for ihh, hh in enumerate(mA[nm1])])
                if cond1 != 0: mcond1.append([hh if ihh != nm1 else hh-1 for ihh, hh in enumerate(mA[nm1])])

        #mCond1 has a lot of element which are repeat each other
        if mcond1 != []:
            for hh in mcond1:
                hh1 = [abs(hh3)  for hh3 in hh if hh3 != 0]
                if hh1 != []:
                    hh2 = [int(hh3/min(hh1)) for hh3 in hh]
                    ihh2 = [ihh3  for ihh3, hh3 in enumerate(hh2) if hh3 != 0]
                    coeff = 1-2*(hh2[ihh2[0]] < 0)
                    hh2 = tuple([coeff*hh3 for hh3 in hh2])
                    if not hh2 in mcond:
                        mcond.append(hh2)
        self.condm = mcond

    def get_j0j2(self,fdirprog):
        fname = "formmag.tab"
        fid = open(os.path.join(fdirprog,fname),"r")
        lcontentH = fid.readlines()
        fid.close()
        lcontent = [cfunc.splitlinewminuses(line) for line in lcontentH if line[0] == "F"]
        self.coeff0 = 7*[0.]
        self.coeff2 = 7*[0.]

        for vals in lcontent:
            if ((vals[1] == self.type.strip())&(vals[2] == 0)):
                self.coeff0 = vals[3:10]
            elif ((vals[1] == self.type.strip())&(vals[2] == 2)):
                self.coeff2 = vals[3:10]


class ccore_ref(object):
    """
    have information of input parameter for refinement
    """
    def __init__(self, cag_ref = None):
        self.output = ""
        self.refin = True
        self.sigmas = True
        self.paramref = None
        self.paramrefval = None
        self.paramconstr = None
        self.paramreferrors = None
        self.hessian = None
        self.ihessian = None
        self.corrM = None
        self.fdirxml = None
        if cag_ref != None:
            self.take_from_agent(cag_ref)

    def set_fdirxml(self, fdirxml):
        self.fdirxml = fdirxml
        
    def take_from_agent(self, cag_ref):
        self.output = cag_ref.output
        self.refin = cag_ref.refin
        self.sigmas = cag_ref.sigmas

    def out_to_string(self):
        """
        prepare information for output file
        """
        llsout=[]
        if (self.paramreferrors!=None):
            llsout.append("\n")
            llsout.append(70*"*")
            llsout.append("Refined parameters /values and errorbars")
            llsout.extend(["{:3}. {:40}{:12.7f}{:12.7f}".format(inum+1,name,val,err) for inum,name,val,err in zip(range(len(self.paramref)),self.paramref,self.paramrefval,self.paramreferrors)])
        if (self.hessian!=None):
            llsout.append("\nHessian of refined parameters")
            llsout.append(4*" "+"".join(["{:9}.".format(inum+1) for inum in range(len(self.hessian))]))
            for inum,hh in enumerate(self.hessian):
                llsout.append("{:3}.".format(inum+1)+"".join(["{:10.0f}".format(hh2) for hh2 in hh]))
        if (self.ihessian!=None):
            llsout.append("\nInversed hessian of refined parameters")
            llsout.append(4*" "+"".join(["{:9}.".format(inum+1) for inum in range(len(self.ihessian))]))
            for inum,hh in enumerate(self.ihessian):
                llsout.append("{:3}.".format(inum+1)+"".join(["{:10.7f}".format(hh2) for hh2 in hh]))
        if (self.corrM!=None):
            llsout.append("\nCorrelation matrix of refined parameters")
            llsout.append(4*" "+"".join(["{:6}.".format(inum+1) for inum in range(len(self.corrM))]))
            for inum,hh in enumerate(self.corrM):
                llsout.append("{:3}.".format(inum+1)+"".join(["{:7.3f}".format(hh2) for hh2 in hh]))
        llsout.append(70*"*")
        return llsout



