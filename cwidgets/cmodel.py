"""
define main class cagent

subclasses:
    ag_ph, ag_epx, ag_ph_at, ag_exp_ph, ag_ref

it is contained the information which are presented in the parameter xml file.
It is a bridge between GUI and core parameters

"""

#import sys
#from PyQt5 import QtWidgets
#from PyQt5 import QtGui
#from PyQt5 import QtCore

import os
import xml
import xml.etree.ElementTree

import widg_min
import read_cif
import load_bscat_from_tab
import xml2var
import ccore

class cmodel(widg_min.cmodel_min):
    """
    minimal model for cwidget
    """
    def __init__(self):
        super(cmodel, self).__init__()
        self.ph = []
        self.exp = []
        self.ref =  [cmodel_ref()]
        self.fdirprog = None
        self.fdirxml = None
        self.fxml = None
        self.valchi2 = 0.
        self.n_total = 1.
    def add_ph(self):
        ph = cmodel_ph()
        ph.fdirprog = self.fdirprog
        if (not ph in self.ph):
            self.ph.append(ph)
            self.mod_changed()
    def del_ph(self, ph):
        self.ph.remove(ph)
        self.mod_changed()
    def add_exp(self):
        exp = cmodel_exp()
        if (not exp in self.exp):
            self.exp.append(exp)
            self.set_fdirprog(self.fdirprog)
            self.set_fdirxml(self.fdirxml)
            self.set_fxml(self.fxml)
            self.mod_changed()
    def del_exp(self, exp):
        self.ph.remove(exp)
        self.mod_changed()
    def del_obs(self, obs):
        """
        to add to definition of observer
        """
        super(cmodel, self).del_obs(obs)
        for ph in self.ph:
            ph.del_obs(obs)
        for exp in self.exp:
            exp.del_obs(obs)
        for ref in self.ref:
            ref.del_obs(obs)
    def del_builder(self, build):
        super(cmodel, self).del_builder(build)
        for ph in self.ph:
            ph.del_builder(build)
        for exp in self.exp:
            exp.del_builder(build)
        for ref in self.ref:
            ref.del_builder(build)


    def take_from_core(self, ccore):
        lexp = []
        for hh in ccore.exp:
            exp = cmodel_exp()
            exp.take_from_core(hh)
            lexp.append(exp)
        self.exp = lexp

        lph = []
        for hh in ccore.ph:
            ph = cmodel_ph()
            ph.take_from_core(hh)
            lph.append(ph)
        self.ph = lph

        lref = []
        for hh in ccore.ref:
            ref = cmodel_ref()
            ref.take_from_core(hh)
            lref.append(ref)
        self.ref = lref

        self.set_fdirxml(ccore.fdirxml)
        self.set_fdirprog(ccore.fdirprog)
        self.set_fxml(ccore.fxml)
        self.valchi2 = ccore.valchi2
        self.n_total = ccore.n_total

    def load_xml(self,fdirxml,fxml):
        """
        load the data from xml to file.
        Should be deleted
        """
        self.fdirxml = fdirxml
        self.fxml = fxml

        #structure of xml file
        sparameters = [['exp', 'list'], ['ph', 'list'], ['ref', 'list']]
        sexp = [['name', 'atr'], ['mode', 'atr'], ['output', 'atr'],
                ['input', 'atr'], ['bkgr', 'atr'], ['chi2', 'atr'], ['excl2theta', 'atr'],
                ['zeroshift', 'val'],['zshift_a', 'val'],
                ['p_asym_1', 'val'], ['p_asym_2', 'val'],['p_asym_3', 'val'], ['p_asym_4', 'val'],
                ['pup', 'val'], ['pdown', 'val'], ['x', 'val'], ['y', 'val'],
                ['U', 'val'], ['V', 'val'], ['W', 'val'], ['exp_ph', 'list']]
        sphaseexp = [['name', 'atr'], ['scale', 'val'],['igsize', 'val'], ['ext11', 'val'], ['ext22', 'val'],
                      ['ext33', 'val'], ['ext12', 'val'], ['ext13', 'val'], ['ext23', 'val']]
        sref = [['refin', 'atr'], ['output', 'atr'],['sigmas', 'atr']]
        sph = [['name', 'atr'], ['spgr', 'atr'], ['cella', 'val'], ['cellb', 'val'],
                  ['cellc', 'val'], ['cellalpha', 'val'], ['cellbeta', 'val'], ['cellgamma', 'val'],
                  ['atom', 'list']]
        satom = [['modedw', 'atr'], ['modechi', 'atr'], ['name', 'atr'], ['type', 'atr'], ['coeff0', 'atr'],
                 ['coeff2', 'atr'], ['magnet', 'atr'], ['bscat', 'val'], ['occ', 'val'],
                 ['coordx', 'val'], ['coordy', 'val'], ['coordz', 'val'], ['lfactor', 'val'],
                 ['kappa', 'val'],
                 ['chi11', 'val'], ['chi22', 'val'], ['chi33', 'val'], ['chi12', 'val'],
                 ['chi13', 'val'], ['chi23', 'val'], ['chiiso', 'val'], ['biso', 'val'],
                 ['mag_ia', 'val'], ['mag_ib', 'val'], ['mag_ic', 'val'], ['al_field', 'val'],
                 ['tol_ising', 'val'],
                 ['beta11', 'val'], ['beta22', 'val'], ['beta33', 'val'], ['beta12', 'val'],
                 ['beta13', 'val'], ['beta23', 'val']]
        dstruct = {'parameters':sparameters, 'exp':sexp, 'ref':sref,
                   'ph':sph, 'exp_ph':sphaseexp, 'atom':satom}

        self.dstruct = dstruct

        tree = xml.etree.ElementTree.parse(os.path.join(fdirxml, fxml))
        root = tree.getroot()
        dout = xml2var.readxml(root, self.dstruct, root)


        ag_ref = cmodel_ref()
        if "ref" in dout.keys():
            drefinement = dout['ref'][0]
            ag_ref.load_dict(drefinement)
        #cref.fdirxml = fdirxml

        lag_exp, lag_ph = [], []
        if "ph" in dout.keys():
            for dphase in dout['ph']:
                ag_ph = cmodel_ph()
                ag_ph.load_dict(dphase)
                ag_ph.fdirprog = self.fdirprog
                lag_ph.append(ag_ph)

        if "exp" in dout.keys():
            for dexperiment in dout['exp']:
                ag_exp = cmodel_exp()
                ag_exp.load_dict(dexperiment)
                ag_exp.fdirxml = fdirxml
                ag_exp.fdirprog = self.fdirprog
                ag_exp.fxml = fxml

                lag_exp.append(ag_exp)

        self.ph = lag_ph
        self.exp = lag_exp
        self.ref = [ag_ref]


    def save_to_xml(self, fxml):
        """
        save files to xml file
        """
        self.fxml = fxml
        fpar = os.path.join(self.fdirxml, fxml)
        stype = 'parameters'
        xml_elem = xml.etree.ElementTree.Element(stype)
        xml_refin = self.ref[0].to_xml_elem()
        xml_elem.insert(0, xml_refin)
        for ag_ph in reversed(self.ph):
            xml_phase = ag_ph.to_xml_elem()
            xml_elem.insert(0, xml_phase)
        for ag_exp in reversed(self.exp):
            xml_exp = ag_exp.to_xml_elem()
            xml_elem.insert(0, xml_exp)
        xml_tree = xml.etree.ElementTree.ElementTree(xml_elem)
        xml_tree.write(fpar)

    def set_fxml(self, fxml):
        self.fxml = fxml
        for exp in self.exp:
            exp.fxml = fxml

    def set_fdirxml(self, fdirxml):
        self.fdirxml = fdirxml
        for exp in self.exp:
            exp.fdirxml = fdirxml
        for ref in self.ref:
            ref.fdirxml = fdirxml
            

    def set_fdirprog(self, fdirprog):
        self.fdirprog = fdirprog
        for exp in self.exp:
            exp.fdirprog = fdirprog
        for ph in self.ph:
            ph.fdirprog = fdirprog

    def check_phaseexp(self):
        """
        If experiment does not content information about phases but there is phase, then the first one is added to phaseexp
        """
        if self.ph == []:
            return
        lname = [hh.name for hh in self.ph]
        for ag_exp in self.exp:
            if ag_exp.exp_ph == []:
                ag_phaseexp = cmodel_exp_ph()
                ag_phaseexp.name = lname[0]
                ag_exp.exp_ph.append(ag_phaseexp)
        for ag_exp in self.exp:
            for iphexp, ag_phaseexp  in enumerate(ag_exp.exp_ph):
                if (not (ag_phaseexp.name in lname)):
                    iph = min([iphexp, len(lname)-1])
                    ag_phaseexp.name = lname[iph]

    def get_constr_ucp(self):
        lconstr = []
        for ph in self.ph:
            lconstr_ph = ph.get_constr_ucp()
            lconstr.extend(lconstr_ph)
        return lconstr
    
    def get_constr_adp(self):
        lconstr = []
        for ph in self.ph:
            lconstr_ph = ph.get_constr_adp()
            lconstr.extend(lconstr_ph)
        return lconstr
    
    def get_link_on_constr(self):
        l_link = []
        core = ccore.ccore()
        core.take_from_agent(self)
        l_link = core.take_param_const()
        return l_link

class cmodel_ph(widg_min.cmodel_min):
    """
    minimal model for cwidget phase
    """
    def __init__(self):
        super(cmodel_ph, self).__init__()
        self.name = "phase"
        self.spgr = "P1"
        self.spgr_code = "1"
        self.cella = [1., False, ""]
        self.cellb = [1., False, ""]
        self.cellc = [1., False, ""]
        self.cellalpha = [90.,False, ""]
        self.cellbeta = [90.,False, ""]
        self.cellgamma = [90.,False, ""]
        self.atom = []
        self.fdirprog = None
    def del_obs(self, obs):
        super(cmodel_ph, self).del_obs(obs)
        for at in self.atom:
            at.del_obs(obs)
    def del_builder(self, build):
        super(cmodel_ph, self).del_builder(build)
        for at in self.atom:
            at.del_builder(build)
    def add_at(self):
        at = cmodel_at()
        if (not at in self.atom):
            self.atom.append(at)
            self.mod_changed()
    def del_at(self, lat):
        for at in lat:
            self.atom.remove(at)
        self.mod_changed()


    def take_from_core(self,ccore_ph):
        self.name = ccore_ph.name
        self.spgr = ccore_ph.spgr
        self.spgr_code = ccore_ph.spgr_code
        self.cella = [hh for hh in ccore_ph.cella]
        self.cellb = [hh for hh in ccore_ph.cellb]
        self.cellc = [hh for hh in ccore_ph.cellc]
        self.cellalpha = [hh for hh in ccore_ph.cellalpha]
        self.cellbeta = [hh for hh in ccore_ph.cellbeta]
        self.cellgamma = [hh for hh in ccore_ph.cellgamma]
        latom = []
        for hh in ccore_ph.atom:
            atom = cmodel_at()
            atom.take_from_core(hh)
            latom.append(atom)
        self.atom = latom
        self.fdirprog = ccore_ph.fdirprog


    def load_dict(self,dphase):
        """
        load info about phases from dictionary (take from xml file) to variable cphase
        """
        if 'name' in dphase.keys(): self.name = dphase["name"]
        if 'spgr' in dphase.keys(): self.spgr = dphase["spgr"]
        if 'cella' in dphase.keys(): self.cella = dphase['cella']
        if 'cellb' in dphase.keys(): self.cellb = dphase['cellb']
        if 'cellc' in dphase.keys(): self.cellc = dphase['cellc']
        if 'cellalpha' in dphase.keys(): self.cellalpha = dphase['cellalpha']
        if 'cellbeta' in dphase.keys(): self.cellbeta = dphase['cellbeta']
        if 'cellgamma' in dphase.keys(): self.cellgamma = dphase['cellgamma']

        lag_ph_at = []
        if 'atom' in dphase.keys():
            for datom in dphase['atom']:
                ag_ph_at = cmodel_at()
                ag_ph_at.load_dict(datom)
                lag_ph_at.append(ag_ph_at)
        self.atom = lag_ph_at

    def to_xml_elem(self):
        xml_elem = xml.etree.ElementTree.Element('ph')
        xml_elem.attrib["name"]=self.name
        xml_elem.attrib["spgr"]=self.spgr

        def tempfunc4val(xmlelem,valname,val):
             xmlSubElem=xml.etree.ElementTree.SubElement(xmlelem,valname)
             xmlSubElem.attrib["refined"]="false" if ((not(val[1]))|(val[2]!="")) else "true"
             xmlSubElem.attrib["constr"]="true" if (val[2]!="") else "false"
             xmlSubElem.text=val[2] if (val[2]!="") else "{}".format(val[0])
             return

        tempfunc4val(xml_elem,'cella',self.cella)
        tempfunc4val(xml_elem,'cellb',self.cellb)
        tempfunc4val(xml_elem,'cellc',self.cellc)
        tempfunc4val(xml_elem,'cellalpha',self.cellalpha)
        tempfunc4val(xml_elem,'cellbeta', self.cellbeta)
        tempfunc4val(xml_elem,'cellgamma',self.cellgamma)
        for ag_ph_at in self.atom:
             xml_elem_at = ag_ph_at.to_xml_elem()
             xml_elem.insert(-1, xml_elem_at)
        return xml_elem

    def load_from_cif(self, fcif):
        """
        call internal program to load info about crystall structure from file .cif
        """
        lphase_cif = read_cif.read_cif(fcif)
        if len(lphase_cif) == 0:
            print "phase is not found"
        phase_cif = lphase_cif[0]
        if phase_cif.name == None:
            phase_cif.name = os.path.basename(fcif)
        self.name = phase_cif.name
        self.spgr = "".join(phase_cif.spgr.split(" "))
        self.cella = [phase_cif.cella, False, ""]
        self.cellb = [phase_cif.cellb, False, ""]
        self.cellc = [phase_cif.cellc, False, ""]
        self.cellalpha = [phase_cif.cellalpha, False, ""]
        self.cellbeta = [phase_cif.cellbeta, False, ""]
        self.cellgamma = [phase_cif.cellgamma, False, ""]
        latom = []
        for hh in phase_cif.atom:
            atom = cmodel_at()
            atom.load_from_atom_cif(hh)
            latom.append(atom)
        self.atom = latom
        try: 
            self.load_bscat()
        except:
            pass
        self.mod_changed()

    def load_bscat(self):
        """
        read bscatterings from file
        """
        lname = [hh.name for hh in self.atom]
        lbscat = load_bscat_from_tab.read_bscat_from_tab(lname)
        for hh, bscat in zip(self.atom, lbscat):
            hh.bscat = [bscat, False, ""]
        self.mod_changed()

    def get_constr_ucp(self):
        core_ph = ccore.ccore_ph()
        core_ph.take_from_agent(self)
        lconstr = core_ph.get_constr_ucp()
        return lconstr
    
    def get_constr_adp(self):
        core_ph = ccore.ccore_ph()
        core_ph.take_from_agent(self)
        lconstr = core_ph.get_constr_adp()
        return lconstr


class cmodel_exp(widg_min.cmodel_min):
    """
    minimal model for cwidget experiment
    """
    def __init__(self):
        super(cmodel_exp, self).__init__()
        self.name = "exp"
        self.input = "exp.dat"
        self.output = "exp.out"
        self.bkgr = "exp.bkgr"
        self.powder = True
        self.pup = [1., False, ""]
        self.pdown = [1., False, ""]
        self.U = [0., False, ""]
        self.V = [0., False, ""]
        self.W = [0.1, False, ""]
        self.x = [0., False, ""]
        self.y = [0., False, ""]
        self.zeroshift = [0.,False,""]
        self.zshift_a  = [0.,False,""]
        self.p_asym_1  = [0.,False,""]
        self.p_asym_2  = [0.,False,""]
        self.p_asym_3  = [0.,False,""]
        self.p_asym_4  = [0.,False,""]
        self.excl2theta = "(0.0,1.0)"
        self.exp_ph = []
        self.modechi2_up = True
        self.modechi2_diff = False
        self.modechi2_down = True
        self.valchi2 = 0.
        self.valchi2up = 0.
        self.valchi2down = 0.
        self.valchi2diff = 0
        self.n_total = 0
        self.fdirprog = None
        self.fdirxml = None
        self.fxml = None
        self.tth_min = "0.1"
        self.tth_max = "179.9"
        self.phi_min = "-89.9"
        self.phi_max = "89.9"
        self.mode_2dpd = False         
        
        
    def del_obs(self, obs):
        super(cmodel_exp, self).del_obs(obs)
        for exp_ph in self.exp_ph:
            exp_ph.del_obs(obs)
    def del_builder(self, build):
        super(cmodel_exp, self).del_builder(build)
        for exp_ph in self.exp_ph:
            exp_ph.del_builder(build)
    def add_exp_ph(self):
        exp_ph = cmodel_exp_ph()
        if (not exp_ph in self.exp_ph):
            self.exp_ph.append(exp_ph)
            self.mod_changed()
    def del_exp_ph(self, lexp_ph):
        for exp_ph in lexp_ph:
            self.exp_ph.remove(exp_ph)
        self.mod_changed()

    def take_from_core(self, ccore_exp):
        self.name = ccore_exp.name
        self.input = ccore_exp.input
        self.output = ccore_exp.output
        self.bkgr = ccore_exp.bkgr
        self.powder = ccore_exp.powder
        self.pup = [hh for hh in ccore_exp.pup]
        self.pdown = [hh for hh in ccore_exp.pdown]
        self.U = [hh for hh in ccore_exp.U]
        self.V = [hh for hh in ccore_exp.V]
        self.W = [hh for hh in ccore_exp.W]
        self.x = [hh for hh in ccore_exp.x]
        self.y = [hh for hh in ccore_exp.y]
        self.zeroshift = [hh for hh in ccore_exp.zeroshift]
        self.zshift_a  = [hh for hh in ccore_exp.zshift_a]
        self.p_asym_1  = [hh for hh in ccore_exp.p_asym_1]
        self.p_asym_2  = [hh for hh in ccore_exp.p_asym_2]
        self.p_asym_3  = [hh for hh in ccore_exp.p_asym_3]
        self.p_asym_4  = [hh for hh in ccore_exp.p_asym_4]
        self.excl2theta = ccore_exp.excl2theta
        self.tth_min = ccore_exp.tth_min
        self.tth_max = ccore_exp.tth_max
        self.phi_min = ccore_exp.phi_min
        self.phi_max = ccore_exp.phi_max
        self.mode_2dpd = ccore_exp.mode_2dpd         
        l_exp_ph = []
        for hh in ccore_exp.exp_ph:
            exp_ph = cmodel_exp_ph()
            exp_ph.take_from_core(hh)
            l_exp_ph.append(exp_ph)
        self.exp_ph = l_exp_ph
        self.modechi2_up = ccore_exp.modechi2_up
        self.modechi2_down = ccore_exp.modechi2_down
        self.modechi2_diff = ccore_exp.modechi2_diff
        self.valchi2 = ccore_exp.valchi2
        self.valchi2up = ccore_exp.valchi2up
        self.valchi2down = ccore_exp.valchi2down
        self.valchi2diff = ccore_exp.valchi2diff
        self.n_total = ccore_exp.n_total
        self.fdirprog = ccore_exp.fdirprog
        self.fdirxml = ccore_exp.fdirxml
        self.fxml = ccore_exp.fxml

    def load_dict(self,dexperiment):
        if 'name' in dexperiment.keys(): self.name = dexperiment["name"]
        if 'input' in dexperiment.keys(): self.input = dexperiment["input"]
        if 'output' in dexperiment.keys(): self.output = dexperiment["output"]
        if 'mode' in dexperiment.keys(): self.powder = (dexperiment["mode"] == "powder")

        if 'pup' in dexperiment.keys(): self.pup = dexperiment["pup"]
        if 'pdown' in dexperiment.keys(): self.pdown = dexperiment["pdown"]

        if not self.powder:
            pass
        else:
            if 'bkgr' in dexperiment.keys(): self.bkgr = dexperiment["bkgr"]
            if 'U' in dexperiment.keys(): self.U = dexperiment["U"]
            if 'V' in dexperiment.keys(): self.V = dexperiment["V"]
            if 'W' in dexperiment.keys(): self.W = dexperiment["W"]
            if 'x' in dexperiment.keys(): self.x = dexperiment["x"]
            if 'y' in dexperiment.keys(): self.y = dexperiment["y"]
            if 'zeroshift' in dexperiment.keys(): self.zeroshift = dexperiment["zeroshift"]
            if 'zshift_a' in dexperiment.keys(): self.zshift_a = dexperiment["zshift_a"]
            if 'p_asym_1' in dexperiment.keys(): self.p_asym_1 = dexperiment["p_asym_1"]
            if 'p_asym_2' in dexperiment.keys(): self.p_asym_2 = dexperiment["p_asym_2"]
            if 'p_asym_3' in dexperiment.keys(): self.p_asym_3 = dexperiment["p_asym_3"]
            if 'p_asym_4' in dexperiment.keys(): self.p_asym_4 = dexperiment["p_asym_4"]
            if 'excl2theta' in dexperiment.keys(): self.excl2theta = dexperiment["excl2theta"]
            lhelp = []
            if 'chi2' in dexperiment.keys(): lhelp = dexperiment["chi2"].strip().split(",")
            lcalcchi2 = [("up" in lhelp), ("down" in lhelp), ("diff" in lhelp)]
            self.modechi2_up = lcalcchi2[0]
            self.modechi2_down = lcalcchi2[1]
            self.modechi2_diff = lcalcchi2[2]
        lphexp = []
        if 'exp_ph' in dexperiment.keys():
            for dphaseexp in dexperiment["exp_ph"]:
                ag_exp_ph = cmodel_exp_ph()
                ag_exp_ph.load_dict(dphaseexp)
                lphexp.append(ag_exp_ph)
        self.exp_ph = lphexp

    def to_xml_elem(self):
        xml_elem = xml.etree.ElementTree.Element('exp')
        xml_elem.attrib["name"]=self.name
        xml_elem.attrib["input"]=self.input
        xml_elem.attrib["output"]=self.output
        xml_elem.attrib["mode"]="powder" if (self.powder) else "mono"

        def tempfunc4val(xml_elem,valname,val):
             xmlSubElem=xml.etree.ElementTree.SubElement(xml_elem,valname)
             xmlSubElem.attrib["refined"]="false" if ((not(val[1]))|(val[2]!="")) else "true"
             xmlSubElem.attrib["constr"]="true" if (val[2]!="") else "false"
             xmlSubElem.text=val[2] if (val[2]!="") else "{}".format(val[0])
             return

        tempfunc4val(xml_elem, 'pup', self.pup)
        tempfunc4val(xml_elem, 'pdown', self.pdown)

        if (not(self.powder)):
            for ag_ph_exp in self.exp_ph:
                xml_ph_exp = ag_ph_exp.to_xml_elem()
                xml_elem.insert(-1, xml_ph_exp)
        else:
            xml_elem.attrib["bkgr"]=self.bkgr
            xml_elem.attrib["excl2theta"]=self.excl2theta
            shelp=""
            for hh1,hh2 in zip([self.modechi2_up, self.modechi2_down, self.modechi2_diff],['up,','down,','diff,']):
                if (hh1): shelp+=hh2
            xml_elem.attrib["chi2"]=shelp[:-1]
            tempfunc4val(xml_elem,'U',self.U)
            tempfunc4val(xml_elem,'V',self.V)
            tempfunc4val(xml_elem,'W',self.W)
            tempfunc4val(xml_elem,'x',self.x)
            tempfunc4val(xml_elem,'y',self.y)
            tempfunc4val(xml_elem,'zeroshift',self.zeroshift)
            tempfunc4val(xml_elem,'zshift_a',self.zshift_a)
            tempfunc4val(xml_elem,'p_asym_1',self.p_asym_1)
            tempfunc4val(xml_elem,'p_asym_2',self.p_asym_2)
            tempfunc4val(xml_elem,'p_asym_3',self.p_asym_3)
            tempfunc4val(xml_elem,'p_asym_4',self.p_asym_4)

            for ag_ph_exp in self.exp_ph:
                xml_ph_exp = ag_ph_exp.to_xml_elem()
                xml_elem.insert(-1, xml_ph_exp)
        return xml_elem







class cmodel_at(widg_min.cmodel_min):
    """
    minimal model for cwidget phase
    """
    def __init__(self):
        super(cmodel_at, self).__init__()
        self.name = "O"
        self.type = 'O2'
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
        self.biso = [0.0000, False, ""]
        self.beta11, self.beta22, self.beta33 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.beta12, self.beta13, self.beta23 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.kappa = [1.0000, False, ""]
        self.chiiso = [0.0000, False, ""]
        self.chi11, self.chi22, self.chi33 = [0., False, ""], [0., False, ""], [0., False, ""]
        self.chi12, self.chi13, self.chi23 = [0., False, ""], [0., False, ""], [0., False, ""]

        self.mag_ia, self.mag_ib, self.mag_ic= [0., False, ""], [0., False, ""], [0., False, ""]
        self.al_field, self.tol_ising = [1., False, ""], [0.2, False, ""]

    def take_from_core(self,ccore_ph_at):
        self.name = ccore_ph_at.name
        self.type = ccore_ph_at.type
        self.type_n = ccore_ph_at.type_n
        self.modedw = ccore_ph_at.modedw
        self.modechi = ccore_ph_at.modechi
        self.modemagn = ccore_ph_at.modemagn
        self.coordx = [hh for hh in ccore_ph_at.coordx]
        self.coordy = [hh for hh in ccore_ph_at.coordy]
        self.coordz = [hh for hh in ccore_ph_at.coordz]
        self.bscat = [hh for hh in ccore_ph_at.bscat]
        self.occ = [hh for hh in ccore_ph_at.occ]
        self.lfactor = [hh for hh in ccore_ph_at.lfactor]
        self.biso = [hh for hh in ccore_ph_at.biso]
        self.beta11, self.beta22, self.beta33 = [hh for hh in ccore_ph_at.beta11], [hh for hh in ccore_ph_at.beta22], [hh for hh in ccore_ph_at.beta33]
        self.beta12, self.beta13, self.beta23 = [hh for hh in ccore_ph_at.beta12], [hh for hh in ccore_ph_at.beta13], [hh for hh in ccore_ph_at.beta23]
        self.kappa = [hh for hh in ccore_ph_at.kappa]
        self.chiiso = [hh for hh in ccore_ph_at.chiiso]
        self.chi11, self.chi22, self.chi33 = [hh for hh in ccore_ph_at.chi11], [hh for hh in ccore_ph_at.chi22], [hh for hh in ccore_ph_at.chi33]
        self.chi12, self.chi13, self.chi23 = [hh for hh in ccore_ph_at.chi12], [hh for hh in ccore_ph_at.chi13], [hh for hh in ccore_ph_at.chi23]
        self.mag_ia, self.mag_ib, self.mag_ic= [hh for hh in ccore_ph_at.mag_ia], [hh for hh in ccore_ph_at.mag_ib], [hh for hh in ccore_ph_at.mag_ic]
        self.al_field, self.tol_ising = [hh for hh in ccore_ph_at.al_field], [hh for hh in ccore_ph_at.tol_ising]

    def load_dict(self, datom):
        """
        load info about atoms from dictionary (take from xml file) to variable CAtom
        """
        if "name" in datom.keys(): self.name = datom['name']
        if "type" in datom.keys(): self.type = datom['type']
        if "type_n" in datom.keys(): self.type_n = datom['type_n']
        if 'bscat' in datom.keys(): self.bscat = datom['bscat']

        if 'occ' in datom.keys(): self.occ = datom['occ']
        if 'coordx' in datom.keys(): self.coordx = datom['coordx']
        if 'coordy' in datom.keys(): self.coordy = datom['coordy']
        if 'coordz' in datom.keys(): self.coordz = datom['coordz']

        if 'modedw' in datom.keys():
            self.modedw = datom['modedw']
            #if datom['modedw'] == 'aniso':
            if 'beta11' in datom.keys(): self.beta11 = datom['beta11']
            if 'beta22' in datom.keys(): self.beta22 = datom['beta22']
            if 'beta33' in datom.keys(): self.beta33 = datom['beta33']
            if 'beta12' in datom.keys(): self.beta12 = datom['beta12']
            if 'beta13' in datom.keys(): self.beta13 = datom['beta13']
            if 'beta23' in datom.keys(): self.beta23 = datom['beta23']
            #elif datom['modedw'] == 'iso':
            if 'biso' in datom.keys(): self.biso = datom['biso']
            #else:
            #    print "Mistake in subroutine loadparam"
        if 'magnet' in datom.keys():
            self.modemagn = (datom['magnet'] == "true")
            if datom['magnet'] == "true":
                if 'lfactor' in datom.keys(): self.lfactor = datom['lfactor']
                if 'kappa' in datom.keys(): self.kappa = datom['kappa']
                if 'modechi' in datom.keys():
                    self.modechi = datom['modechi']
                    if datom['modechi'] == 'aniso':
                        if 'chi11' in datom.keys(): self.chi11 = datom['chi11']
                        if 'chi22' in datom.keys(): self.chi22 = datom['chi22']
                        if 'chi33' in datom.keys(): self.chi33 = datom['chi33']
                        if 'chi12' in datom.keys(): self.chi12 = datom['chi12']
                        if 'chi13' in datom.keys(): self.chi13 = datom['chi13']
                        if 'chi23' in datom.keys(): self.chi23 = datom['chi23']
                    elif datom['modechi'] == 'iso':
                        if 'chiiso' in datom.keys():  self.chiiso = datom['chiiso']
                    else:
                        print "Mistake in load_dict subroutine of CAtom"
                if 'mag_ia' in datom.keys(): self.mag_ia = datom['mag_ia']
                if 'mag_ib' in datom.keys(): self.mag_ib = datom['mag_ib']
                if 'mag_ic' in datom.keys(): self.mag_ic = datom['mag_ic']
                if 'al_field' in datom.keys(): self.al_field = datom['al_field']
                if 'tol_ising' in datom.keys(): self.tol_ising = datom['tol_ising']

    def to_xml_elem(self):
        """
        to construct xml object for atom
        """
        xml_elem = xml.etree.ElementTree.Element('atom')
        xml_elem.attrib["name"] = self.name
        xml_elem.attrib["type"] = self.type

        def tempfunc4val(xml_elem, val_name, val):
            """
            temp func
            """
            xml_sub_elem = xml.etree.ElementTree.SubElement(xml_elem, val_name)
            xml_sub_elem.attrib["refined"] = "false" if ((not val[1])|(val[2] != "")) else "true"
            xml_sub_elem.attrib["constr"] = "true" if (val[2] != "") else "false"
            xml_sub_elem.text = val[2] if (val[2] != "") else "{}".format(val[0])
            return

        tempfunc4val(xml_elem, 'coordx', self.coordx)
        tempfunc4val(xml_elem, 'coordy', self.coordy)
        tempfunc4val(xml_elem, 'coordz', self.coordz)
        tempfunc4val(xml_elem, 'bscat', self.bscat)
        tempfunc4val(xml_elem, 'occ', self.occ)

        xml_elem.attrib["modedw"] = self.modedw
        #if self.modedw == "aniso":
        tempfunc4val(xml_elem, 'beta11', self.beta11)
        tempfunc4val(xml_elem, 'beta22', self.beta22)
        tempfunc4val(xml_elem, 'beta33', self.beta33)
        tempfunc4val(xml_elem, 'beta12', self.beta12)
        tempfunc4val(xml_elem, 'beta13', self.beta13)
        tempfunc4val(xml_elem, 'beta23', self.beta23)
        #elif self.modedw == "iso":
        tempfunc4val(xml_elem, 'biso', self.biso)
        #else:
        #    print "Mistake in CAtom, method to_xml_elem."
        xml_elem.attrib["magnet"] = "true" if (self.modemagn) else "false"
        xml_elem.attrib["modechi"] = self.modechi
        if self.modemagn:
            tempfunc4val(xml_elem, 'lfactor', self.lfactor)
            tempfunc4val(xml_elem, 'kappa', self.kappa)
            if self.modechi == "aniso":
                tempfunc4val(xml_elem, 'chi11', self.chi11)
                tempfunc4val(xml_elem, 'chi22', self.chi22)
                tempfunc4val(xml_elem, 'chi33', self.chi33)
                tempfunc4val(xml_elem, 'chi12', self.chi12)
                tempfunc4val(xml_elem, 'chi13', self.chi13)
                tempfunc4val(xml_elem, 'chi23', self.chi23)
            elif self.modechi == "iso":
                tempfunc4val(xml_elem, 'chiiso', self.chiiso)
            else:
                print "Mistake in cmodel_at, method to_xml_elem."
            tempfunc4val(xml_elem, 'mag_ia', self.mag_ia)
            tempfunc4val(xml_elem, 'mag_ib', self.mag_ib)
            tempfunc4val(xml_elem, 'mag_ic', self.mag_ic)
            tempfunc4val(xml_elem, 'al_field', self.al_field)
            tempfunc4val(xml_elem, 'tol_ising', self.tol_ising)
            
        return xml_elem

    def load_from_atom_cif(self, atom_cif):
        """
        load information from atom_cif to cag_ph_at
        """
        self.name = atom_cif.name
        self.type = atom_cif.type
        self.type_n = atom_cif.type
        self.coordx = [atom_cif.coordx, False, ""]
        self.coordy = [atom_cif.coordy, False, ""]
        self.coordz = [atom_cif.coordz, False, ""]
        self.occ = [atom_cif.occ, False, ""]
        self.biso = [atom_cif.biso, False, ""]
        self.beta11 = [atom_cif.beta11, False, ""]
        self.beta22 = [atom_cif.beta22, False, ""]
        self.beta33 = [atom_cif.beta33, False, ""]
        self.beta12 = [atom_cif.beta12, False, ""]
        self.beta13 = [atom_cif.beta13, False, ""]
        self.beta23 = [atom_cif.beta23, False, ""]
        



class cmodel_exp_ph(widg_min.cmodel_min):
    """
    minimal model for cwidget experiment
    """
    def __init__(self):
        super(cmodel_exp_ph, self).__init__()
        self.name = "phase"
        self.scale = [1.00000, False, ""]
        self.igsize = [0.00000, False, ""]
        self.ext11 = [0.00000, False, ""]
        self.ext22 = [0.00000, False, ""]
        self.ext33 = [0.00000, False, ""]
        self.ext12 = [0.00000, False, ""]
        self.ext13 = [0.00000, False, ""]
        self.ext23 = [0.00000, False, ""]

    def take_from_core(self,ccore_exp_ph):
        self.name = ccore_exp_ph.name
        self.scale = [hh for hh in ccore_exp_ph.scale]
        self.igsize = [hh for hh in ccore_exp_ph.igsize]
        self.ext11 = [hh for hh in ccore_exp_ph.ext11]
        self.ext22 = [hh for hh in ccore_exp_ph.ext22]
        self.ext33 = [hh for hh in ccore_exp_ph.ext33]
        self.ext12 = [hh for hh in ccore_exp_ph.ext12]
        self.ext13 = [hh for hh in ccore_exp_ph.ext13]
        self.ext23 = [hh for hh in ccore_exp_ph.ext23]

    def load_dict(self,dphexp):
        if 'name' in dphexp.keys(): self.name = dphexp["name"]
        if 'scale' in dphexp.keys(): self.scale = dphexp["scale"]
        if 'igsize' in dphexp.keys(): self.igsize = dphexp["igsize"]
        if "ext11" in dphexp.keys():
            self.ext11 = dphexp["ext11"]
            self.ext22 = dphexp["ext22"]
            self.ext33 = dphexp["ext33"]
            self.ext12 = dphexp["ext12"]
            self.ext13 = dphexp["ext13"]
            self.ext23 = dphexp["ext23"]

    def to_xml_elem(self):
        xml_elem=xml.etree.ElementTree.Element('exp_ph')
        xml_elem.attrib["name"]=self.name

        def tempfunc4val(xml_elem,valname,val):
             xml_sub_elem = xml.etree.ElementTree.SubElement(xml_elem,valname)
             xml_sub_elem.attrib["refined"]="false" if ((not(val[1]))|(val[2]!="")) else "true"
             xml_sub_elem.attrib["constr"]="true" if (val[2]!="") else "false"
             xml_sub_elem.text=val[2] if (val[2]!="") else "{}".format(val[0])
             return

        tempfunc4val(xml_elem,'scale',self.scale)
        tempfunc4val(xml_elem,'igsize',self.igsize)
        tempfunc4val(xml_elem,'ext11',self.ext11)
        tempfunc4val(xml_elem,'ext22',self.ext22)
        tempfunc4val(xml_elem,'ext33',self.ext33)
        tempfunc4val(xml_elem,'ext12',self.ext12)
        tempfunc4val(xml_elem,'ext13',self.ext13)
        tempfunc4val(xml_elem,'ext23',self.ext23)
        return  xml_elem





class cmodel_ref(widg_min.cmodel_min):
    """
    minimal model for cwidget experiment
    """
    def __init__(self):
        super(cmodel_ref, self).__init__()
        self.name = "ref"
        self.output = "output.lis"
        self.refin = True
        self.sigmas = True
        self.fdirxml = None
        self.paramref = None
        self.paramrefval = None
        self.paramreferrors = None

    def take_from_core(self,ccore_ref):
        self.output = ccore_ref.output
        self.refin = ccore_ref.refin
        self.sigmas = ccore_ref.sigmas
        self.fdirxml = ccore_ref.fdirxml
        
        if ccore_ref.paramref != None:
            self.paramref = [hh for hh in ccore_ref.paramref]
        else:
            self.paramref = None
        if ccore_ref.paramrefval != None:
            self.paramrefval = [hh for hh in ccore_ref.paramrefval]
        else:
            self.paramrefval = None
        if ccore_ref.paramreferrors != None:
            self.paramreferrors = [hh for hh in ccore_ref.paramreferrors]
        else:
            self.paramreferrors = None

    def load_dict(self,dref):
        if 'output' in dref.keys(): self.output = dref["output"]
        if 'refin' in dref.keys(): self.refin = (dref["refin"] == "true")
        if 'sigmas' in dref.keys(): self.sigmas = (dref["sigmas"] == "true")

    def to_xml_elem(self):
        xml_elem = xml.etree.ElementTree.Element('ref')
        xml_elem.attrib["refin"] = "true" if (self.refin) else "false"
        xml_elem.attrib["output"] = self.output
        xml_elem.attrib["sigmas"] = "true" if (self.sigmas) else "false"
        return xml_elem

