"""
transform inforamtion from cif file to class phase
"""
import os
import sys
import math
import numpy

class cphase(object):
    def __init__(self):
        self.name = None
        self.spgr = None
        self.cella = None
        self.cellb = None
        self.cellc = None
        self.cellalpha = None
        self.cellbeta =  None
        self.cellgamma = None
        self.atom = []
        self.volume = None
        self.rucp = None
    def calc_ucp_recip(self):
        """Calculate reciprocal unit cell parameters\n
        Input parameters:
        ucp -- unit cell parameters [a,b,c,alpha,beta,gamma]"""
        self.calc_ucp_vol()
        vol = self.volume
        a, b, c, alpha, beta, gamma = self.cella, self.cellb, self.cellc, self.cellalpha, self.cellbeta, self.cellgamma
        alpha,beta,gamma=math.radians(alpha),math.radians(beta),math.radians(gamma)
        alphar=numpy.arccos((numpy.cos(beta)*numpy.cos(gamma)-numpy.cos(alpha))/(numpy.sin(beta)*numpy.sin(gamma)))
        betar=numpy.arccos((numpy.cos(gamma)*numpy.cos(alpha)-numpy.cos(beta))/(numpy.sin(alpha)*numpy.sin(gamma)))
        gammar=numpy.arccos((numpy.cos(alpha)*numpy.cos(beta)-numpy.cos(gamma))/(numpy.sin(alpha)*numpy.sin(beta)))
        ar = b*c*numpy.sin(alpha)/vol
        br = c*a*numpy.sin(beta)/vol
        cr = a*b*numpy.sin(gamma)/vol
        rucp = [ar,br,cr,math.degrees(alphar),math.degrees(betar),math.degrees(gammar)]
        self.rucp = rucp
    def calc_ucp_vol(self):
        """ Calculate volume of unit cell\n
        Input parameter: ucp -- list of unit cell parameters [a,b,c,alpha(degree),beta(degree),gamma(degree)]\n
        Output parameters is Volume of unit cell"""
        ucp = [self.cella, self.cellb, self.cellc, self.cellalpha, self.cellbeta, self.cellgamma]
        rad=math.pi/180
        vol=ucp[0]*ucp[1]*ucp[2]*(1-math.cos(ucp[3]*rad)**2-math.cos(ucp[4]*rad)**2-math.cos(ucp[5]*rad)**2+2*math.cos(ucp[3]*rad)*math.cos(ucp[4]*rad)*math.cos(ucp[5]*rad))**0.5
        self.volume = vol
    def trans_u_to_b_beta(self):
        self.calc_ucp_recip()
        rucp = self.rucp
        for atom in self.atom:
            atom.trans_u_iso_to_biso()
            atom.trans_u_aniso_to_beta(rucp)
    def print_phase(self):
        print "unit cell is {:} {:} {:} {:} {:} {:}".format(self.cella, self.cellb, self.cellc, self.cellalpha, self.cellbeta, self.cellgamma)
        print "atoms:"
        print "label type     x       y       z     occ adp_type biso  beta11  beta22  beta33  beta12  beta13  beta23"
        for atom in self.atom:
            atom.print_atom()
        
        

class catom(object):
    def __init__(self):
        self.name = ""
        self.type = ""
        self.coordx = 0.
        self.coordy = 0.
        self.coordz = 0.
        self.bscat = 0.
        self.occ = 1.
        self.lfactor = 1.
        self.biso = 0.
        self.beta11, self.beta22, self.beta33 = 0., 0., 0.
        self.beta12, self.beta13, self.beta23 = 0., 0., 0.
        self.chi11, self.chi22, self.chi33 = 0., 0., 0.
        self.chi12, self.chi13, self.chi23 = 0., 0., 0.
        self.u_iso = 0.
        self.adp_type = None
        self.u11, self.u22, self.u33 = 0., 0., 0.
        self.u12, self.u13, self.u23 = 0., 0., 0.
    def trans_u_iso_to_biso(self):
        #it should be checked
        self.biso = self.u_iso * 8. * math.pi**2
    def trans_u_aniso_to_beta(self, rucp):
        ra, rb, rc = rucp[0], rucp[1], rucp[2]
        lrucp = [ra*ra,rb*rb,rc*rc,ra*rb,ra*rc,rb*rc]
        u_aniso = [self.u11, self.u22, self.u33, self.u12, self.u13, self.u23]
        beta = [2*math.pi**2*rel*u_el for rel, u_el in zip(lrucp, u_aniso)]
        [self.beta11, self.beta22, self.beta33, self.beta12, self.beta13, self.beta23] = beta
        #sbeta=[2*math.pi**2*rel*sUel for rel,sUel in map(None,lrucp,sU)]
    def print_atom(self):
        print "{:5}{:3}{:8.5f}{:8.5f}{:8.5f}{:8.5f} {:5}{:8.5f}{:8.5f}{:8.5f}{:8.5f}{:8.5f}{:8.5f}{:8.5f}".format(self.name, self.type, self.coordx, self.coordy, self.coordz, self.occ, self.adp_type, self.biso, self.beta11, self.beta22, self.beta33, self.beta12, self.beta13, self.beta23)
    


class ccif(object):
    def __init__(self):
        self.help_info = "cif file"
        
        
def read_cif(fname):
    cif = ccif()
    
    fid = open(fname,"r")
    lcontent = fid.readlines()
    fid.close()
    
    lphase = convert_cif_to_phases(lcontent)
    
    for phase in lphase:
        phase.trans_u_to_b_beta()
    
    return lphase

def smart_spleet(str):
    """
    split string like:
    "C C 0.0033 0.0016 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4'"
    in the list like:
    ['C', 'C', '0.0033', '0.0016', 'International Tables Vol C Tables 4.2.6.8 and 6.1.1.4']
    """
    lflag = []
    flag_in = False
    lval, val = [], []
    for hh in str.strip():
        if (hh == " ")&(not flag_in):
            if val != []:
                lval.append("".join(val))
            val = []
        elif (hh == " ")&(flag_in):
            val.append(hh)
        elif hh == "'":
            flag_in = not flag_in
        else:
            val.append(hh)
    if val != []:
        lval.append("".join(val))
    return lval
    
def smart_numb(str):
    """
    convert string of the format '5.7(3)' into number 5.7
    """
    ind = str.find("(")
    if ind == -1:
        try:
            numb = float(str)
        except:
            numb = 0.
    else:
        numb = float(str[:ind])
    return numb
    
def loop_to_dict(lcontent):
    """
    convert list of string presented as loop to dictionary
    """
    dict = {}
    lname_str = [hh for hh in lcontent if hh[0] == "_"]
    lname_tuple = [tuple(hh[1:].strip().split("_")) for hh in lname_str]
    for hh in lname_tuple:
        dict[hh] = []
    ldata_str = [hh for hh in lcontent if hh[0] != "_"]
    for hh in ldata_str:
        lval = smart_spleet(hh)
        if len(lval) != len(lname_tuple):
            print "problem with line: \n{:}".format(hh)
            print len(lval), len(lname_tuple)
        for name, val in zip(lname_tuple, lval):
            dict[name].append(val)
    return dict

def lstr_to_dict(lcontent):
    """
    convert list of string presented as parameter with one value to dictionary
    """
    dict = {}
    lflag = [True if hh.startswith("_") else False for hh in lcontent[1:]]
    lflag.append(False)
    val = []
    for flag, hh in zip(lflag, lcontent):
        if hh.startswith("_"):
            if val != []:
                dict[name_tuple] = "\n".join(val)
                val = []
            lhelp = smart_spleet(hh)
            name_str = lhelp[0]
            name_tuple = tuple(name_str[1:].strip().split("_"))
            if flag:
                if len(lhelp) >= 1:
                    dict[name_tuple] = " ".join(lhelp[1:])
                else:
                    dict[name_tuple] = ""
        else:
            if (not hh.startswith(";")):
                val.append(hh.strip())
    return dict

def serch_in_dict(dict, l_to_search, only_it = False):
    """
    keys of dictionary is tuples.
    If all keay words from the list l_to_search are in name of key then value of such variable is added to lval
    """
    ltuple, lval = [], []
    lname_tuple = dict.keys()
    for name_tuple in lname_tuple:
        name_tuple_upper = tuple([hh.upper() for hh in name_tuple])
        cond = all([(elem.upper() in name_tuple_upper) for elem in l_to_search])
        if only_it:
            cond &= (len(name_tuple)==len(l_to_search))
        if cond:
            ltuple.append(name_tuple)
            lval.append(dict[name_tuple])
    return ltuple, lval

    
def dict_to_phase(dict):
    """
    Transfer dictionary to class cphase
    """
    phase = cphase()
    
    def cfunc(l_to_search):
        ltuple, lval = serch_in_dict(dict, l_to_search)
        if len(lval) > 1:
            print "Several solution is found. The first one is taken.\nList of the solutions is"
            for tup, val in zip(ltuple, lval):
                print tup, val
            val = lval[0]
        elif len(lval) == 0:
            val = None
            print "The parameters '"+" ".join(l_to_search)+"' are not found."
        else:
            val = lval[0]
        return val
    
    phase.cella = smart_numb(cfunc(["length","a"]))
    phase.cellb = smart_numb(cfunc(["length","b"]))
    phase.cellc = smart_numb(cfunc(["length","c"]))
    phase.cellalpha = smart_numb(cfunc(["angle","alpha"]))
    phase.cellbeta = smart_numb(cfunc(["angle","beta"]))
    phase.cellgamma = smart_numb(cfunc(["angle","gamma"]))
    phase.spgr = cfunc(["symmetry","space","group","name","H-M"])
    return phase

def dict_to_atoms(dict):
    """
    Transfer dictionary to class catom
    """
    latom = []
    def cfunc(l_to_search, only_it = False, numb = 1):
        ltuple, lval = serch_in_dict(dict, l_to_search, only_it)
        if len(lval) == 0:
            print "The parameters '"+" ".join(l_to_search)+"' are not found."
            val = numb*[None]
        else:
            val = lval[0]
        return val
    lname = cfunc(["atom","site","label"], True)
    numb = len(lname)
    ltype = cfunc(["atom","site","type","symbol"], True, numb)
    lcoordx = cfunc(["atom","site","fract","x"], True, numb)
    lcoordy = cfunc(["atom","site","fract","y"], True, numb)
    lcoordz = cfunc(["atom","site","fract","z"], True, numb)
    locc = cfunc(["atom","site","occupancy"], True, numb)
    lu_iso = cfunc(["atom","site","U","iso","or","equiv"], True, numb)
    ladp_type = cfunc(["atom","site","adp","type"],True, numb)
    for name, type, coordx, coordy, coordz, occ, u_iso, adp_type in zip(lname, ltype, lcoordx, lcoordy, lcoordz, locc, lu_iso, ladp_type):
        atom = catom()
        if name!=None: atom.name = name.strip()
        if type!=None: atom.type = type.strip()
        if coordx!=None: atom.coordx = smart_numb(coordx)
        if coordy!=None: atom.coordy = smart_numb(coordy)
        if coordz!=None: atom.coordz = smart_numb(coordz)
        if occ!=None: atom.occ = smart_numb(occ)
        if u_iso!=None: atom.u_iso = smart_numb(u_iso)
        if adp_type!=None: atom.adp_type = adp_type.strip()
        latom.append(atom)
    return latom

def add_dict_to_atoms(latom,dict_to_add,atom_label,dict_tuple):
    """
    only for U_11, ..., U_23
    """

    def cfunc(l_to_search, only_it = False):
        ltuple, lval = serch_in_dict(dict_to_add, l_to_search, only_it)
        if len(lval) == 0:
            print "The parameters '"+" ".join(l_to_search)+"' are not found."
        return lval[0]
    
    llabel = cfunc(list(dict_tuple), True)
    lu_11 = cfunc(["atom", "site", "aniso", "U", "11"], True)
    lu_22 = cfunc(["atom", "site", "aniso", "U", "22"], True)
    lu_33 = cfunc(["atom", "site", "aniso", "U", "33"], True)
    lu_12 = cfunc(["atom", "site", "aniso", "U", "12"], True)
    lu_13 = cfunc(["atom", "site", "aniso", "U", "13"], True)
    lu_23 = cfunc(["atom", "site", "aniso", "U", "23"], True)
    for atom in latom:
        for label, u_11, u_22, u_33, u_12, u_13, u_23 in zip(llabel, lu_11, lu_22, lu_33, lu_12, lu_13, lu_23):
            if atom.__dict__[atom_label] == label:
                atom.u11 = smart_numb(u_11)
                atom.u22 = smart_numb(u_22)
                atom.u33 = smart_numb(u_33)
                atom.u12 = smart_numb(u_12)
                atom.u13 = smart_numb(u_13)
                atom.u23 = smart_numb(u_23)
                break
    
def convert_cif_to_phases(lcontent):
    lphase= []
    flag = True
    lflag = [False if (hh.strip().startswith("#")|(hh.strip() == "")) else True for hh in lcontent]
    #lflag = [False if ((hh[0] == "#")|((hh[0]==" ")|(hh[0]=="\n"))) 
    #         else True for hh in lcontent]
    lcontent_1 = [hh.strip() for hh, flag in zip(lcontent,lflag) if flag]
    lcontent_1 = [hh[:hh.find("#")] if hh.find("#") != -1 else hh for hh in lcontent_1]
    sword = "data_"
    lnumb_b = [ihh for ihh, hh in enumerate(lcontent_1) if hh.startswith(sword)]
    lnumb_e = []
    if len(lnumb_b) > 1:
        lnumb_e = lnumb_b[1:]
    lnumb_e.append(len(lcontent_1))
    for numb_b, numb_e in zip(lnumb_b, lnumb_e):
        data_name = lcontent_1[numb_b].strip().split("_")[1]
        lcontent_2 = lcontent_1[(numb_b+1):numb_e]
        lflag = []
        flag = True
        sword = "loop_"
        for hh in lcontent_2:
            cond_1 = hh.startswith(sword) 
            if cond_1:
                flag = False
                lflag.append(flag)
            elif ((hh[0] != "_")&(not cond_1)):
                flag = True
                lflag.append(lflag[-1])
            elif (hh[0] == "_"):
                lflag.append(flag)
            else:
                print "Problem with line:\n{:}".format(hh)
        #for flag, hh in zip(lflag, lcontent_2):
        #    print flag, "          ", hh[:-1]
        #print ""
        lcontent_2_single_param = [hh for flag, hh in zip(lflag, lcontent_2) if flag]
        lcontent_2_loops = [hh for flag, hh in zip(lflag, lcontent_2) if (not flag)]
        
        dict_one_val = lstr_to_dict(lcontent_2_single_param)
        phase = dict_to_phase(dict_one_val)
        phase.name = data_name

        sword = "loop_"
        lnumb_b_2 = [ihh for ihh, hh in enumerate(lcontent_2_loops) if hh.startswith(sword)]
        lnumb_e_2 = []
        if len(lnumb_b_2) > 1:
            lnumb_e_2 = lnumb_b_2[1:]
        lnumb_e_2.append(len(lcontent_2_loops))
        latom = []
        for numb_b_2, numb_e_2 in zip(lnumb_b_2, lnumb_e_2):
            lcontent_2_one_loop = lcontent_2_loops[(numb_b_2+1):numb_e_2]
            dict_loop = loop_to_dict(lcontent_2_one_loop)
            l_to_search = ["atom","site","fract","x"]
            ltuple, lval = serch_in_dict(dict_loop, l_to_search)
            if lval!= []:
                latom = dict_to_atoms(dict_loop)
        #add anisotropic dispacement parameters
        for numb_b_2, numb_e_2 in zip(lnumb_b_2, lnumb_e_2):
            lcontent_2_one_loop = lcontent_2_loops[(numb_b_2+1):numb_e_2]
            dict_loop = loop_to_dict(lcontent_2_one_loop)
            l_to_search = ["atom", "site", "aniso", "U", "11"]
            ltuple, lval = serch_in_dict(dict_loop, l_to_search)
            if lval!= []:
                add_dict_to_atoms(latom,dict_loop,"name",("atom", "site", "aniso", "label"))
        phase.atom = latom
        lphase.append(phase)
    return lphase
        
def main(larg):
    if len(larg) > 1:
        fname = larg[1]
    else:
        fname = raw_input("What is the name of the 'cif' file?\n")
    lphase = read_cif(fname)
    print "Readed phases:"
    for phase in lphase:
        print ""
        phase.trans_u_to_b_beta()
        phase.print_phase()

if __name__ == "__main__":
    main(sys.argv)