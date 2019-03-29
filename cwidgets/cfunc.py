#define number of functions for general algebra
import math
import cmath
import numpy
import numpy.linalg
from PyQt5 import QtWidgets, QtGui
from PyQt5.QtCore import Qt

def float_to_string(numb,digits = 2):
    """
    transform any number to string:
        12278790 = 1.23*10^7
        0.000060 = 6.00*10^-5
    """
    if numb == 0.:
        snumb = "0.00"
        return snumb
    elif numb > 0:
        numb_power = int(math.log10(numb))
    else:        
        numb_power = -int(-1.*math.log10(numb))
    numb_main = round(numb*10^(-numb_power),digits)
    snumb = ""
    return snumb

def save_gui_to_val(gui_u):
    lhelp = gui_u[1].statusTip().split(":")
    if len(lhelp) == 2:
        sconstr = "{:}".format(lhelp[1])
    else:
        sconstr = ""
    refined = (gui_u[0].checkState() != 0)
    val = float(gui_u[1].text())
    return [val,refined,sconstr]


def make_gui_val(label, val):
    
    
    layout=QtWidgets.QHBoxLayout()
    if isinstance(val,list):
        gui_cb = QtWidgets.QCheckBox(label)
        gui_cb.setCheckState(2*val[1])
        gui_le = QtWidgets.QLineEdit("{0[0]:.5f}".format(val))
        if (val[2]!=""):
            gui_le.setReadOnly(True)
            gui_le.setToolTip("constraint: {}".format(val[2]))
            gui_le.setStatusTip("restriction: {:}".format(val[2]))
            gui_le.setStyleSheet("color: rgb(156, 156, 156);")
        else:
            gui_le.setReadOnly(False)
            gui_le.setStyleSheet("color: rgb(0, 0, 0);")
            #self.ui.lineEdit_5.setStyleSheet("color: rgb(28, 43, 255);")

            #gui_le.setToolTip("no constraint on the parameter")
            gui_le.setStatusTip("no restrictions (press right mouse button to set up it)")
        gui_le.setContextMenuPolicy(Qt.CustomContextMenu)
        gui_le.customContextMenuRequested.connect(lambda x: cont_menu_le(gui_le,val))
    else:
        gui_cb = QtWidgets.QLabel("{}".format(label))
        gui_le = QtWidgets.QLineEdit("{}".format(val))
    gui_le.setFrame(False)
    gui_le.setAlignment(Qt.AlignRight)
    gui_le.setTextMargins(0,0,0,0)
    layout.addWidget(gui_cb)
    layout.addStretch(1)
    layout.addWidget(gui_le)
    #layout.setSpacing(0)
    return [gui_cb,gui_le],layout


def cont_menu_le(obj1,val):
    cursor = QtGui.QCursor()
    cmenu = QtWidgets.QMenu()
    act_constr_def = cmenu.addAction("Set as constraint")
    #act_constr_def.setStatusTip("press to set constraint on the parameter")
    act_constr_del = cmenu.addAction("Delete constraint")
    action = cmenu.exec_(cursor.pos())
    if action == act_constr_del:
        obj1.setStatusTip("no restrictions (press right mouse button to set up it)")
        #val[2] = ""
        obj1.setReadOnly(False)
    elif action == act_constr_def:
        text,ok = QtWidgets.QInputDialog.getText(obj1,'Input Dialog',
            'Enter constraint:')
        if ok:
            obj1.setStatusTip("restriction: {:}".format(str(text)))
            #val[2]=str(text)
            #val[2]="{}".format(obj1.text())
            obj1.setReadOnly(True)
        #obj1.setText("{:.5f}".format(val[0]))

def readexpdata(finp):
    ddata={}
    fid=open(finp,'r')
    lcontentH=fid.readlines()
    fid.close()
    lparam=[line[1:].strip() for line in lcontentH if line.startswith('#')]
    if (len(lparam)>1):
        for line in lparam:
            lhelp=splitlinewminuses(line)
            if (len(lhelp)>2):
                ddata[lhelp[0]]=lhelp[1:]
            elif (len(lhelp)==2):
                ddata[lhelp[0]]=lhelp[1]
            else:
                print "Mistake in experimental file '{}' in line:\n{}".format(finp,line)
                print "The program is stopped."
                quit()
    lnames=lparam[-1].split()
    for name in lnames:
        ddata[name]=[]
    lcontent=[line for line in lcontentH if line[0]!='#']
    for line in lcontent:
        for name,val in zip(lnames,splitlinewminuses(line)):
            ddata[name].append(val)
    return ddata

def read_2dmatrices_np(finp):
    """
    Read 2d dimension data separated by empty lines
    """
    ll_int, l_ang1, l_ang2 = [], [], [] 

    fid = open(finp, 'r')
    lcont = fid.readlines()
    fid.close()
    
    lcont_1 = [hh[1:] for hh in lcont if (hh.startswith("#"))]
    dinput = {}
    for hh in lcont_1:
        lhelp = hh.split()
        if lhelp[0] == "wavelength":
            dinput["wavelength"] = float(lhelp[1])
        elif lhelp[0] == "field":
            dinput["field"] = [float(hh2) for hh2 in lhelp[1:]]
            
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
    return ll_int, l_ang1, l_ang2, dinput

def save_2dmatrices_np(finp, np_phi, np_tth, np_int_e_u, np_int_e_su, np_int_e_d, np_int_e_sd, np_int_e_diff, np_int_e_sdiff):
    line = "       {:7}   ".format(len(np_phi)) + " ".join(["{:10.2f}".format(hh) for hh in np_tth])
    lsout_u = [line]
    lsout_su = [line]
    lsout_d = [line]
    lsout_sd = [line]
    lsout_diff = [line]
    lsout_sdiff = [line]
    
    for lint_u_prof, lint_su_prof, lint_d_prof, lint_sd_prof, lint_diff_prof, lint_sdiff_prof, phi in zip(
        np_int_e_u, np_int_e_su, np_int_e_d, np_int_e_sd, np_int_e_diff, np_int_e_sdiff, np_phi):
                        
        line_u = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_u_prof])
        line_su = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_su_prof])
        line_d = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_d_prof])
        line_sd = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_sd_prof])
        line_diff = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_diff_prof])
        line_sdiff = "    {:10.5f}    ".format(phi) + " ".join(["{:10.2f}".format(hh) for hh in lint_sdiff_prof])
        lsout_u.append(line_u)
        lsout_su.append(line_su)
        lsout_d.append(line_d)
        lsout_sd.append(line_sd)
        lsout_diff.append(line_diff)
        lsout_sdiff.append(line_sdiff)
                    

    fid = open(finp, "w")
    fid.write("\n".join(lsout_u)+3*"\n"+"\n".join(lsout_su)+3*"\n"+
              "\n".join(lsout_d)+3*"\n"+"\n".join(lsout_sd)+3*"\n"+
              "\n".join(lsout_diff)+3*"\n"+"\n".join(lsout_sdiff)+3*"\n")
    fid.close()
    return
    

def read_el_cards(fitables):
    """
    reading information about space grooupe from file fitables to list of cards ldcard
    Info in file fitables:
    
    1 P1               Triclinic
    choice: 1
    centr: false
    pcentr: 0, 0, 0
    symmetry: X,Y,Z
    
    2 P-1              Triclinic
    ...
    """
    fid = open(fitables, "r")
    lcontent = fid.readlines()
    fid.close()

    lcontent = [hh.strip() for hh in lcontent if hh.strip() != ""]
    ldcard = []
    dcard = None
    for hh in lcontent:
        lhelp = hh.split()
        if lhelp[0].isdigit():
            if dcard != None:
                ldcard.append(dcard)
            dcard = {"number":lhelp[0], "name": lhelp[1], "syngony": lhelp[2]}
        else:
            lhelp = hh.split(":")
            if (lhelp[0].strip() in dcard.keys()):
                dcard[lhelp[0].strip()].append(lhelp[1].strip())
            else:
                dcard[lhelp[0].strip()] = [lhelp[1].strip()]
    return ldcard 

def calcsthovl(hkl,ucp):
    """ calculate sin(theta)/lambda from Miller indices h,k,l
    ucp is unit cell parameters listed in a list or 1 dimensional ndarray as [a(Angst),b(Angst),c(Angst),alpha(degree),beta(degree),gamma(degree)]"""
    A=( 1 - math.cos(math.radians(ucp[3]))**2 - math.cos(math.radians(ucp[4]))**2 - math.cos(math.radians(ucp[5]))**2 + 2*math.cos(math.radians(ucp[3]))*math.cos(math.radians(ucp[4]))*math.cos(math.radians(ucp[5])))
    B1 = ((math.sin(math.radians(ucp[3]))*hkl[0]*1./ucp[0])**2+(math.sin(math.radians(ucp[4]))*hkl[1]*1./ucp[1])**2+(math.sin(math.radians(ucp[5]))*hkl[2]*1./ucp[2])**2)
    B2 = 2.*(hkl[1]*hkl[2]*math.cos(math.radians(ucp[3])))/(ucp[1]*ucp[2])+2.*(hkl[0]*hkl[2]*math.cos(math.radians(ucp[4])))/(ucp[0]*ucp[2])+2.*(hkl[0]*hkl[1]*math.cos(math.radians(ucp[5])))/(ucp[0]*ucp[1])
    #it should be checked, I am not sure
    B = B1-B2
    invd=(B*1./A)**0.5
    return (invd>0)*0.5*invd**((invd>0))


def calc_u_b_ucp_from_ub(mub):
    """
    calculate matrix U and matrix B from UB matrix
    UB = U * B
    UB = [[ub_11, ub_12, ub_13], 
          [ub_21, ub_22, ub_23],
          [ub_31, ub_32, ub_33]]
    
    """
    mu = [[0,0,0],[0,0,0],[0,0,0]]
    mb = [[0,0,0],[0,0,0],[0,0,0]]
    [ub_11, ub_12, ub_13] = mub[0]
    [ub_21, ub_22, ub_23] = mub[1]
    [ub_31, ub_32, ub_33] = mub[2]
    ia = [ub_11, ub_21, ub_31]
    n_ia = sum([hh1**2 for hh1 in ia])**0.5
    ib = [ub_12, ub_22, ub_32]
    n_ib = sum([hh1**2 for hh1 in ib])**0.5
    ic = [ub_13, ub_23, ub_33]
    n_ic = sum([hh1**2 for hh1 in ic])**0.5
    igamma = math.acos(sum([hh1*hh2 for hh1, hh2 in zip(ia, ib)])/(n_ia*n_ib))
    ibeta = math.acos(sum([hh1*hh2 for hh1, hh2 in zip(ia, ic)])/(n_ia*n_ic))
    ialpha = math.acos(sum([hh1*hh2 for hh1, hh2 in zip(ib, ic)])/(n_ib*n_ic))
    ucp = calcrucp([n_ia, n_ib, n_ic, math.degrees(ialpha), math.degrees(ibeta), math.degrees(igamma)])
    c, alpha = ucp[2], math.radians(ucp[3])
    mb = [[n_ia, n_ib*math.cos(igamma), n_ic*math.cos(ibeta)],
         [   0, n_ib*math.sin(igamma), -n_ic*math.sin(ibeta)*math.cos(alpha)],
         [   0,                    0, 1/c]]
    
    mib = numpy.linalg.inv(mb)
    mu = [[float(hh2) for hh2 in hh1] for hh1 in multMAT(mub,mib)]
    return mu, mb, ucp
    
def calcrucp(ucp):
    """Calculate reciprocal unit cell parameters\n
    Input parameters:
    ucp -- unit cell parameters [a,b,c,alpha,beta,gamma]"""
    vol=calcVolOfUC(ucp)
    [a,b,c,alpha,beta,gamma]=ucp
    alpha,beta,gamma=math.radians(alpha),math.radians(beta),math.radians(gamma)
    alphar=numpy.arccos((numpy.cos(beta)*numpy.cos(gamma)-numpy.cos(alpha))/(numpy.sin(beta)*numpy.sin(gamma)))
    betar=numpy.arccos((numpy.cos(gamma)*numpy.cos(alpha)-numpy.cos(beta))/(numpy.sin(alpha)*numpy.sin(gamma)))
    gammar=numpy.arccos((numpy.cos(alpha)*numpy.cos(beta)-numpy.cos(gamma))/(numpy.sin(alpha)*numpy.sin(beta)))
    ar=b*c*numpy.sin(alpha)/vol
    br=c*a*numpy.sin(beta)/vol
    cr=a*b*numpy.sin(gamma)/vol
    return [ar,br,cr,math.degrees(alphar),math.degrees(betar),math.degrees(gammar)]

def calcmT(kloc):
    """define rotation matrix to have new z axis along kloc
    Rotation matrix is defined by Euler angles
    """
    be=math.acos(kloc[2])
    sb=math.sin(be)
    if (sb!=0.):
        sa1=kloc[0]*1./sb
        ca2=-1*kloc[1]*1./sb
        if (sa1>1):  sa1=1.
        if (sa1<-1): sa1=-1.
        if (ca2>1):  ca2=1.
        if (ca2<-1): ca2=-1.
        al1=math.asin(sa1)
        al2=math.acos(ca2)
        if (abs(al2-al1)<0.00001):
            al=al1
        elif (sa1>0.):
            al=al2
        elif (sa1<=0.):
            al=2.*math.pi-al2
    else:
        al=0.
    ga=0.
    ca,cb,cg=math.cos(al),math.cos(be),math.cos(ga)
    sa,sb,sg=math.sin(al),math.sin(be),math.sin(ga)
    mT=[[ca*cg-sa*cb*sg, -ca*sg-sa*cb*cg,  sa*sb],
        [sa*cg+ca*cb*sg, -sa*sg+ca*cb*cg, -ca*sb],
        [         sb*sg,           sb*cg,     cb]]
    if (((sa*sb-kloc[0])**2+(-ca*sb-kloc[1])**2+(cb-kloc[2])**2)>0.0001):
        if tuple(kloc)==(0., 0., 0.):
            return mT
        print "Mistake with kloc: {:.4f} {:.4f} {:.4f}".format(kloc[0],kloc[1],kloc[2])
        print sa1,ca2,math.degrees(al)/180.
        print math.sin(al),math.cos(al)
        print "Program is stopped"
        quit()
    return mT


def calc_np_shape_hkl(np_ttheta, ttheta_hkl, zeroshift, UVWIgxy, zshift_a, p_asym):
    if ttheta_hkl == 0.:
        return 0.*np_ttheta
    
    U, V, W, Ig, x, y = UVWIgxy[0], UVWIgxy[1], UVWIgxy[2], UVWIgxy[3], UVWIgxy[4], UVWIgxy[5]
    #function for assymetry
    func_fa = lambda z: 2*z*numpy.exp(-z**2)
    func_fb = lambda z: 2*(2*z**2-3)* func_fa(z)

    p1, p2, p3, p4 = p_asym[0], p_asym[1], p_asym[2], p_asym[3]

    def calc_assym(z,  tth_hkl):
        return 1.+(p1*func_fa(z)+p2*func_fb(z))*1./numpy.tanh(0.5*tth_hkl)+(p3*func_fa(z)+p4*func_fb(z))*1./numpy.tanh(tth_hkl)

    gauss_p = lambda ag, bg, np_ttheta : ag*numpy.exp(-bg*np_ttheta**2)
    lor_p = lambda al, bl, np_ttheta : al*1./(1.+bl*np_ttheta**2)


    np_zshift = zeroshift + zshift_a*numpy.cos(0.5*np_ttheta)
    np_theta_rad = (0.5*math.pi/180)*np_ttheta
    
    np_hg2 = abs(U*numpy.tan(np_theta_rad)**2+V*numpy.tan(np_theta_rad)+W+Ig*1.0/((numpy.cos(np_theta_rad))**2))
    np_hl = x*numpy.tan(np_theta_rad) + y*1./numpy.cos(np_theta_rad)
    np_hg = np_hg2**0.5

    np_fwhm = (np_hg**5 + 2.69269*(np_hg**4)*(np_hl) + 
               2.42843*(np_hg**3)*(np_hl**2) + 4.47163*(np_hg**2)*(np_hl**3) + 
               0.07842*np_hg*(np_hl**4) + np_hl**5)**0.2
               
    np_hh = np_hl*1./np_fwhm
    np_eta = 1.36603*np_hh - 0.47719*np_hh**2 + 0.11116*np_hh**3
    np_ag = (2./np_fwhm)*(math.log(2.)/math.pi)**0.5
    np_bg = 4*math.log(2)/(np_fwhm**2)
    np_al = 2./(math.pi*np_fwhm)
    np_bl = 4./(np_fwhm**2)
    np_LorFact=(1.0/(numpy.sin(np_theta_rad)*numpy.sin(2.*np_theta_rad)))

    np_shape_hkl = numpy.where(abs(np_ttheta-(ttheta_hkl+np_zshift))<10., 
       np_LorFact*((np_eta*lor_p(np_al, np_bl, np_ttheta - (ttheta_hkl+ np_zshift)) + 
       (1.-np_eta)*gauss_p(np_ag, np_bg, np_ttheta-(ttheta_hkl+np_zshift)))*
        calc_assym((np_ttheta-(ttheta_hkl+np_zshift))*1./np_fwhm, ttheta_hkl)),
         0.0)
    return np_shape_hkl

def calcprofile(ltthetahkl, lInthkl, lttheta, zeroshift, UVWIgxy, zshift_a, p_asym):
    """calculate profile for reflections
    """
    U, V, W, Ig, x, y = UVWIgxy[0], UVWIgxy[1], UVWIgxy[2], UVWIgxy[3], UVWIgxy[4], UVWIgxy[5]
    lInt=[]
    
    #function for assymetry
    func_fa = lambda z: 2*z*math.exp(-z**2)
    func_fb = lambda z: 2*(2*z**2-3)* func_fa(z)
    
    p1, p2, p3, p4 = p_asym[0], p_asym[1], p_asym[2], p_asym[3]

    def calc_assym(z,  tth_hkl):
        if tth_hkl !=0. :
            res = 1.+(p1*func_fa(z)+p2*func_fb(z))*1./math.tanh(0.5*tth_hkl)+(p3*func_fa(z)+p4*func_fb(z))*1./math.tanh(tth_hkl)
        else:
            res = 1. 
        return res

    gauss_p = lambda ag, bg, ttheta : ag*math.exp(-bg*ttheta**2)
    lor_p = lambda al, bl, ttheta : al*1./(1.+bl*ttheta**2)
    
    for ttheta in lttheta:
        zshift = zeroshift+zshift_a*math.cos(0.5*ttheta)
        rtheta = math.radians(0.5*ttheta)
        hg2 = abs(U*math.tan(rtheta)**2+V*math.tan(rtheta)+W+Ig*1.0/((math.cos(rtheta))**2))
        hl = x*math.tan(rtheta)+y*1./math.cos(rtheta)
        hg = hg2**0.5

        fwhm = (hg**5 + 2.69269*(hg**4)*(hl) + 2.42843*(hg**3)*(hl**2) + 4.47163*(hg**2)*(hl**3) + 0.07842*hg*(hl**4) + hl**5)**0.2
        hh = hl*1./fwhm
        eta = 1.36603*hh - 0.47719*hh**2 + 0.11116*hh**3
        ag = (2./fwhm)*(math.log(2.)/math.pi)**0.5
        bg = 4*math.log(2)/(fwhm**2)
        al = 2./(math.pi*fwhm)
        bl = 4./(fwhm**2)
        LorFact=(1.0/(math.sin(rtheta)*math.sin(2*rtheta)))

        lres1=[0.]
        #ttheta, tthetahkl, zeroshift are in degrees
        lres1.extend([Inthkl*(eta*lor_p(al,bl,ttheta-(tthetahkl+zshift))+(1.-eta)*gauss_p(ag,bg,ttheta-(tthetahkl+zshift)))*calc_assym((ttheta-(tthetahkl+zshift))*1./fwhm,tthetahkl) if abs(ttheta-(tthetahkl+zshift))<10. else 0. for tthetahkl,Inthkl in zip(ltthetahkl,lInthkl)])
        lInt.append(LorFact*sum(lres1))
    return lInt

def chiLOC(mchi,ucp):
    """
    representation of chi in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
    expressions are taken from international tables
    X = M x, x = iM X
    xT*CHI*x = XT iMT CHI iM X

    output chiLOC = iMT CHI iM
    """
    #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
    #[a,b,c,alpha,beta,gamma]=ucp
    rucp=calcrucp(ucp)
    [b1,b2,b3]=rucp[:3]
    beta2,beta3=math.radians(rucp[4]),math.radians(rucp[5])
    a3,alpha1=ucp[2],math.radians(ucp[3])
    x1 = b1
    x2 = b2*math.sin(beta3)
    x3 = 1./a3
    x4 = b2*math.cos(beta3)
    x5 = b3*math.cos(beta2)
    x6 = -b3*math.sin(beta2)*math.cos(alpha1)
    #B=[[x1,x4,x5],
    #   [0.,x2,x6],
    #   [0.,0.,x3]]
    #it shuld be checked
    #iB=numpy.linalg.inv(B)
    y1 = 1./x1
    y2 = 1./x2
    y3 = 1./x3
    y4 = -1*x4*1./(x1*x2)
    y6 = -1*x6*1./(x2*x3)
    y5 = (x4*x6-x2*x5)*1./(x1*x2*x3)
    iB=[[y1,y4,y5],
        [0.,y2,y6],
        [0.,0.,y3]]
    iBn=[[iB[0][0]*b1,iB[0][1]*b2,iB[0][2]*b3],
         [iB[1][0]*b1,iB[1][1]*b2,iB[1][2]*b3],
         [iB[2][0]*b1,iB[2][1]*b2,iB[2][2]*b3]]
    iBTn=transposeM(iBn)
    chiLOC=multMAT(multMAT(iBTn,mchi),iBn)
    #alphar,betar,gammar=math.radians(alpha),math.radians(beta),math.radians(gamma)
    #ca,cb,cg=math.cos(alphar),math.cos(betar),math.cos(gammar)
    #sa,sb,sg=math.sin(alphar),math.sin(betar),math.sin(gammar)
    #ta=sa/ca
    #casq,cbsq,cgsq=ca**2,cb**2,cg**2
    #sasq,sbsq,sgsq=sa**2,sb**2,sg**2
    #phi=(1-casq-cbsq-cgsq+2*ca*cb*cg)**0.5
    #M=[[a*phi/sa,          0.,0.],
    #   [a*(cg-ca*cb)/sa, b*sa,0.],
    #   [a*cb,            b*ca, c]]
    #modX=((a*phi/sa)**2+(a*(cg-ca*cb)/sa)**2+(a*cb)**2)**0.5
    #modY=b
    #modZ=c
    #iM=[[sa/(a*phi),                   0.,  0.],
    #    [(ca*cb-cg)/(b*phi*sa), 1./(b*sa),  0.],
    #    [(ca*cg-cb)/(c*phi*sa),-1./(c*ta),1./c]]
    #I am not sure about it
    #iMn=[[           sa/(a*phi)*modX,        0.*modY,  0.*modZ],
    #     [(ca*cb-cg)/(b*phi*sa)*modX, 1./(b*sa)*modY,  0.*modZ],
    #     [(ca*cg-cb)/(c*phi*sa)*modX,-1./(c*ta)*modY,1./c*modZ]]
    #hh1=len(iM)
    #hh2=len(iM[0])
    #iMT=[[ctemp*iM[ih1][ih2] for ih1 in range(hh1)]for ih2 in range(hh2)]
    #iMTn=[[iMn[ih1][ih2] for ih1 in range(hh1)]for ih2 in range(hh2)]
    #not sure
    #chiLOC=multMAT(multMAT(iMTn,mchi),iMn)
    #chiLOC=multMAT(multMAT(M,mchi),iM)
    return chiLOC

def magLOC(mag,ucp):
    """
    representation of mag in crystallographic coordinate system defined as x||a*, z||c, y= [z x] (right handed)
    expressions are taken from international tables
    X = M x, x = iM X
    xT*CHI*x = XT iMT CHI iM X

    output chiLOC = iMT CHI iM
    """
    #mchi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
    #[a,b,c,alpha,beta,gamma]=ucp
    rucp=calcrucp(ucp)
    [b1,b2,b3]=rucp[:3]
    beta2,beta3=math.radians(rucp[4]),math.radians(rucp[5])
    a3,alpha1=ucp[2],math.radians(ucp[3])
    x1 = b1
    x2 = b2*math.sin(beta3)
    x3 = 1./a3
    x4 = b2*math.cos(beta3)
    x5 = b3*math.cos(beta2)
    x6 = -b3*math.sin(beta2)*math.cos(alpha1)
    #B=[[x1,x4,x5],
    #   [0.,x2,x6],
    #   [0.,0.,x3]]
    #it shuld be checked
    #iB=numpy.linalg.inv(B)
    y1 = 1./x1
    y2 = 1./x2
    y3 = 1./x3
    y4 = -1*x4*1./(x1*x2)
    y6 = -1*x6*1./(x2*x3)
    y5 = (x4*x6-x2*x5)*1./(x1*x2*x3)
    iB=[[y1,y4,y5],
        [0.,y2,y6],
        [0.,0.,y3]]
    iBn=[[iB[0][0]*b1,iB[0][1]*b2,iB[0][2]*b3],
         [iB[1][0]*b1,iB[1][1]*b2,iB[1][2]*b3],
         [iB[2][0]*b1,iB[2][1]*b2,iB[2][2]*b3]]
    iBTn=transposeM(iBn)
    #totaly not sure !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #it should be necessary checked !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    mag_loc = multMv(iBTn,mag)
    return mag_loc



def chirot(chi,elsymm):
    """
    calculate R*chi*RT
    rotation of chi by element of symmetry
    """
    [b1,r11,r12,r13,b2,r21,r22,r23,b3,r31,r32,r33]=elsymm
    mR= [[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]]
    mRT=[[r11,r21,r31],[r12,r22,r32],[r13,r23,r33]]
    mChi=[[chi[0],chi[3],chi[4]],[chi[3],chi[1],chi[5]],[chi[4],chi[5],chi[2]]]
    mChir=multMAT(multMAT(mR,mChi),mRT)
    return mChir

def magrot(mag,elsymm):
    """
    calculate R*mag
    rotation of mag by element of symmetry
    """
    [b1,r11,r12,r13,b2,r21,r22,r23,b3,r31,r32,r33]=elsymm
    mR= [[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]]
    #mRT=[[r11,r21,r31],[r12,r22,r32],[r13,r23,r33]]
    magr=multMv(mR,mag)
    return magr

def multMAT(A,B):
    """
    multiplication of two matrix
    """
    BT=transposeM(B)
    C=[[sum([hh3*hh4 for hh3,hh4 in zip(hh1,hh2)]) for hh2 in BT] for hh1 in A]
    return C

def multMv(M,v):
    """multiplication of matrix on vector
    """
    vres=[sum([hh1*hh2 for hh1,hh2 in zip(hM,v)]) for hM in M]
    return vres

def transposeM(M):
    """matrix transpose
    """
    nhh1=len(M[0])
    MT=[[hh1[ihh2] for hh1 in M] for ihh2 in range(nhh1)]
    return MT

def vecprod(a,b):
    [a0,a1,a2]=a
    [b0,b1,b2]=b
    c=[a2*b1-a1*b2, a0*b2-a2*b0, a1*b0-a0*b1]
    return c

def unitv(v):
    normv=([hh**2 for hh in v])**0.5
    uv=[hh*1./normv for hh in v]
    return uv,normv

def els4pos(dsymm,xyz):
    """
    give the lelements of symmetry which transfer atom to the same atom
    """
    lelsymm = dsymm["elsymm"]
    lorig = dsymm["orig"]
    centr = dsymm["centr"]
    pcentr = dsymm["pcentr"]
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


def calcVolOfUC(ucp):
    """ Calculate volume of unit cell\n
    Input parameter: ucp -- list of unit cell parameters [a,b,c,alpha(degree),beta(degree),gamma(degree)]\n
    Output parameters is Volume of unit cell"""
    rad=math.pi/180
    vol=ucp[0]*ucp[1]*ucp[2]*(1-math.cos(ucp[3]*rad)**2-math.cos(ucp[4]*rad)**2-math.cos(ucp[5]*rad)**2+2*math.cos(ucp[3]*rad)*math.cos(ucp[4]*rad)*math.cos(ucp[5]*rad))**0.5
    return vol

def tr2elsymm(str1):
    """transform string to element of symmetry: (x,y,-z) -> 0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
    """
    str2="".join(str1.split(" "))
    lhelp1,lhelp2,lhelp3=[],[],[]
    lhelp1=[hh for hh in str2.split('(') if hh!=""]
    [lhelp2.extend(hh.split(')')) for hh in lhelp1 if hh!=""]
    [lhelp3.extend(hh.split(',')) for hh in lhelp2 if hh!=""]
    lAx=['x','y','z']
    lelsymm=[]
    for hh in lhelp3:
        elsymmh=[0.0,0,0,0]
        strh=hh
        for inum,Ax in enumerate(lAx):
            if (strh.find(Ax)!=-1):
                if (strh.find("+"+Ax)!=-1):
                    elsymmh[inum+1]=1
                    strh="".join(strh.split("+"+Ax))
                elif (strh.find("-"+Ax)!=-1):
                    elsymmh[inum+1]=-1
                    strh="".join(strh.split("-"+Ax))
                else:
                    elsymmh[inum+1]=1
                    strh="".join(strh.split(Ax))
        if (strh==""):
            pass
        elif (strh.find("/")!=-1):
            lhelp1=strh.split("/")
            elsymmh[0]=float(lhelp1[0])/float(lhelp1[1])
        else:
            elsymmh[0]=float(strh)
        lelsymm.append(elsymmh)
    elsymm=[]
    [elsymm.extend(hh) for hh in lelsymm]
    return elsymm

def constrsite(chi,lelsymm):
    """
    give limitation on thermal parameters or susceptibility of atoms

    chi is list of elements chi [11,22,33,12,13,23]
    chic is list of constraind elements chi [11,22,33,12,13,23]
    lelsymm is list of element of symmetry for given site [t1,r11,r12,r13,t2,r21,r22,r23,t3,r31,r32,r33]
    """
    chic=[hh for hh in chi]
    [c1,c2,c3,c4,c5,c6]=chi
    mCond=[]
    for els in lelsymm:
        m1=(0,1,2,3,4,5)
        m2=m1
        kl1=((0,0),(1,1),(2,2),(0,1),(0,2),(1,2))
        kl2=kl1
        [b1,r11,r12,r13,b2,r21,r22,r23,b3,r31,r32,r33]=els
        mR=[[r11,r12,r13],[r21,r22,r23],[r31,r32,r33]]#elemets of rotation matrix is integers: -1,0,1 and nothing else.
        miR=[[int(hh2) for hh2 in hh1] for hh1 in numpy.linalg.inv(mR)]
        mA=[[0 for hh2 in m2] for hh1 in m1]
        for nm1,nkl1 in zip(m1,kl1):
            i1,j1=nkl1
            for nm2,nkl2 in zip(m2,kl2):
                k1,l1=nkl2
                mA[nm1][nm2]=mR[i1][k1]*miR[l1][j1]
                if ((nm1>2)|(nm2>2)):
                    mA[nm1][nm2]+=mR[i1][l1]*miR[k1][j1]
        for nm1 in m1:
            cond1=sum([abs(hh) for hh in mA[nm1]])
            if (cond1==0):chic[nm1]=0
            for nm2 in m2:
                if ((nm1!=nm2)&(mA[nm1][nm2]!=0)):
                    chic[nm1]=chic[nm2]*mA[nm1][nm2]
                elif (nm1==nm2)&(mA[nm1][nm2]==-1):
                    chic[nm1]=0
    return chic


def splitlinewminuses(line):
    """Split line even if between numbers there is no spaces and only minus.\n
    If it is possible transform numbers to float or integer it will be done.\n
    """
    lbefore1,lbefore2='',''
    n1=''
    splitedline=[]
    for letter in line:
        if letter.isdigit():
            n1+=letter
        elif ((letter=='.') and (lbefore1.isdigit())):
            n1+=letter
        elif ((letter=='-') and (lbefore1=='E') and (lbefore2.isdigit())):
            n1+=letter
        elif ((letter=='-') and (lbefore1!='E')):
            if ((n1!=' ') and (n1!='')):
                splitedline.append(n1)
            n1=letter
        elif ((letter==' ') and (not((n1==' ') or (n1=='')))):
            splitedline.append(n1)
            n1=''
        else:
            n1+=letter
        lbefore2=lbefore1
        lbefore1=letter
        if (n1==' '):
            n1=''
    if ((n1!=' ') and (n1!='')):
        splitedline.append(n1)
        n1=''
    fsplitedline=[]
    for elsplitedline in splitedline:
        try:
            if elsplitedline[0]=='-':
                hhh=elsplitedline[1:]
            else:
                hhh=elsplitedline
            if hhh.isdigit():
                hh=int(elsplitedline)
            else:
                hh=float(elsplitedline)
            fsplitedline.append(hh)
        except:
            fsplitedline.append(elsplitedline)
    if (len(fsplitedline)==1):
        fsplitedline=fsplitedline[0]
    return fsplitedline

def listhkl(sthovlMIN,sthovlMAX,ucp,dsymm):
    """Give list of h,k,l reflection which are in the range from sthovlMIN to sthovlMAX\n
    ucp is unit cell parameters"""
    lhkl,lmult=[],[]
    lhklres=[]
    hklmax=[int(2*ucp[ih]*sthovlMAX) for ih in range(3)]
    hmin,hmax=-1*hklmax[0],hklmax[0]
    kmin,kmax=-1*hklmax[1],hklmax[1]
    lmin,lmax=-1*hklmax[2],hklmax[2]
    #if (dsymm["centr"]):
    hmin=0
    lorig=dsymm["orig"]
    lsymm=dsymm["elsymm"]
    for h in range(hmin,hmax+1,1):
        for k in range(kmin,kmax+1,1):
            for l in range(lmin,lmax+1,1):
                flag=(abs(sum([cmath.exp(2*math.pi*1j*(orig[0]*h+orig[1]*k+orig[2]*l)) for orig in lorig]))>0.00001)
                #flag=True
                if (flag):
                    lhkls=[(h*symm[1]+k*symm[5]+l*symm[9], h*symm[2]+k*symm[6]+l*symm[10], h*symm[3]+k*symm[7]+l*symm[11]) for symm in lsymm]
                    lhkls.extend([(-hkl[0],-hkl[1],-hkl[2]) for hkl in lhkls])
                    lhkls.sort(key=lambda x:10000*x[0]+100*x[1]+x[2])
                    if (not(lhkls[-1] in lhkl)):
                        lhkl.append(lhkls[-1])
                        lmult.append(len(set(lhkls)))
    lhklsthovl=[(hkl, calcsthovl(hkl,ucp), mult) for hkl, mult in zip(lhkl, lmult)]
    lhklsthovl.sort(key=lambda x: x[1])
    lhklres = [hklsthovl[0] for hklsthovl in lhklsthovl if ((hklsthovl[1]>sthovlMIN)&(hklsthovl[1]<sthovlMAX))]
    lmultres = [hklsthovl[2] for hklsthovl in lhklsthovl if ((hklsthovl[1]>sthovlMIN)&(hklsthovl[1]<sthovlMAX))]
    return lhklres,lmultres

def calck(hkl,UB):
    if tuple(hkl) == (0, 0, 0):
        return [0., 0., 0.]
    kh=[]
    for hUB in UB:
        kh.append(sum([hh1*hh2 for hh1,hh2 in zip(hUB,hkl)]))
    normkh=(sum([hh**2 for hh in kh]))**0.5
    k=[hh*1./normkh for hh in kh]
    return k

def calcFR(FN,msfperp,Pu,Pd):
    """
    msfperp is  component of magnetic structure factor perpendicular to k_s-k_i
    pU,pD is spin direction of inicident neutrons for flipper UP and DOWN
    FN is nuclear structure factor
    Important: msfperp, pU and pD should be given in the same Cartesian coordinate system (but it does not matter in which global or local)
    """
    FNsq=abs(FN)**2
    FMsqperp=sum([abs(hh1)**2 for hh1 in msfperp])
    PuFM=sum([hh1*hh2 for hh1,hh2 in zip(Pu,msfperp)])
    PdFM=sum([hh1*hh2 for hh1,hh2 in zip(Pd,msfperp)])
    Iu=FNsq+2.*(FN.real*PuFM.real+FN.imag*PuFM.imag)+FMsqperp
    Id=FNsq+2.*(FN.real*PdFM.real+FN.imag*PdFM.imag)+FMsqperp
    FR=Iu*1./Id
    return FR

def calcCorrM(invHessian):
    """
    calculate the correlation matrix.
    """
    lcorrM=[]
    lrdiag=[abs(hh1[ihh1])**0.5 for ihh1,hh1 in enumerate(invHessian)]
    lcorrM=[[hh3*1./(hh4*hh2) for hh3,hh4 in zip(hh1,lrdiag)] for hh1,hh2 in zip(invHessian,lrdiag)]
    return lcorrM


#EXTINCTION PART, BEGINNING
def calc_fr_extinc(hkl,fn,fm_perp,e_up,ext,p_up,p_down,ucp,wavelength):
    #qhkl = calcqhkl(hkl,ext)
    fn_sq = abs(fn)**2
    fm_perp_sq = sum([abs(hh)**2 for hh in fm_perp])
    fm_perp_e_up = sum([hh1*hh2 for hh1,hh2 in zip(fm_perp,e_up)])
    fm_perp_e_up_sq = abs(fm_perp_e_up)**2
    
    #fm_perp_z = [fm_perp_e_up*hh1  for hh1 in e_up]
    #fm_perp_xy = [hh1-hh2  for hh1,hh2 in zip(fm_perp,fm_perp_z)]
    #fm_perp_xy_sq = sum([abs(hh)**2 for hh in fm_perp_xy])
    fnp = (fm_perp_e_up*fn.conjugate()+fm_perp_e_up.conjugate()*fn).real
    fp_sq = fn_sq+fm_perp_sq+fnp
    fm_sq = fn_sq+fm_perp_sq-fnp
    #fpm_sq = fm_perp_xy_sq
    fpm_sq = fm_perp_sq - fm_perp_e_up_sq

    """
    yp = calcyext(hkl, qhkl, fp_sq, wavelength, ucp)
    ym = calcyext(hkl, qhkl, fm_sq, wavelength, ucp)
    ypm = calcyext(hkl, qhkl, fpm_sq, wavelength, ucp)
    
    """
    yp = calcyext_bl_spher(hkl, fp_sq, wavelength, ucp, ext[0], ext[1], "gauss")
    ym = calcyext_bl_spher(hkl, fm_sq, wavelength, ucp, ext[0], ext[1], "gauss")
    ypm = calcyext_bl_spher(hkl, fpm_sq, wavelength, ucp, ext[0], ext[1], "gauss")
    

    pppl = 0.5*((1+p_up)*yp+(1-p_up)*ym)
    ppmin= 0.5*((1-p_down)*yp+(1+p_down)*ym)
    pmpl = 0.5*((1+p_up)*yp-(1-p_up)*ym)
    pmmin= 0.5*((1-p_down)*yp-(1+p_down)*ym)

    I_p = (fn_sq+fm_perp_e_up_sq)*pppl + pmpl*fnp + ypm*fpm_sq
    I_m = (fn_sq+fm_perp_e_up_sq)*ppmin + pmmin*fnp + ypm*fpm_sq

    fratio = I_p*1./I_m
    #if (hkl==(1,1,1)):
    #    print I_p,I_m,fratio
    #    print pppl,ppmin,pmpl,pmmin
    #    print yp, ym, ypm 
    #    print fpm_sq
    return fratio

def calc_int_extinc_powder(hkl,fn,fm_perp_e_up,fm_perp_sq,ext,p_up,p_down,ucp,wavelength):
    """
    supposed that q2 = 1
    """
    
    #qhkl = calcqhkl(hkl,ext)
    fn_sq = abs(fn)**2

    fnp = (fm_perp_e_up*fn.conjugate()+fm_perp_e_up.conjugate()*fn).real
    fp_sq = fn_sq+fm_perp_sq+fnp
    fm_sq = fn_sq+fm_perp_sq-fnp
    fpm_sq = 0.0
    """
    yp = calcyext(hkl, qhkl, fp_sq, wavelength, ucp)
    ym = calcyext(hkl, qhkl, fm_sq, wavelength, ucp)
    ypm = calcyext(hkl, qhkl, fpm_sq, wavelength, ucp)
    """
    yp = calcyext_bl_spher(hkl, fp_sq, wavelength, ucp, ext[0], ext[1], "gauss")
    ym = calcyext_bl_spher(hkl, fm_sq, wavelength, ucp, ext[0], ext[1], "gauss")
    ypm = calcyext_bl_spher(hkl, fpm_sq, wavelength, ucp, ext[0], ext[1], "gauss")    

    pppl = 0.5*((1+p_up)*yp+(1-p_up)*ym)
    ppmin= 0.5*((1-p_down)*yp+(1+p_down)*ym)
    pmpl = 0.5*((1+p_up)*yp-(1-p_up)*ym)
    pmmin= 0.5*((1-p_down)*yp-(1+p_down)*ym)

    I_p = pppl*(fn_sq+fm_perp_sq) + pmpl*fnp + ypm*fpm_sq
    I_m = ppmin*(fn_sq+fm_perp_sq) + pmmin*fnp + ypm*fpm_sq
    #print "{0[0]:2}{0[1]:2}{0[2]:2} {1:6.2f} {2:6.2f} {3:6.2f} {4:6.2f} {5:6.2f}".format(hkl, fn_sq, fm_perp_sq, pmpl, fnp, fpm_sq)
    #print "{:8.2f}{:8.2f}".format(I_p,I_m)
    return I_p, I_m




def calcyext_bl_spher(hkl, f_sq, wavelength, ucp, r, g, model):
    #f_sq in 10-12cm
    #extinction for spherical model
    K = 1.
    vol = calcVolOfUC(ucp)
    sthovl = calcsthovl(hkl, ucp)
    stheta = sthovl * wavelength
    s2theta = 2. * stheta * (1. - stheta**2)**0.5
    c2theta = 1. - 2. * stheta**2

    q = (f_sq*K/vol**2)*(wavelength**3)*1./s2theta

    t = 1.5*r
    alpha = 1.5*r*s2theta*1./wavelength
    x = 2./3*q*alpha*t

    A = 0.20 + 0.45 * c2theta
    B = 0.22 - 0.12 * (0.5-c2theta)**2
    yp = (1+2*x+(A*x**2)*1./(1.+B*x))**(-0.5)
    
    if alpha == 0.:
        ag, al = 0., 0.
    else:
        ag = alpha*g*(g**2+0.5*alpha**2)**(-0.5)
        al = alpha*g*1./(g+alpha*2./3.)
    
    
    if model == "gauss":
        xs = 2./3.*q*ag*t
        A = 0.58 + 0.48 * c2theta + 0.24 * c2theta**2
        B = 0.02 - 0.025 * c2theta
        #print "A, B", A, B
        ys = (1+2.12*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
    elif model == "lorentz":
        xs = 2./3.*q*al*t
        A = 0.025 + 0.285 * c2theta
        if c2theta>0:
            B = 0.15 - 0.2 * (0.75-c2theta)**2
        else:
            B = -0.45 * c2theta
        ys = (1+2*xs+(A*xs**2)*1./(1+B*xs))**(-0.5)
    else:
        ys = 1.
    #print "ys", ys
    yext = yp * ys
    return yext


def calcyext(hkl, qhkl, f_sq, wavelength, ucp):
    """Calculate extinction correction accordin Ext-Model=4 in FullProff programm.\n
    lqhkl is tensor defined by six variables to calculate lqhkl see function calcq(lq,lh,lk)
    lFsq is list of absolute value of structure factor, in (10-12cm)**2,\n
    wavelength is wavelength in Angstrem\n
    lttheta is list of 2 theta in degrees"""
    sthovl = calcsthovl(hkl,ucp)
    coeff = 0.001/(4)
    if ((sthovl*wavelength < 1)&(sthovl!=0.)):
        sttheta = math.sin(2*math.asin(sthovl*wavelength))
        yext = (1. + coeff * qhkl * f_sq * (wavelength**3) * 1./(sttheta*sthovl**2))**(-0.5)
    else:
        yext = 0.0
    return yext

def calcqhkl(hkl,cext):
    """Calculate anisotropic extinction parameter: qhkl=c1 h**2+c2 k**2+c3 l**2+c4 h*k+c5 h*l+c6 k*l, where
    cext=[c1,c2,c3,c4,c5,c6]
    """
    [c1,c2,c3,c4,c5,c6] = cext[0:6]
    [h,k,l] = hkl[0:3]
    qhkl = c1*h**2+c2*k**2+c3*l**2+c4*h*k+c5*h*l+c6*k*l
    return qhkl

#to delete from extinction part

def calcextparameters(ucp,cext,pUP,pDOWN,lhkl,c4FN,c4FM,wavelength,lq2,lFN,lFM):
    """ calculate extintiction together with intermediate parameters like in FullProf Ext-Mod=4
    """
    lqhkl,lyp,lym,lypm,lpppl,lppmin,lpmpl,lpmmin,lIpl,lImin,lRext=[],[],[],[],[],[],[],[],[],[],[]

    lFNc=[h*c4FN for h in lFN]
    lFMc=[h*c4FM for h in lFM]

    lqhkl=calcqhkl(lhkl,cext)

    lFpsq,lFmsq,lFpmsq=calcsqFpFmFpm(lFNc,lFMc,lq2)

    lyp=calcyext(lhkl,lqhkl,lFpsq,wavelength,ucp)
    lym=calcyext(lhkl,lqhkl,lFmsq,wavelength,ucp)
    lypm=calcyext(lhkl,lqhkl,lFpmsq,wavelength,ucp)

    lpppl = [0.5*((1+pUP  )*yp+(1-pUP  )*ym) for yp,ym in map(None,lyp,lym)]
    lppmin= [0.5*((1-pDOWN)*yp+(1+pDOWN)*ym) for yp,ym in map(None,lyp,lym)]
    lpmpl = [0.5*((1+pUP  )*yp-(1-pUP  )*ym) for yp,ym in map(None,lyp,lym)]
    lpmmin= [0.5*((1-pDOWN)*yp-(1+pDOWN)*ym) for yp,ym in map(None,lyp,lym)]

    lIpl,lImin,lRext=[],[],[]
    for FN,FM,q2,pppl,pmpl,ppmin,pmmin,ypm  in map(None,lFNc,lFMc,lq2,lpppl,lpmpl,lppmin,lpmmin,lypm):
        FNsq=(FN*FN.conjugate()).real
        FMsq=(FM*FM.conjugate()).real
        FNM=(FN*FM.conjugate()+FM*FN.conjugate()).real
        Ipl  = (FNsq+FMsq*q2**2)*pppl + FNM*pmpl*q2 + FMsq*ypm*q2*(1-q2)
        Imin = (FNsq+FMsq*q2**2)*ppmin+ FNM*pmmin*q2+ FMsq*ypm*q2*(1-q2)
        lIpl.append(Ipl)
        lImin.append(Imin)
        lRext.append(Ipl/Imin)
    return lqhkl,lyp,lym,lypm,lpppl,lppmin,lpmpl,lpmmin,lIpl,lImin,lRext


def calcsqFpFmFpm(lFNc,lFMc,lq2):
    """
    calculate Fm,Fp and Fmp for extinction correction
    """
    lFpsq,lFmsq,lFpmsq=[],[],[]
    for FN,FM,q2 in map(None,lFNc,lFMc,lq2):
        Fpsq = (abs(FN)**2) + (FN*FM.conjugate()+FM*FN.conjugate()).real*q2 + (abs(FM)**2)*q2
        Fmsq = (abs(FN)**2) - (FN*FM.conjugate()+FM*FN.conjugate()).real*q2 + (abs(FM)**2)*q2
        Fpmsq= (abs(FM)**2)*q2*(1-q2)
        lFpsq.append(Fpsq)
        lFmsq.append(Fmsq)
        lFpmsq.append(Fpmsq)
    return lFpsq,lFmsq,lFpmsq

#EXTINCTION PART, ENDING
