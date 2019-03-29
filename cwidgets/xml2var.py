#module which allow to read and write xml files to the dictionary in accordance to structure

import xml
import xml.etree.ElementTree

def readxml(root,dstruct,rootmain):
    dout={}
    lstruc=dstruct[root.tag]
    for struc in lstruc:
        if (struc[1]=='atr'):
            val=root.get(struc[0])
            if (val!=None):
                dout[struc[0]]=root.get(struc[0])
        elif (struc[1]=='val'):
            child=root.find(struc[0])
            if (child!=None):
                try:
                    val=float(child.text)
                    fval=(child.get("refined")=="true")
                    cval=''
                except:
                    sconstr=child.text
                    sfunc=sconstr.split("[")[0]
                    lval=[]
                    for hh1 in sconstr.split("[")[1:]:
                        sway=hh1.split("]")[0].strip()
                        val=float(getxmlparam(rootmain,sway,"text"))
                        lval.append(val)
                    for inumbx,val in enumerate(lval):
                        exec("x{} = {}".format(inumbx+1,val))
                    val=eval(sfunc)
                    fval=False#fval should be false for refined values
                    cval=str(child.text)
                dout[struc[0]]=[val,fval,cval]
        elif (struc[1]=='list'):
            ldout=[]
            for child in root:
                if (struc[0]==child.tag):
                    dout2=readxml(child,dstruct,rootmain)
                    ldout.append(dout2)
            if (ldout!=[]):
                dout[struc[0]]=ldout
    return dout


def getxmlparam(xmlobject,sway,formatansw):
    res=None
    if (sway==""):
        if (formatansw=="text"):
            res=xmlobject.text
        elif (formatansw=="attrib"):
            res=xmlobject.attrib  
        elif (formatansw=="tag"):
            lres=[]
            try:
                for section in xmlobject:
                    lres.append(section.tag)
                res=tuple(lres)
            except:
                res=None
        else:
            print "MISTAKE."
        return res
    numb1,numb2=sway.find('('),sway.find(')')
    sparam1=sway[(numb1+1):numb2]
    sparam2=sway[(numb2+1):].strip()
    params=sparam1.split(",")
    nparams=len(params)
    lres=[]
    for section in xmlobject:
        dattrib=section.attrib
        if (section.tag==params[0]):
            if (nparams==2):
                if (dattrib["name"]==params[1]):
                    res=getxmlparam(section,sparam2,formatansw)
                    lres.append(res)
            elif (nparams==1):
                res=getxmlparam(section,sparam2,formatansw)
                lres.append(res)
            else:
                print "MISTAKE in getxmlparam"
                res=None   
    if (len(lres)==1):      
        lres=lres[0]
    return lres


def func4constr(root,lway):
    lway1=lway
    lhelp=getxmlparam(root,lway1,"text")
    try:
        val=float(lhelp)
    except:
        valsign=float(lhelp[:3])
        lway1=eval(lhelp[3:])
        val=valsign*float(getxmlparam(root,lway1,"text"))
    fval=(getxmlparam(root,lway1,"attrib")["refined"]=="true")
    return val, fval


