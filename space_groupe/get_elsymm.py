def tr2elsymm(str1):
    """transform string to element of symmetry: (x,y,-z) -> 0.0 1 0 0  0.0 0 1 0  0.0 0 0 -1
    """
    str2="".join(str1.split(" "))
    lhelp1,lhelp2,lhelp3=[],[],[]
    lhelp1=[hh for hh in str2.split('(') if hh!=""]
    [lhelp2.extend(hh.split(')')) for hh in lhelp1 if hh!=""]
    [lhelp3.extend(hh.split(',')) for hh in lhelp2 if hh!=""]
    lAx=['X','Y','Z']
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

def read_el_cards(fitables):
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

def get_symm(ldcard, number_or_name_space_origin):
    dsymm = {}
    lhelp = (number_or_name_space_origin.strip()).split()
    if len(lhelp) == 1:
        n_choise = "1"
    else:
        n_choise = lhelp[1]

    try:
        spgr_n = lhelp[0]
        spgr_name = ""
    except:
        spgr_n = ""
        spgr_name = lhelp[0]
        
    flag = False
    for dcard in ldcard:
        if (((dcard["number"] == spgr_n)|(dcard["name"] == spgr_name))&(dcard["choice"][0] == n_choise)):
            flag = True
            break
    if (not flag):
        return dsymm

    lelsymm = []
    for ssymm in dcard["symmetry"]:
        lelsymm.append(tr2elsymm(ssymm))
    centr = dcard["centr"][0]=="true"
    pcentr = [float(hh) for hh in dcard["pcentr"][0].split(",")]
    fletter = dcard["name"][0]
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
    dsymm={"centr":centr,"elsymm":lelsymm,"orig":lorig,"pcentr":pcentr}
    return dsymm

import os
os.getcwd()

fitables = "itables.txt"
number_or_name_space_origin = "227 2"
ldcard = read_el_cards(fitables)
get_symm(ldcard, number_or_name_space_origin)