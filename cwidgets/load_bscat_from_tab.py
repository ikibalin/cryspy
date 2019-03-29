def reduce_label(label):
    llabel_reduced = []
    flag = False
    for hh in label.strip():
        if ((hh.isdigit()|(hh == "_"))&flag):
            break
        else:
            llabel_reduced.append(hh)
            if (not (hh.isdigit()|(hh == "_"))):
                flag = True
    return "".join(llabel_reduced)

def read_bscat_from_tab(lname):
    fname = "bscat.tab"
    fid = open(fname,'r')
    lcontent = fid.readlines()
    fid.close
    ldata = [line.split() for line in lcontent]
    lbscat = []
    for label in lname:
        bscat =0.0
        lab_reduced = reduce_label(label)
        for hh in ldata:
            if hh[0].upper() == lab_reduced.upper():
                print " ".join(hh)
                try:
                    bscat = float(hh[2])*0.1
                except:
                    bscat = 0.
                lbscat.append(bscat) 
                break
    return lbscat

def read_tab(fname):
    dtab = {}
    return dtab

def main():
    fname = "bscat.tab"
    dtab = read_tab(fname)
    

if __name__ == "__main__":
    main()