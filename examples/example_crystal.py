# Work with sublibrary crystal of library NeuPy

from neupy import (Cell, 
                                SpaceGroup,
                                AtomType,
                                AtomSite, 
                                Extinction, 
                                Crystal)

#define unit cell
cell = Cell(a=7.2,alpha=62,beta=84,gamma=76,b=8)
cell.set_val(c=6.3)
cell.list_vals()
#inversed volume for unit cell
cell.get_val("ivol")
#calculation of sin(theta)/lambda for reflection (1,0,2)
cell.calc_sthovl(1,0,2)



#define unit cell: Fd-3m second origin
spgr = SpaceGroup("P1", "1")
spgr.list_vals()
#equivalent reflections for (2,1,0) in Fd-3m space group
h_eq, k_eq, l_eq, mult = spgr.calc_hkl_equiv(2, 1, 0)
ls_out = ["{:4}{:4}{:4}".format(h_1, h_2, h_3) for h_1, h_2, h_3 in zip(h_eq, k_eq, l_eq)]
print("\n".join(ls_out))


#define atoms
atom_type_1 = AtomType(type_n = "Fe", type_m="Fe3", x=0.2, y=0.1, z=0.4, flag_m=True)
atom_type_1.set_val(chi_type="cani", chi_11=2.0, chi_22=1.0, chi_33=1.5, chi_12=7.5, chi_13=0.5, chi_23=3.5)

#define atoms
atom_type_2 = AtomType(type_n = "Fe", type_m="Fe3", x=0.05, y=0.72, z=0.43, flag_m=True)
atom_type_2.set_val(chi_type="cani", chi_11=2.0, chi_22=1.0, chi_33=1.5, chi_12=7.5, chi_13=0.5, chi_23=3.5)


#define atom site in crystal
atom_site = AtomSite()
atom_site.add_atom(atom_type_1)
atom_site.add_atom(atom_type_2)


#define extinction
extinction = Extinction(domain_radius=100, mosaicity=20, model="gauss")


#define crystal
crystal = Crystal(cell=cell, space_group=spgr, atom_site=atom_site, extinction=Extinction)
# calculation of nuclear structure factor and structure factor tensor defined in Cartesian 
# crystallographic system for reflection (1,2,3)
# the values are given in 10**-12 cm**-1.
f_nucl, sft_11, sft_12, sft_13, sft_21, sft_22, sft_23, sft_31, sft_32, sft_33, d_info = crystal.calc_sf(1, 2, 3, f_print=True)


#singony of cell is changed to the "Cubic" as space group is Fd-3m
print(cell)
