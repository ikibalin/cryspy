Usage in python scripts
===========================================

There are several ways to use **CrysPy** library.

The `Jupyter Notebook <https://jupyter.org/>`_ is an open-source 
web application 
that allows you to create and share documents that contain 
live code, equations, visualizations and narrative text.
It becames a popular instrument for the data analysis.
In the next sections the restricted number of example shows the general methods which allows you
to create and to manipulate with given objects. 


Initialization 
-------------------

The given below python script demonstrate the application of 
the library to main manipulation with the objects:
**(i)** direct creating of the object; **(ii)** import object into string;
**(iii)** access to the object attributes; **(iv)** methods of the object::

    import numpy
    import cryspy 

    # creating of the object by defined class
    cell = cryspy.Cell(length_a=8.3, 
                       length_b=8.3,
                       length_c=8.3,
                       angle_alpha=90.,
                       angle_beta=90.,
                       angle_gamma=90.)

    # import to string
    str_cell = cell.to_cif()
    print(str_cell)

    # or
    print(cell)

    # the object cell has such internal attributes as 
    # volume of unit cell, reciprocal space, B matrix
    print("volume: ", cell.volume)
    print("reciprocal_length_a: ", cell.reciprocal_length_a)
    print("reciprocal_angle_alpha: ", cell.reciprocal_angle_alpha)
    print("matrix_B: ", cell.m_B)

    # to create object from string use the class method 'from_cif'
    new_cell = cryspy.Cell.from_cif(str_cell)

    # cell object has predefined methods
    str_report = cell.report_cell()
    print(str_report)
    np_h = numpy.array([2, 2, 1], dtype=int)
    np_k = numpy.array([0, 2, 1], dtype=int)
    np_l = numpy.array([0, 0, 1], dtype=int)
    sthovl = cell.calc_sthovl(h=np_h, k=np_k, l=np_l)
    print("sin(theta)/lambda for reflections [(2,0,0), (2,2,0), (1,1,1)]", sthovl)

    # the full list of predefined methods can be easily found in CrysPy editor software
    # or in help
    help(cell)

    # all classes and functions defined in library are given in several lists:
    print(cryspy.L_ITEM_CONSTR_CLASS)
    print(cryspy.L_LOOP_CONSTR_CLASS)
    print(cryspy.L_DATA_CONSTR_CLASS)
    print(cryspy.L_GLOBAL_CONSTR_CLASS)
    print(cryspy.L_FUNCTION)

    # for each of them you can get help
    help(cryspy.L_ITEM_CONSTR_CLASS[0])

    # another example with SpaceGroup
    space_group = cell.SpaceGroup(it_number = 87)
    print(space_group.to_cif())

    # elements of symmetry and wyckoff position for given space group
    print(space_group.space_group_wyckoff)
    print(space_group.full_space_group_symop)

    print(space_group.report_space_group())

    # another example with AtomSite
    atom_fe = cryspy.AtomSite(label="Fe", symbol_type = "Fe3+", 
                       fract_x=0., fract_y=0.2, fract_z=0.2,
                       occupancy=1.0, adp_type="Uiso", u_iso_or_equiv=0.) 

    atom_o = cryspy.AtomSite(label="O", symbol_type = "O2-", 
                       fract_x=0., fract_y=0.0, fract_z=0.0,
                       occupancy=1.0, adp_type="Uiso", u_iso_or_equiv=0.) 

    atom_site_loop = AtomSiteL(item=[atom_fe, atom_o])
    print(atom_site_loop.to_cif())

Nuclear structure factor calculations 
--------------------------------------

The crystal structure are described as a container of corresponding objects::

    cryspy.Cell(data_name="magnetit", 
                cell=cell, 
                space_group=space_group, 
                atom_site=atom_site_loop)

    print(cryspy.to_cif())

    # calculation of nuclear structure factor    
    refln_obj = cryspy.calc_refln(h=np_h, k=np_k, l=np_l)
    print(refln_obj.to_cif())

    # or    
    np_fn = cryspy.calc_f_nucl(np_h, np_k, np_l)
    print("numpy array of nuclear structure factor", np_fn)



Susceptibility approach 
------------------------

in progress

Polarized neutron diffraction on single crystal 
------------------------------------------------

in progress

1D polarized neutron powder diffraction 
----------------------------------------

in progress

2D polarized neutron powder diffraction 
----------------------------------------

in progress
