"""
This is a realization of symmetry description given in 
International Tables for Crystallography (2006) Vol.G, 
Chapter 4.7, pp459-466 "Symmetry dictionary (symCIF) 
by I.D. Brown

It includes classes:
      
 :SpaceGroup:             SPACE_GROUP category
 :SpaceGroupL:            looped SPACE_GROUP category
 :SpaceGroupSymop:       SPACE_GROUP_SYMOP category
 :SpaceGroupSymopL:      looped SPACE_GROUP_SYMOP category
 :SpaceGroupWyckoff:       SPACE_GROUP_WYCKOFF category
 :SpaceGroupWyckoffL:      looped SPACE_GROUP_WYCKOFF category

and additional files:

 :it_tables.txt:             a file contains symmetry elements for each space group with standard cell definition 
 :CONSTANT_AND_FUNCTIONS.py: a set of constants and functions are save in the python code. It can be used as separate console handbook.

Use help function to get more detailed description about used classes::

    import cryspy 
    help(cryspy.SpaceGroup)

"""
__version__ = '1.0.0'
__author__ = 'Iurii Kibalin'
__contact__ = 'iurii.kibalin@cea.fr'
