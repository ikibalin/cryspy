## Constitution of the library

The library could be separadet on eight categories. 
The functions and objects from each category can 
import objects predefiened earlier.

1. Function, used the standard modules.
2. Common object used as parent classes
3 Item/Loop objects based on the parent classes with defined methods.
4. Functions, which use the Item/Loop objects.
5. Data objects, which includes Item/Loop objects
6. Functions, which use the Data objects
7. Global objects, which includes Item/Loop objects and Data objects
8. Functions, which use the Global objects

Import of the objects is performed consecutively from 1. until 8.


## Naming

The functions are defined in files named as `functions_#level_#tag.py` where `#tag` 
is short description of defined functions, `#level` marks dependences from object defined 
in the same category (it is 1 if there are no dependencies, it is 2 if there is dependence
from functions defined with level 1, it is 3 if there is dependence from functions 
defined with level 2 and so on)

The class for each Item/Loop/Data/Global object is defined in files named as 
`cl_#level_#class_name.py`, where `#class_name.py` is the name of the class.

## Directory

A_functions_base
B_parent_classes
C_item_loop_classes
D_functions_item_loop
E_data_classes
F_functions_data
G_global_classes
H_functions_global

## Function names

The function should be named as `#verb_#output_by_#input`