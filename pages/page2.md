---
layout: page
title: Content
permalink: /content/
---

Please use the following links to get detailed information about different objects and procedures of the CrysPy library:

# Objects
1. Item/Loop classes
    - Classes describing the crystal structures
    - Classes describing experiments
2. Data containers to describe
    - crystal structure: class `Crystal`
    - magnetic crystal structure: class `MagCrystal`
    - single diffraction experiment with polarized neutrons: class `Diffrn`
    - unpolarized and polarized neutron diffraction experiment (in equatorial plane): class `Pd`
    - unpolarized and polarized neutron diffraction experiment (2D patterns): class `Pd2d`
    - unpolarized and polarized neutron diffraction T.O.F. experiment (in equatorial plane): class `TOF`
3. Global containers
    - Data refinement of described experiments: class `RhoChi`
    - Reconstraction of magnetization density by Maximum Entropy Method: class `MEM`


# Procedures and functions
1. Base functions
    - Algebra
    - Crystallography base
    - Extinction
    - Flip ratio
    - Matrices
2. Functions operating with predefined items and loops
    - Procedure for flip ratios
3. Functions operating with data containers
    - Procedures for Maximum Entropy calculations
4. Functions operating with global containers
    - str_to_globaln, file_to_globaln, str_to_items

# How to create own object
1. Basic rules
2. The extension of the library by own objects