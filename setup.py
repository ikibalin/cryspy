#!/usr/bin/env python

# Copyright (c) 2018-2023 Iurii Kibalin   
# https://github.com/ikibalin/cryspy  
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
# 
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# 
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# 
# * Neither the name "CrysPy" nor the names of its contributors may
#   be used to endorse or promote products derived from this software
#   without specific prior written permission.
# 
# This software is provided by the copyright holders and contributors "as
# is" and any express or implied warranties, including, but not limited
# to, the implied warranties of merchantability and fitness for a
# particular purpose are disclaimed. In no event shall the copyright owner
# or contributors be liable for any direct, indirect, incidental, special,
# exemplary, or consequential damages (including, but not limited to,
# procurement of substitute goods or services; loss of use, data, or
# profits; or business interruption) however caused and on any theory of
# liability, whether in contract, strict liability, or tort (including
# negligence or otherwise) arising in any way out of the use of this
# software, even if advised of the possibility of such damage.

import os
import os.path
from setuptools import setup, find_packages

f_name = os.path.join(os.path.dirname(__file__), "readme.rst")
with open(f_name, 'r') as f:
    long_description = f.read()

setup(
    name='cryspy',
    version='0.7.7',
    description='PNPD data analysis',
    long_description = long_description,
    author='Iurii Kibalin',
    author_email='yurikibalin@outlook.com',
    url = 'https://github.com/ikibalin/cryspy',
    license          = 'MIT License',
    keywords         = 'Polarized Neutron Diffraction Data Analysis',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
    ],    
    include_package_data=True,
    packages=find_packages(),  #same as name
    install_requires=[
        'numpy', 
        'scipy', 
        'pycifstar',
        'matplotlib',
    ] 
)