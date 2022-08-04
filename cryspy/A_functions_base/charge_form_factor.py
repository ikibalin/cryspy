import os
import numpy
import json
import pycifstar
from cryspy.A_functions_base.orbital_functions import calc_jl

f_dir = os.path.dirname(__file__)

with open(os.path.join(f_dir, "atom_weight.txt")) as fid:
    data = fid.read()

D_ATOM_WEIGHT = json.loads(data)

with open(os.path.join(f_dir, "atom_electron_configuration.txt")) as fid:
    data = fid.read()

D_ATOM_ELECTRON_CONFIGURATION = json.loads(data)

f_atom_rho_orbital_radial_slater = os.path.join(f_dir, "atom_rho_orbital_radial_slater.rcif")
DATA_ATOM_RHO_ORBITAL_RADIAL_SLATER = pycifstar.to_data(f_atom_rho_orbital_radial_slater)

D_ATOM_RHO_ORBITAL_RADIAL_SLATER = {}
for loop in DATA_ATOM_RHO_ORBITAL_RADIAL_SLATER.loops:
    s_atom, s_shell = loop.name.strip().split("_")
    name_sh = [name.split("_")[-1] for name in loop.names]
    d_shell = {}
    D_ATOM_RHO_ORBITAL_RADIAL_SLATER[(s_atom.capitalize(), s_shell)] = d_shell
    for name in loop.names:
        name_sh = name.split("_")[-1]
        if name_sh == "n0":
            data_type = int
        else:
            data_type = float
        d_shell[name_sh] = numpy.array(loop[name], dtype=data_type)



f_scattering_amplitude = os.path.join(f_dir, "table_6.1.1.4.txt")
D_SCATTERING_AMPLITUDE = {}
for line in numpy.loadtxt(f_scattering_amplitude, skiprows=1, dtype=str):
    D_SCATTERING_AMPLITUDE[str(line[0])] = line[2:11].astype(float)

f_dispersion = os.path.join(f_dir, "table_4.2.6.8.txt")
D_DISPERSION = {}
with open(f_dispersion, "r") as fid:
    D_DISPERSION["table_wavelength"] = (numpy.array(fid.readline().strip().split()[2:], dtype=float))[::-1]
for line in numpy.loadtxt(f_dispersion, skiprows=1, dtype=str):
    D_DISPERSION[str(line[0])] = (line[2:12].astype(float) + 1j*line[13:23].astype(float))[::-1]


def get_n_zeta_coeff_for_atom_with_orbital(atom_name: str, orbital_name: str):
    shell = orbital_name[1:].lower()
    d_atom = D_ATOM_RHO_ORBITAL_RADIAL_SLATER[(atom_name.capitalize(), shell)]
    n = d_atom["n0"]
    zeta = d_atom["zeta0"]
    coeff = d_atom[orbital_name.lower()]
    return n, zeta, coeff


def get_atom_name_ion_charge_shell(ion_name: str):
    ion_name_sh = ion_name.strip().lower()
    flag_isotope, flag_charge = False, False
    isotope_number = ""
    charge_position = None

    for i_letter, letter in enumerate(ion_name_sh):
        if letter.isdigit():
            if i_letter == len(isotope_number):
                isotope_number += letter
                flag_isotope = True
            else:
                flag_charge = True
                charge_position = i_letter
                break

    if flag_isotope:
        isotope_position = len(isotope_number)
        isotope_number = int(isotope_number)
    else:
        isotope_position = 0
        isotope_number = None

    if flag_charge:
        atom_name = ion_name_sh[isotope_position:(charge_position+isotope_position)].capitalize()
        ion_charge = ion_name_sh[(charge_position+isotope_position):]
        ion_charge = int(f"{ion_charge[-1]:}{ion_charge[:-1]:}")
    else:
        atom_name = ion_name_sh[isotope_position:].capitalize()
        ion_charge = 0

    full_shell = get_full_shell(atom_name, ion_charge)
    return atom_name, ion_charge, full_shell

def add_ion_charge_to_shell(shell: str, ion_charge: int):
    shells = ["s", "p", "d", "f"]
    max_electrons = [2, 6, 10, 14]
    if ion_charge == 0:
        shell_out = shell
    else:
        l_orbital = shell.strip().split()
        if ion_charge > 0:
            l_orbital_new = []
            for orbital in reversed(l_orbital):
                n_shell, orbit_name, n_el = get_shell_orbit_electrons(orbital)
                if n_el > ion_charge:
                    orbital_new = f"{n_shell:}{orbit_name:}{n_el-ion_charge:}"
                    ion_charge = 0
                    l_orbital_new.append(orbital_new)
                elif n_el <= ion_charge:
                    ion_charge -= n_el
                elif ion_charge == 0:
                    l_orbital_new.append(orbital)
            l_orbital_new.reverse()
        else:
            l_orbital_new = []
            for orbital in l_orbital:
                n_shell, orbit_name, n_el = get_shell_orbit_electrons(orbital)
                n_max = max_electrons[shells.index(orbit_name)]
                n_add_max = n_max-n_el
                if n_add_max > 0:
                    if abs(ion_charge) <= n_add_max:
                        orbital_new = f"{n_shell:}{orbit_name:}{n_el-ion_charge:}"
                        ion_charge = 0
                    else:
                        orbital_new = f"{n_shell:}{orbit_name:}{n_max:}"
                        ion_charge += n_add_max
                    l_orbital_new.append(orbital_new)
                else:
                    l_orbital_new.append(orbital)
        shell_out = " ".join(l_orbital_new)
    if ion_charge != 0:
        raise UserWarning("Too much electrons for the valence shell")
    return shell_out


def get_shell_orbit_electrons(orbital:str):
    n_shell = int(orbital[0])
    orbit_name = orbital[1]
    n_el = int(orbital[2:])
    return n_shell, orbit_name, n_el


def get_full_shell(atom_name: str, ion_charge: int = 0):
    core_sh = D_ATOM_ELECTRON_CONFIGURATION[atom_name]
    l_shell = core_sh.split("]")
    if len(l_shell) == 1:
        shell_full = l_shell[0].strip()
        shell_full = add_ion_charge_to_shell(shell_full, ion_charge)
    else:
        atom_name_core = l_shell[0][1:]
        shell_start_full = get_full_shell(atom_name_core)
        shell_end = l_shell[1].strip()
        shell_end = add_ion_charge_to_shell(shell_end, ion_charge)
        shell_full = shell_start_full + " " + shell_end
    return shell_full

L_PREDEFINED_NAME = [loop.name for loop in DATA_ATOM_RHO_ORBITAL_RADIAL_SLATER.loops]

def calc_jl_for_shell_of_ion(sthovl, ion_name: str, shell: str, kappa:float = 1., l_max:int=5):
    n_shell, orbit_name = shell[0], shell[1]
    atom_name, ion_charge, full_shell = get_atom_name_ion_charge_shell(ion_name)
    n_total_el = sum([get_shell_orbit_electrons(orbital)[2] for orbital in full_shell.split()])
    jl_res = None
    for orbital in shell.strip().split():
        n_shell, orbit_name, n_el = get_shell_orbit_electrons(orbital)
        s_name = f"{atom_name.capitalize()}_{orbit_name.lower()}"
        shell = f"{n_shell:}{orbit_name:}"
        try:
            d_atom = D_ATOM_RHO_ORBITAL_RADIAL_SLATER[atom_name.capitalize(), orbit_name.lower()]
            n = d_atom["n0"]
            zeta = d_atom["zeta0"]
            coeff = d_atom[shell]
            if numpy.any(numpy.isnan(coeff)):
                raise KeyError
        except KeyError: #FIXME
            n = numpy.array([n_shell,], dtype=int)
            zeta = numpy.array([n_total_el,], dtype=float)
            coeff = numpy.array([1,], dtype=float)
            raise UserWarning(f"Parameters for '{s_name:}' are not found")

        jl = calc_jl(sthovl, coeff, n, zeta, kappa=kappa, l_max=l_max)
        if jl_res is None:
            jl_res = jl*n_el
        else:
            jl_res += jl*n_el
    return jl_res


def calc_jl_for_ion(sthovl, ion_name: str, kappa:float = 1., l_max:int=5):
    atom_name, ion_charge, full_shell = get_atom_name_ion_charge_shell(ion_name)
    jl_res = calc_jl_for_shell_of_ion(sthovl, ion_name, full_shell, kappa = kappa, l_max=l_max)
    return jl_res


def calc_scattering_amplitude_tabulated(ion_name: str, sthovl, flag_sthovl: bool = False):
    params = D_SCATTERING_AMPLITUDE[ion_name]
    a1, b1 = params[0], params[1]
    a2, b2 = params[2], params[3]
    a3, b3 = params[4], params[5]
    a4, b4 = params[6], params[7]
    c = params[8]
    
    sthovl_sq = numpy.square(sthovl)
    scattering_amplitude = (
        a1*numpy.exp(-b1*sthovl_sq) + 
        a2*numpy.exp(-b2*sthovl_sq) + 
        a3*numpy.exp(-b3*sthovl_sq) + 
        a4*numpy.exp(-b4*sthovl_sq) + 
        c
    )
    dder = {}
    if flag_sthovl:
        dder["sthovl"] = -2.*sthovl*(
            a1*b1*numpy.exp(-b1*sthovl_sq) + 
            a2*b2*numpy.exp(-b2*sthovl_sq) + 
            a3*b3*numpy.exp(-b3*sthovl_sq) + 
            a4*b4*numpy.exp(-b4*sthovl_sq)
        )
    return scattering_amplitude, dder
