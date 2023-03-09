import numpy
from cryspy.A_functions_base.database import DATABASE
from cryspy.A_functions_base.charge_form_factor import  get_atom_name_ion_charge_shell
from cryspy.B_parent_classes.cl_4_global import GlobalN

from cryspy.A_functions_base.structure_factor import calc_bulk_susceptibility_by_dictionary
from cryspy.A_functions_base.unit_cell import calc_volume_uc_by_unit_cell_parameters
from cryspy.E_data_classes.cl_1_crystal import Crystal
from cryspy.E_data_classes.cl_2_pd import Pd
from cryspy.E_data_classes.cl_2_pd2d import Pd2d
from cryspy.E_data_classes.cl_2_diffrn import Diffrn
from cryspy.E_data_classes.cl_2_tof import TOF

def report_powder_experiments(rcif_object: GlobalN):
    d_elements = DATABASE["Elements"]
    ls_out = []
    rcif_dict = rcif_object.get_dictionary()
    rcif_dict_keys = rcif_dict.keys()
    l_crystal_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("crystal_")]
    l_pd_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("pd_")]
    l_pd2d_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("pd2d_")]
    
    l_crystal = [item for item in rcif_object.items if isinstance(item, Crystal)]
    l_crystal_name = [item.data_name.lower() for item in l_crystal]
    
    for item in rcif_object.items:
        f_powder = isinstance(item, (Pd, Pd2d, TOF))
        f_single = isinstance(item, (Diffrn, ))
        if f_powder or f_single:
            phase = item.phase

        if f_single:
            l_phase_items = [phase, ]
        elif f_powder:
            l_phase_items = phase.items

        if f_powder or f_single:
            field = item.setup.field
            ls_out.append(f"In experiment '{item.data_name:}' \n")

            
            m_chi_total = numpy.zeros((3,3), dtype=float)
            l_s_v_sq, l_v_uc, l_m_uc = [], [], []
            l_m_chi, l_s_v_m = [], []
            l_m_chi_ion = []


            for phase_item in l_phase_items:
                phase_label = phase_item.label.lower()
                if f_single:
                    phase_scale = 1.
                else:
                    phase_scale = phase_item.scale
                ind_phase = l_crystal_name.index(phase_label)
                crystal = l_crystal[ind_phase]

                crystal.apply_constraints()
                atom_site = crystal.atom_site
                atom_multiplicity = numpy.array(atom_site.multiplicity, dtype=int)
                atom_occupancy = numpy.array(atom_site.occupancy, dtype=float)
                l_atom_type_symbol = [get_atom_name_ion_charge_shell(ion_name)[0].capitalize() for ion_name in atom_site.type_symbol]
                atomic_weight = numpy.array([ float(d_elements[("Atomic weight", atom_type)]) for atom_type in l_atom_type_symbol], dtype=float)
                
                unit_cell_weight = numpy.sum(atom_multiplicity*atom_occupancy*atomic_weight)
                coeff = 6.022*0.927*1000/unit_cell_weight

                crystal_dict = crystal.get_dictionary()
                crystal_dict_keys = crystal_dict.keys()
                if "atom_para_sc_chi" in crystal_dict_keys:
                    dict_in_out = {}
                    flag_use_precalculated_data = False
                    bulk_susceptibility, dder = calc_bulk_susceptibility_by_dictionary(
                        crystal_dict, dict_in_out, flag_use_precalculated_data = flag_use_precalculated_data)
                    m_chi = numpy.array([
                        [bulk_susceptibility[0], bulk_susceptibility[1], bulk_susceptibility[2]],
                        [bulk_susceptibility[3], bulk_susceptibility[4], bulk_susceptibility[5]],
                        [bulk_susceptibility[6], bulk_susceptibility[7], bulk_susceptibility[8]]], dtype=float)
                    atom_para_index = atom_para_index = crystal_dict["atom_para_index"]
                    n_ions = numpy.sum((atom_multiplicity*atom_occupancy)[atom_para_index])
                else:
                    m_chi = numpy.zeros((3,3), dtype=float)
                    n_ions = 1
                
                v_uc, dder = calc_volume_uc_by_unit_cell_parameters(
                    crystal_dict["unit_cell_parameters"], flag_unit_cell_parameters=False)

                s_v_sq = phase_scale * v_uc * v_uc
                s_v_m = phase_scale * v_uc * unit_cell_weight
                l_m_uc.append(unit_cell_weight)
                l_v_uc.append(v_uc)
                l_s_v_sq.append(s_v_sq)
                l_s_v_m.append(s_v_m)
                
                l_m_chi.append(m_chi*coeff)
                l_m_chi_ion.append(m_chi/n_ions)
            sum_s_v_m = numpy.sum(l_s_v_m)
            sum_s_v_sq = numpy.sum(l_s_v_sq)
            ls_out.append(f"For polycrystalline sample at field {field:.2f}T:")
            ls_out.append("-----------------------------------------------------------")
            ls_out.append("Phase      Density      wF  mass M      vF       M   Moment")
            ls_out.append("            g/cm^3           emu/g        emu/cm^3 mu_B/ion")
            ls_out.append("-----------------------------------------------------------")
            
            m_chi_mixture = 0.
            m_volume_mixture = 0.
            for phase_item, v_uc, m_uc, s_v_sq, m_chi, s_v_m, m_chi_ion in zip(l_phase_items, l_v_uc, l_m_uc, l_s_v_sq, l_m_chi, l_s_v_m, l_m_chi_ion):
                volume_fraction = s_v_sq/sum_s_v_sq
                weight_fraction = s_v_m/sum_s_v_m
                m_chi_aver = field*(m_chi[0, 0]+m_chi[1, 1]+m_chi[2, 2])/3
                density = m_uc/(v_uc*0.6022) # g/cm^3
                m_volume_aver = m_chi_aver*density
                m_chi_ion_aver = field*(m_chi_ion[0, 0]+m_chi_ion[1, 1]+m_chi_ion[2, 2])/3
                
                ls_out.append(f"{phase_item.label:10} {density:7.3f}  {weight_fraction*100:5.1f}% {m_chi_aver:7.0f}  {volume_fraction*100:5.1f}% {m_volume_aver:7.0f} {m_chi_ion_aver:8.3f}")
                m_chi_mixture += weight_fraction*m_chi_aver
                m_volume_mixture += volume_fraction*m_volume_aver
            ls_out.append("-----------------------------------------------------------")
            ls_out.append("  wF is weight-fraction of given phase in a mixture;")
            ls_out.append("  mass M is mass magnetization of powder sample;")
            ls_out.append("  vF is volume-fraction of given phase in a mixture;")
            ls_out.append("   M is volume magnetization of powder sample.\n")
            ls_out.append(f"Mass magnetization of a mixture is {m_chi_mixture:.2f} [emu/g];")
            ls_out.append(f"Volume magnetization of a mixture is {m_volume_mixture:.2f} [emu/cm^3].\n")
            
            ls_out.append("----------------------------------------------------")      
            ls_out.append("Phase           Mass susceptibility tensor x 10 000")
            ls_out.append("----------------------------------------------------")      
            for phase_item, m_chi in zip(l_phase_items, l_m_chi):
                ls_out.append(f"{phase_item.label:10} {m_chi[0, 0]:6.2f} {m_chi[1, 1]:6.2f} {m_chi[2, 2]:6.2f} {m_chi[0, 1]:6.2f} {m_chi[0, 2]:6.2f} {m_chi[1, 2]:6.2f}")
            ls_out.append("----------------------------------------------------")
            ls_out.append("The tensor is defined in Cartezian coordinate system")
            ls_out.append("such as axis Z || c, axis X || a* for each phase in ")
            ls_out.append("[emu/(Oe g)].")
    return "\n".join(ls_out)
