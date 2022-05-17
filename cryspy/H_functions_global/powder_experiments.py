
import numpy
from cryspy.B_parent_classes.cl_4_global import GlobalN
from cryspy.A_functions_base.structure_factor import calc_bulk_susceptibility_by_dictionary
from cryspy.A_functions_base.unit_cell import calc_volume_uc_by_unit_cell_parameters

def report_powder_experiments(rcif_object: GlobalN):
    ls_out = []
    rcif_dict = rcif_object.get_dictionary()
    rcif_dict_keys = rcif_dict.keys()
    l_crystal_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("crystal_")]
    l_pd_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("pd_")]
    l_pd2d_dict_keys = [hh for hh in rcif_dict_keys if hh.startswith("pd2d_")]
    for pd_dict_key in l_pd_dict_keys+l_pd2d_dict_keys:
        pd_dict = rcif_dict[pd_dict_key]

        ls_out.append(f"In experiment '{pd_dict['name']:}'\n")
        phase_scale = pd_dict["phase_scale"]
        phase_name = pd_dict["phase_name"]
        m_chi_total = numpy.zeros((3,3), dtype=float)
        l_s_v_sq, l_m_s_v, l_v_uc = [], [], []
        l_m_chi = []
        for name_sh, scale in zip(phase_name, phase_scale):
            name = f"crystal_"+name_sh
            crystal_dict = rcif_dict[name]
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
            else:
                m_chi = numpy.zeros((3,3), dtype=float)
            
            v_uc, dder = calc_volume_uc_by_unit_cell_parameters(
                crystal_dict["unit_cell_parameters"], flag_unit_cell_parameters=False)
            s_v_sq = scale * numpy.square(v_uc)
            m_s_v = m_chi * scale*v_uc
            l_v_uc.append(v_uc)
            l_s_v_sq.append(s_v_sq)
            l_m_s_v.append(m_s_v)
            l_m_chi.append(m_chi)
        sum_s_v_sq = numpy.sum(l_s_v_sq)
        ls_out.append("Phase       Volume  Fraction  Aver. moment")
        ls_out.append("                    (volume)        per UC ")

        for p_name, v_uc, s_v_sq, m_chi in zip(phase_name, l_v_uc, l_s_v_sq, l_m_chi):
            volume_fraction = s_v_sq/sum_s_v_sq
            m_chi_aver = (m_chi[0, 0]+m_chi[1, 1]+m_chi[2, 2])/3
            ls_out.append(f"{p_name:10} {v_uc:7.2f}    {volume_fraction*100:5.1f}% {m_chi_aver:6.1f} mu_B/T")

        ls_out.append("\n\nPhase          Bulk susceptibility tensor per phase")
        ls_out.append("             (11 , 22, 33, 12, 13, 23 mu_B/T per UC)")
        for p_name, m_chi in zip(phase_name, l_m_chi):
            ls_out.append(f"{p_name:10} {m_chi[0, 0]:6.1f} {m_chi[1, 1]:6.1f} {m_chi[2, 2]:6.1f} {m_chi[0, 1]:6.1f} {m_chi[0, 2]:6.1f} {m_chi[1, 2]:6.1f}")

        m_bulk = numpy.zeros((3,3), dtype=float)
        for p_name, m_s_v in zip(phase_name, l_m_s_v):
            moment = m_s_v/sum_s_v_sq*l_v_uc[0]
            m_bulk += moment
        
        m_bulk_aver = (m_bulk[0, 0] + m_bulk[1, 1] + m_bulk[2, 2])/3
        ls_out.append(f"\nBulk powder susceptibility:")
        ls_out.append(f"{m_bulk_aver:.2f} mu_B/T per volume of '{phase_name[0]:}'")
        # ls_out.append(f"\nBulk tensor susceptibility:")
        # ls_out.append(f"(mu_B/T per volume of '{phase_name[0]:}')")
        # ls_out.append(f" {m_bulk[0, 0]:5.2f} {m_bulk[0, 1]:5.2f} {m_bulk[0, 2]:5.2f}")
        # ls_out.append(f" {m_bulk[1, 0]:5.2f} {m_bulk[1, 1]:5.2f} {m_bulk[1, 2]:5.2f}")
        # ls_out.append(f" {m_bulk[2, 0]:5.2f} {m_bulk[2, 1]:5.2f} {m_bulk[2, 2]:5.2f}")

    print("\n".join(ls_out))
    return "\n".join(ls_out)
