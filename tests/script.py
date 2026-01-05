import cryspy

f_name = "/Users/iuriikibalin/Desktop/2k6tMnf2.rcif"
rcif_obj = cryspy.load_file(f_name)

cryspy.mempy_magnetization_density_reconstruction(rcif_obj)