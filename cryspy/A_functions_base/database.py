import os
import pickle


f_dir = os.path.dirname(__file__)

with open(os.path.join(f_dir, "database.pickle"), "rb") as fid:
    DATABASE = pickle.load(fid)