import glob
import os
import pickle as pickle
import numpy as np
import hdf5storage
from scipy.io import loadmat
from os import mkdir, path
from matplotlib import image


def save_variable(variable, name):
    final_name = name + '.pickle'
    with open(final_name, 'wb') as f:
        pickle.dump(variable, f, pickle.HIGHEST_PROTOCOL)

        
def load_pickle_variable(variable_name):
    filename = variable_name +  '.pickle'
    with open(filename, 'rb') as f:
        variable = pickle.load(f)
    return variable