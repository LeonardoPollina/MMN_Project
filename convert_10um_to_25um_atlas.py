import nrrd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# atlas_path = "/gpfs/bbp.cscs.ch/data/project/proj62/csaba/atlas/bbp_prod_files/2022/annotation_10.nrrd"
atlas_path = "data/atlas/Allen_CCFv3_2022/annotation_10.nrrd"
atlas, header = nrrd.read(atlas_path)

# Define the scaling factor
scaling_factor = 2.5

# Compute the new shape of the array
new_shape = tuple(int(dim/scaling_factor) for dim in atlas.shape)

# Generate the evenly spaced indices along each dimension of the array
indices = [np.linspace(0, dim-1, new_dim) for dim, new_dim in zip(atlas.shape, new_shape)]

# Use the indices to select a subarray of the original array
new_atlas = atlas[np.ix_(indices[0].astype(int), indices[1].astype(int), indices[2].astype(int))]

# Print the shapes of the original and new arrays
print("Original shape:", atlas.shape)
print("New shape:", new_atlas.shape)

# load 25um atlas
atlas_path25 = "data/atlas/25um/annotation_25_ccf2017.nrrd"
atlas25, header25 = nrrd.read(atlas_path25)

# save converted atlas with 25um header
converted_atlas_path = "data/atlas/25um/annotation_10_to_25.nrrd"
nrrd.write(converted_atlas_path, new_atlas, header25)

# check the save is correct
atlas25_loaded, header25_loaded = nrrd.read(converted_atlas_path)
np.testing.assert_array_equal(atlas25_loaded,new_atlas)