import numpy as np
from voxcell import VoxelData
from voxcell.nexus.voxelbrain import Atlas
from utils import *
import morphio
import neurom as nm
from voxcell.math_utils import voxel_intersection

atlas_dir = 'data/atlas/25um/'
hierarchy_path = f'{atlas_dir}/hierarchy.json'
sample_morph = 'data/morphologies/swc30/AA0998.swc'

atlas = Atlas.open(atlas_dir)
region_map = atlas.load_region_map()
mask = atlas.get_region_mask('CA1')
vd = VoxelData.load_nrrd(f'{atlas_dir}/brain_regions.nrrd')

morph = morphio.Morphology(sample_morph)
m = nm.load_morphology(sample_morph)
m_axon = m.neurites[0] # for this morph only

mask = atlas.get_region_mask('CA1')
vi = voxel_intersection(m_axon.points[:2,:3],mask)
# vi
# array([[272,  96, 272],
#        [273,  96, 272],
#        [273,  96, 271]])

# do this for every pair in first dimension: 0:2, 2:4...
indices, sub_segments = voxel_intersection(
        m_axon.points[0:2,:3], mask, return_sub_segments=True
    )
regions = vd.raw[tuple(indices.T.tolist())]

# can check the name
#region_names = [region_map.get(i,'name') for i in regions] 

# kept_path = [632,533] for DG and region X for ex
# https://bbpgitlab.epfl.ch/neuromath/user/aberchet/axon-synthesis/-/blob/main/axon_synthesis/PCSF/clustering/utils.py#L253

# adrien cut axons tuft from some data
# and use syntheiss to generate new ones

# using wm recipe. get where to connect

# luigi task to load morph folder and gets a pd dataframe

#class clusterTerminal :
# based on proximity or brain region (our project)
# choose from clustering_funcs -> brain regions

# python -m luigi --module axon_synthesis.workflows CreateInputs --local-scheduler
'''
run this in luigi.cfg file containing folder

-m luigi --module axon_synthesis.workflows CreateInputs --local-scheduler
'''