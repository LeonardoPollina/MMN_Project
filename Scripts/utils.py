import json

def get_children(json_file, atlas_id=None, acronym=None):
    """
    Takes a JSON hierarchy file and an atlas_id or acronym as input.
    Returns the list of child structure ids.
    """
    with open(json_file, "r") as f:
        hierarchy = json.load(f)
        
    msg = hierarchy["msg"][0] # only index 0 exists
    
    if atlas_id is not None:
        def search_by_id(node, target_id):
            if node["id"] == target_id:
                return node
            elif "children" in node:
                for child in node["children"]:
                    result = search_by_id(child, target_id)
                    if result is not None:
                        return result
            return None
        
        target_node = search_by_id(msg, atlas_id)
        
    elif acronym is not None:
        def search_by_acronym(node, target_acronym):
            if node["acronym"] == target_acronym:
                return node
            elif "children" in node:
                for child in node["children"]:
                    result = search_by_acronym(child, target_acronym)
                    if result is not None:
                        return result
            return None
        
        target_node = search_by_acronym(msg, acronym)
    
    if target_node is not None and "children" in target_node:
        return [child["id"] for child in target_node["children"]]
    else:
        return []
    

def idx_to_voxel_id(voxeldata,idx,idy,idz):
    return voxeldata.raw[idx,idy,idz]

def idxs_to_voxel_ids(voxeldata,indices):
    id_list = []
    for idx,idy,idz in indices:
        id_list.append(voxeldata.raw[idx,idy,idz])
    return id_list

def voxel_ids_to_region_names(region_map,region_ids,param='name'):
    '''
    gets attr from a list of region ids . attr can be name or acronym etc.
    '''
    attr_list = []
    
    if isinstance(region_ids,int):
        return region_map.get(region_ids,param)
    
    for region_id in region_ids:
        attr_list.append(region_map.get(region_id, param)) # i.e. 632th annotation
    return attr_list


def read_swc(filename):
    """
    Read an SWC file and return the segments as a list of dictionaries.
    """
    segments = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split()
            segment = {
                'id': int(fields[0]),
                'type': int(fields[1]),
                'x': float(fields[2]),
                'y': float(fields[3]),
                'z': float(fields[4]),
                'radius': float(fields[5]),
                'parent': int(fields[6]) - 1  # SWC parent ID is 1-based
            }
            segments.append(segment)
    return segments