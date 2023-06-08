# Project

## Scientific Question :

1. Long projecting axons have different region targets. So far, algebraic topology has been utilized in the dendritic sturctures which pertain in a single brain region. Long-range axons, however, span across multiple regions with different innervation patterns. If we filter out the axons once they are out of WM and calculate topological metrics (persistent images, silhoulettes, landscapes etc) with radial, and path distance metrics, can we see a quantitative difference between target innervation patterns as well or is the pattern conserved in persistent spaces.
    1. To do this, for each neuron in janelia database, we will:
        1. Extract axons of the morphology
        2. Mask axons based on different brain regions. 
            1. Issue 1: Some might have multiple innervation patterns per region so we need to log that and use either atlas based coordinate system to see if the projection patterns differ in AP , DV , ML axes or not.
                1. To see if they are different, we can get persistent image diff and do some kind of statistical test to say its diff or not.
                    1. If they are diff, and we have enough examples of morphs, can we say that along the coordinate system, can we say they conserve the innervation patterns for different cells of the same source region or mtype. Although we dont have mtype information
                    2. If they are the same, we can average them and give it to VGG network to extract feature vector of (around 155x1) 
2. Use of persistent diagrams with multiple brain regions has not been utilized so far in the field. Can we obtain a metric to explain the topography within a brain region using algebraic topology and atlas based features ?
3.    

---

-

1. So first need we need to get is the morphologies in swc format. Then we can get a sample morphology and using the allen atlas (i assume 25um but check), we can get what each neurite correspond to which region.
2. Once we get the region_list from allen per morph axon. we should first investigate how many regions are there for each morph.  
3. Be careful in filtering per leaf region since DGmo DGgr have different allen ids. Use general brain region branches like DG, CA1, Thalamus etc . I did this filtering before for geometric analysis of brain regions
4. Get one morph, and filter its subtrees 
5. Calculate persistent images for the subtrees . Some morphs will have different number of efferent region occupations which will generate variable number of persistent images. Lets say AA0998 has 2 (ipsilateral DG and contralateral DG) but another morph has 5 with some other areas in the play.

- [x]  Check if morph downloaded are v2.5 or v3 . They are 90 perpendicular rotated and could mess up things in filtering
- [x]  Check if they use 25um atlas not 10um
- [x]  Download the swc data from janelia
- [x]  Put folders for morphologies (swc30) and allen 25um v3 atlas in the working directory
- [x]  for each position in swc file, convert position_to_allen_ids
- [x]  plot frequency of number of sections per brain region
- [ ]  (optional) plot dendogram with colorcode of different region to see stems have noisy branches or not
    - [ ]  Too hard to implement in neurom since the code structure is complicated.
- [ ]  (ideally) write a function morph.get_region_mask(acronym,section_type,with_descendants=True). But it needs to be aware of allen atlas type the morph is registered (v3, 25um) .
- [ ]  extract subtrees of morphs based on the region they are in while considering noisy branches still in subtree.
    - [ ]  adrien meeting . use luigi pipeline on sample morphs if necessary since clustering takes time.
- [ ]  Create a processed/ folder for all extracted region by region axons
- [x]  Get persistent images from each efferent brain region per morphology
- [ ]  Group persistent images from same morph in hierarchical format so we also can see if multiple imgs exist per region
    
    e.g. morph : { region1 : [persistent_img1 , persistent_img2], region2: [persistent_img3] … }
    
- [ ]  Compare persistents within regions for mult imgs
- [ ]  download other datasets with long range axon projection:
    - [ ]  SEU data

![download.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/97baeaf5-aeda-4a13-ae2d-77a563d7e64a/download.png)

![download.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/f488cc46-4d3a-47ae-98df-f89f17d4d08b/download.png)

![download.png](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/42e77ab8-7be0-42f1-b406-05a954cb63d5/download.png)

## Outliers

Sometimes axons go into regions due to registration issues to normalized atlas. We need to remove outliers from these cases

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/6676621b-149f-4fb0-9000-4c5c2047db72/Untitled.png)

Above, we see axon from EC going to DG and C1 but along the WM, it touches the postsubiculum. This is an artifact and the axon actuall does not project there. The wm bundle it uses just happens to be next to it 

![Untitled](https://s3-us-west-2.amazonaws.com/secure.notion-static.com/8f765187-60c0-4430-bade-9ec81a4afed1/Untitled.png)

These cells are labeled as granule projection principal cells

### tSNE