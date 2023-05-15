import tmd
import glob,os
from tmd import view 
import numpy as np
import matplotlib.pyplot as plt
from tmd.Topology import analysis
from tmd.view import common as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pickle
import seaborn as sns
from tqdm import tqdm

rc={"axes.labelsize":7, "legend.fontsize":5, "axes.linewidth":0.4,"axes.titlesize":8,
            "xtick.labelsize":10, "xtick.major.size":2, "xtick.major.width":0.5, "xtick.minor.size":1.5, "xtick.minor.width":0.3,
            "ytick.labelsize":10, "ytick.major.size":2, "ytick.major.width":0.5, "ytick.minor.size":1.5, "ytick.minor.width":0.3}

plt.rcParams.update(rc)

def normalize_phs(p1):
    return np.divide(p1, np.max(p1)).tolist()

def load_population(morph_dir:str,sample=False):
    ''' Load any pop of neurons from folder '''
    L = glob.glob(morph_dir)
    np.random.shuffle(L_clones)
    
    if sample:
        pop_tmd = tmd.io.load_population(L, use_morphio=True)
    else:
        pop_tmd = tmd.io.load_population(L[:sample], use_morphio=True)

    phs_neurite = [tmd.methods.get_ph_neuron(n, neurite_type=neurite_type) for n in pop_tmd.neurons]
    phs_neurite_norm = [normalize_phs(p) for p in phs_neurite if len(p)>0]
    return phs_neurite,phs_neurite_norm
    
def compare_clone_strategy(orig_morphs,repaired_morphs,clone_morphs,neurite_type='dendrites',plot=True,save_file_name=False):
    L_rep = glob.glob(orig_morphs)
    pop_rep_tmd = tmd.io.load_population(L_rep, use_morphio=True)

    L_or = glob.glob(repaired_morphs)
    pop_or_tmd = tmd.io.load_population(L_or, use_morphio=True)

    L_clones = glob.glob(clone_morphs)
    np.random.shuffle(L_clones)

    pop_clones_tmd = tmd.io.load_population(L_clones[:500], use_morphio=True)

    phs_den_original = [tmd.methods.get_ph_neuron(n, neurite_type=neurite_type) for n in pop_or_tmd.neurons]
    phs_den_clones = [tmd.methods.get_ph_neuron(n, neurite_type=neurite_type) for n in pop_clones_tmd.neurons]
    phs_den_repaired = [tmd.methods.get_ph_neuron(n, neurite_type=neurite_type) for n in pop_rep_tmd.neurons]

    phs_den_repaired_norm = [normalize_phs(p) for p in phs_den_repaired if len(p)>0]
    phs_den_clones_norm = [normalize_phs(p) for p in phs_den_clones if len(p)>0]
    phs_den_original_norm = [normalize_phs(p) for p in phs_den_original if len(p)>0]

    #xlab = 'End radial distance'
    #ylab = 'Start radial distance'
    xlab = ''
    ylab = ''

    if plot: 
        fig = plt.figure(figsize=(15,15))
        ax1 = fig.add_subplot(331)
        view.plot.diagram(tmd.analysis.collapse(phs_den_original), color='b', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab)
        ax1.set_ylim(-10, 1000)
        ax1.set_xlim(-10, 1000)
        ax1 = fig.add_subplot(334)
        view.plot.diagram(tmd.analysis.collapse(phs_den_repaired), color='g', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab, title='')
        ax1.set_ylim(-10, 1000)
        ax1.set_xlim(-10, 1000)
        ax1 = fig.add_subplot(337)
        view.plot.diagram(tmd.analysis.collapse(phs_den_clones), color='r', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab, title='')
        ax1.set_ylim(-10, 1000)
        ax1.set_xlim(-10, 1000)
        ax1 = fig.add_subplot(332)
        view.plot.diagram(tmd.analysis.collapse(phs_den_original_norm), color='b', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab, title='Persistence diagram (normalized)')
        ax1.set_ylim(-0.1, 1.1)
        ax1.set_xlim(-0.1, 1.1)
        ax1 = fig.add_subplot(335)
        view.plot.diagram(tmd.analysis.collapse(phs_den_repaired_norm), color='g', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab, title='')
        ax1.set_ylim(-0.1, 1.1)
        ax1.set_xlim(-0.1, 1.1)
        ax1 = fig.add_subplot(338)
        view.plot.diagram(tmd.analysis.collapse(phs_den_clones_norm), color='r', alpha=0.9, s=30, new_fig=False, xlabel=xlab, ylabel=ylab, title='')
        ax1.set_ylim(-0.1, 1.1)
        ax1.set_xlim(-0.1, 1.1)
        ax1 = fig.add_subplot(333)
        Z0 = view.plot.persistence_image_average(phs_den_original_norm, new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab)
        ax1 = fig.add_subplot(336)
        Z1 = view.plot.persistence_image_average(phs_den_repaired_norm, new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title='')
        ax1 = fig.add_subplot(339)
        Z2 = view.plot.persistence_image_average(phs_den_clones_norm, new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title='')
        plt.tight_layout()

        if save_file_name:
            plt.savefig(f'{save_file_name}.png')
            #plt.savefig(f'{save_file_name}.svg')

        out_dict = {'orig':phs_den_original_norm,'repaired':phs_den_repaired_norm,'cloned':phs_den_clones_norm}
        return fig , out_dict
    else:
        out_dict = {'orig':phs_den_original_norm,'repaired':phs_den_repaired_norm,'cloned':phs_den_clones_norm}
        return None,out_dict

# make compare_clone_strategy for single neuron group
def get_persistent_pop(morph_path:str,neurite_type:str,normalize=False):
    L_rep = glob.glob(morph_path)
    pop_rep_tmd = tmd.io.load_population(L_rep, use_morphio=True)
    phs_den_original = [tmd.methods.get_ph_neuron(n, neurite_type=neurite_type) for n in pop_or_tmd.neurons]
    if normalize:
        phs_den_original_norm = [normalize_phs(p) for p in phs_den_original if len(p)>0]
        return phs_den_original, phs_den_original_norm
    else:
        return phs_den_original, _

def persistence_image_diff(Z1, Z2, new_fig=True, subplot=111, xlims=None, ylims=None,
                           norm=True, vmin=-1., vmax=1., cmap='bwr',add_colorbar=True, **kwargs):
    """Takes as input two images as exported from the gaussian kernel
       plotting function, and plots their difference. Original function from tmd

    Parameters
    ----------
    Z1 : np.array
        Persistence image 1
    Z2 : np.array
        Persistence image 2
    new_fig : bool  (default: True)
        Whether to create a new figure
    subplot : int (default: 111)
        Subplot to plot in
    xlims : tuple (default: None)   
        X-axis limits
    ylims : tuple (default: None)   
        Y-axis limits   
    norm : bool (default: True)
        Whether to normalize the images
    vmin : float (default: -1.)
        Minimum value for the colormap
    vmax : float (default: 1.)
        Maximum value for the colormap
    cmap : str (default: 'bwr')
        Colormap to use
    add_colorbar : bool (default: True)
        Whether to add a colorbar
    **kwargs : dict 
        Keyword arguments to pass to the plotting function

    Returns
    -------
    cm.plot_style
        Plot style object
    """
    if xlims is None or xlims is None:
        xlims, ylims = ((0, 100), (0, 100))

    difference = analysis.get_image_diff_data(Z1, Z2, normalized=norm)
    fig, ax = cm.get_figure(new_fig=new_fig, subplot=subplot)
    im = ax.imshow(np.rot90(difference), vmin=vmin, vmax=vmax, cmap=cmap,
              interpolation='bilinear', extent=xlims + ylims)
    if add_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

    kwargs['xlim'] = xlims
    kwargs['ylim'] = ylims
    return cm.plot_style(fig=fig, ax=ax, **kwargs)


def get_diff_data(Z0,Z1,norm_factor=None, bw_method=None, xlims=None, ylims=None,**kwargs):
    '''
    To be used after get_persistent_pop to get the persistent image data for 2 persistent images

    Parameters
    ----------
    Z0 : np.array
        Persistence image 1
    Z1 : np.array
        Persistence image 2
    norm_factor : float (default: None)
        Normalization factor for the persistence image
    bw_method : float (default: None)
        Bandwidth method for the gaussian kernel
    xlims : tuple (default: None)
        X-axis limits
    ylims : tuple (default: None)
        Y-axis limits
    **kwargs : dict
        Keyword arguments to pass to the tmd.Topology.analysis.get_persistence_image_data
    '''
    collapsed = tmd.analysis.collapse(Z0)
    #Z0, _ = view.plot.persistence_image(collapsed)
    Z0 = analysis.get_persistence_image_data(collapsed,**kwargs)
    
    collapsed = tmd.analysis.collapse(Z1)
    #Z1, _= view.plot.persistence_image(collapsed)
    Z1 = analysis.get_persistence_image_data(collapsed,**kwargs)

    return Z0,Z1

def plot_persistent_diff(Z0,Z1,savefig=False,**kwargs):
    ''' To be used after get_diff_data to plot the persistent diff plot for 2 persistent images '''
    fig = persistence_image_diff(Z0,Z1,**kwargs)
    
    if savefig:    
        plt.savefig(f'{savedir}.png')
        plt.savefig(f'{savedir}.svg')
    
    return fig

def persistence_image_average(ph_list, new_fig=True, subplot=111, xlims=None,
                              ylims=None, norm_factor=1.0, vmin=None, vmax=None,
                              cmap='jet', weighted=False,add_colorbar=False, **kwargs):
    """
    Merges a list of ph diagrams and plots their respective average image.

    Parameters
    ----------
    ph_list : list
        List of ph diagrams 
    new_fig : bool  (default: True)
        Whether to create a new figure
    subplot : int (default: 111)
        Subplot to plot in  
    xlims : tuple (default: None)
        X-axis limits
    ylims : tuple (default: None)   
        Y-axis limits
    norm_factor : float (default: 1.0)
        Normalization factor for the persistence image
    vmin : float (default: None)
        Minimum value for the colormap
    vmax : float (default: None)
        Maximum value for the colormap
    cmap : str (default: 'jet')
        Colormap to use
    weighted : bool (default: False)
        Whether to weight the average by the number of points in each diagram
    add_colorbar : bool (default: False)
        Whether to add a colorbar
    **kwargs : dict
        Keyword arguments to pass to the plotting function

    Returns
    -------
    av_imgs : np.array
        Average image
    cm.plot_style
        Plot style object
    """
    # pylint: disable=unexpected-keyword-arg
    av_imgs = analysis.get_average_persistence_image(ph_list, xlims=xlims, ylims=ylims,
                                                     norm_factor=norm_factor, weighted=weighted)
    if xlims is None or xlims is None:
        xlims, ylims = analysis.get_limits(ph_list)

    if vmin is None:
        vmin = np.min(av_imgs)
    if vmax is None:
        vmax = np.max(av_imgs)

    fig, ax = cm.get_figure(new_fig=new_fig, subplot=subplot)
    im = ax.imshow(np.rot90(av_imgs), vmin=vmin, vmax=vmax, cmap=cmap,
              interpolation='bilinear', extent=xlims + ylims)
    
    if add_colorbar:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        fig.colorbar(im, cax=cax, orientation='vertical')

    kwargs['xlim'] = xlims
    kwargs['ylim'] = ylims
    kwargs['title'] = kwargs.get('title', 'Average persistence image')
    kwargs['xlabel'] = kwargs.get('xlabel', 'End radial distance from soma')
    kwargs['ylabel'] = kwargs.get('ylabel', 'Start radial distance from soma')

    return av_imgs, cm.plot_style(fig=fig, ax=ax, **kwargs)


if __name__ == '__main__':
    #os.chdir('./')
    #print(f'Current working directory: {os.getcwd()}')

    orig_folder = 'orig_morphs'
    #repaired_folder = 'repaired_morphs'
    #clone_folder = '08_Cloned_morphs'
    morph_format = 'swc'
    
    #xlab = 'End radial distance'
    #ylab = 'Start radial distance'
    xlab = ''
    ylab = ''
    output_dir = './output/persistent_images'
    os.makedirs(output_dir,exist_ok=True)

    file_name_to_dump = f'{output_dir}/persistent_data.pkl'


    if not os.path.exists(file_name_to_dump):
        mtype_dict = {}
        for mt in tqdm(os.listdir('orig_morphs/')):
            _, cur_basal_dict =  compare_clone_strategy(f"{orig_folder}/{mt}/*.{morph_format}",
                                                    f"{repaired_folder}/{mt}/*.{morph_format}",
                                                    f"{clone_folder}/{mt}/*.{morph_format}", 
                                                    neurite_type='basal',
                                                    plot=False,
                                                    save_file_name=f'{output_dir}/persistent_basals_{mt}')
            _, cur_axon_dict = compare_clone_strategy(f"{orig_folder}/{mt}/*.{morph_format}",
                                                    f"{repaired_folder}/{mt}/*.{morph_format}",
                                                    f"{clone_folder}/{mt}/*.{morph_format}", 
                                                    neurite_type='axon',
                                                    plot=False,
                                                    save_file_name=f'{output_dir}/persistent_axons_{mt}')
        
            mtype_dict[mt] = {'basal':cur_basal_dict,'axon':cur_axon_dict}


        _ , mtype_dict['SP_PC']['apical'] = compare_clone_strategy(f"{orig_folder}/{'SP_PC'}/*.{morph_format}",
                                                            f"{repaired_folder}/{'SP_PC'}/*.{morph_format}",
                                                            f"{clone_folder}/{'SP_PC'}/*.{morph_format}",
                                                            neurite_type='apical',
                                                            plot=False,
                                                            save_file_name=f'{output_dir}/persistent_apicals')
        
        with open(file_name_to_dump,'wb') as f:
            pickle.dump(mtype_dict, f)
        print(f'Persistent data saved to {file_name_to_dump}')

    else:
        with open(file_name_to_dump,'rb') as f:
            mtype_dict = pickle.load(f)
        print(f'Persistent data loaded from {file_name_to_dump}')
     

    processed_dict = {} # dor diff data
    for mt in os.listdir('orig_morphs/'):
        mt_basal_orig , mt_basal_clone = get_diff_data(mtype_dict[mt]['basal']['orig'],mtype_dict[mt]['basal']['cloned'])
        processed_dict[mt] = [mt_basal_orig , mt_basal_clone]

    ### BASAL DENDRITES
    fig = plt.subplots(2,6,figsize=(24,8))
    for idx,mt in enumerate(sorted(os.listdir('orig_morphs/'))):
        plt.subplot(2,6,idx+1)
        Z2 = view.plot.persistence_image_average(mtype_dict[mt]['basal']['orig'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title=mt)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_basal_orig.png')
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_basal_orig.svg')

    fig = plt.subplots(2,6,figsize=(24,8))
    for idx,mt in enumerate(sorted(os.listdir('08_Cloned_morphs/'))):
        plt.subplot(2,6,idx+1)
        Z2 = view.plot.persistence_image_average(mtype_dict[mt]['basal']['cloned'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title=mt)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_basal_clone.png')
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_basal_clone.svg')

    # DIFFERENCE
    fig = plt.subplots(2,6,figsize=(24,8))
    for idx,mt in enumerate(sorted(os.listdir('orig_morphs/'))):
        plt.subplot(2,6,idx+1)
        plot_persistent_diff(processed_dict[mt][0] , processed_dict[mt][1],
                        xlims=(-0.1, 1.1), ylims=(-0.1, 1.1),new_fig=False, subplot=131, xlabel=xlab, ylabel=ylab, title=mt,cmap='bwr')
    #plt.suptitle('Original - Clone')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/orig_clone_diff_basal_for_mt.png')
    plt.savefig(f'{output_dir}/orig_clone_diff_basal_for_mt.svg')

    ### APICAL DENDRTIES
    apical_dict =  mtype_dict['SP_PC']['apical']
    fig = plt.subplots(1,3,figsize=(12,4))
    plt.subplot(1,3,1) 
    Z0 = persistence_image_average(apical_dict['orig'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, add_colorbar=True, title='Original')
    plt.subplot(1,3,2)
    Z1 = persistence_image_average(apical_dict['cloned'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab,add_colorbar=True, title='Clone')
    plt.subplot(1,3,3)
    a1,a2 = get_diff_data(apical_dict['orig'],apical_dict['cloned'])
    plot_persistent_diff(a1,a2,add_colorbar=False,title='Original - Clone',new_fig=False, xlabel=xlab, ylabel=ylab)
    #plt.suptitle('Apical Dendrites')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/average_persistent_apical_orig_clone_diff.png')
    plt.savefig(f'{output_dir}/average_persistent_apical_orig_clone_diff.svg')

    ### AXONS
    processed_dict_axon = {}
    for mt in os.listdir('orig_morphs/'):
        mt_axon_orig , mt_axon_clone = get_diff_data(mtype_dict[mt]['axon']['orig'],mtype_dict[mt]['axon']['cloned'])
        processed_dict_axon[mt] = [mt_axon_orig , mt_axon_clone]

    fig = plt.subplots(2,6,figsize=(24,8))
    for idx,mt in enumerate(sorted(os.listdir('orig_morphs/'))):
        plt.subplot(2,6,idx+1)
        plot_persistent_diff(processed_dict_axon[mt][0] , processed_dict_axon[mt][1],
                        xlims=(-0.1, 1.1), ylims=(-0.1, 1.1),new_fig=False, subplot=131, xlabel=xlab, ylabel=ylab, title=mt,cmap='bwr')
    #plt.suptitle('Original - Clone')
    plt.tight_layout()
    plt.savefig(f'{output_dir}/orig_clone_diff_axon_for_mt.png')
    plt.savefig(f'{output_dir}/orig_clone_diff_axon_for_mt.svg')


    fig = plt.subplots(2,6,figsize=(24,8))
    neurite_type= 'axon'
    for idx,mt in enumerate(sorted(os.listdir('08_Cloned_morphs/'))):
        plt.subplot(2,6,idx+1)
        Z2 = view.plot.persistence_image_average(mtype_dict[mt][neurite_type]['cloned'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title=mt)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_{neurite_type}_clone.png')
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_{neurite_type}_clone.svg')


    fig = plt.subplots(2,6,figsize=(24,8))
    neurite_type= 'axon'
    for idx,mt in enumerate(sorted(os.listdir('orig_morphs/'))):
        plt.subplot(2,6,idx+1)
        Z2 = view.plot.persistence_image_average(mtype_dict[mt][neurite_type]['orig'], new_fig=False, xlims=(-0.1, 1.1), ylims=(-0.1, 1.1), xlabel=xlab, ylabel=ylab, title=mt)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_{neurite_type}_orig.png')
    plt.savefig(f'{output_dir}/average_persistent_per_mtype_{neurite_type}_orig.svg')
    