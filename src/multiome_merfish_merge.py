import scanpy as sc
import numpy as np
import pandas as pd
import os
import diopy
import anndata as ad
from harmonypy import run_harmony
import shutil
import hanqing_integrate
import anndata
from scipy.sparse import csr_matrix
from tkinter import simpledialog
from tkinter import Tk
from tkinter import messagebox
from tkinter.filedialog import askopenfilename

#input = string name of experiment, adata to be processed, table of gene name to ID
def modify_adata_before_anchors(experiment_name, a1, gene_to_id_table, lfile):
    a1.obs['experiment'] = [experiment_name for i in range(len(a1.obs.index.tolist()))]
    a1.obs.index = [os.path.basename(lfile).split('.')[0]+str(i) for i in a1.obs.index.tolist()]
    sc.tl.pca(a1)
    gene_name_list = []
    gene_id_ind = gene_to_id_table['gene_id'].tolist()
    for i in a1.var.index.tolist():
        try:
            indexis = gene_id_ind.index(i)
            gname = gene_to_id_table.iloc[indexis]['gene_name']
            gene_name_list.append(gname)
        except:
            None
            gene_name_list.append(None)
    a1.var['gene_name'] = gene_name_list
    return a1, gene_name_list

def create_highly_variable(a1, adata):
    mset = set(a1.var.gene_name.tolist())
    mset.remove(None)
    intersection_list = [i for i in adata.var.index.tolist() if i in mset] 
    highly_var = []
    for i in adata.var.index.tolist():
        if i in set(intersection_list):
            highly_var.append(True)
        else:
            highly_var.append(False)
            
    adata.var['highly_variable'] = highly_var
    highly_var = []
    for i in a1.var.gene_name.tolist():
        if i in set(intersection_list):
            highly_var.append(True)
        else:
            highly_var.append(False)
    a1.var['highly_variable'] = highly_var
    adata_subset = adata[:,adata.var.highly_variable == True]
    a1 = a1[:,a1.var.highly_variable == True]

    return intersection_list, adata, a1, adata_subset

def integrate_motifs(a1):
    motif_cols = [motif for motif in a1.obs.columns if ('MA' == motif[:2]) and ('_' in motif)]
    partial_ad = anndata.AnnData(X = np.hstack([a1.X, np.array(a1.obs[motif_cols])]), var= np.append(a1.var.index.tolist(), a1.obs[motif_cols].columns.tolist()), obs = a1.obs[motif_cols].index.tolist())
    partial_ad.obs = a1.obs[[vh for vh in a1.obs.columns if vh not in motif_cols]]
    partial_ad.uns = a1.uns
    partial_ad.obsm = a1.obsm

    partial_ad.obsp = a1.obsp
    partial_ad.var.index = partial_ad.var[0]
    return motif_cols, partial_ad

def prepare_integration_multiome_object(adata_subset, a1):
    ids =  []
    for i in adata_subset.var.index.tolist():
        ids.append(a1.var.index.tolist().index(i))
    adata_subset = adata_subset[:,ids]
    sc.pp.neighbors(adata_subset)
    m_lis = a1.var.index.tolist()
    a_lis = adata_subset.var.index.tolist()
    total_subscan = None
    for i in range(len(a1.var.index.tolist())):
        current_gene = m_lis[i]
        if i == 0:
            total_subscan = adata_subset[:,a_lis.index(current_gene)]
        else: 
            total_subscan = ad.concat([total_subscan, adata_subset[:,a_lis.index(current_gene)]],join='outer',axis=1)

    df = pd.read_csv(r'C:\Users\amonell\Downloads\chromvar_data.csv', index_col=0)
    total_subscan.obs = pd.concat([adata_subset.obs, df.T], axis = 1)
    return adata_subset, total_subscan

def update_integration_multiome_object(adata, total_subscan, a1, adata_subset):
    #total_subscan.obs = adata_subset.obs
    total_subscan.obsm = adata_subset.obsm
    total_subscan.obsp = adata_subset.obsp
    total_subscan.uns = adata_subset.uns
    sc.tl.pca(adata)
    total_subscan.uns['pca'] = adata.uns['pca']
    total_subscan.uns['umap'] = {'params': {'a': 1.8956058664239412, 'b': 0.8006378441176886}}
    total_subscan.obsm['X_old_umap'] = total_subscan.obsm['X_umap']
    total_subscan.obsm['X_umap'] = total_subscan.obsm['X_wnn.umap']
    total_subscan.obs['batch'] = ['scrna' for i in range(len(total_subscan.obs.index.tolist()))]
    sc.tl.pca(total_subscan)
    a1.obs['batch'] = ['merfish' for i in range(len(a1.obs.index))]
    return adata, total_subscan, a1, adata_subset

def find_anchors_and_transfer_labels(alist):
    integrator = hanqing_integrate.SeuratIntegration()
    anchor = integrator.find_anchor(alist,
                                    k_local=None,
                                    key_local='X_pca',
                                    k_anchor=5,
                                    key_anchor='X',
                                    dim_red='cca',
                                    max_cc_cells=100000,
                                    k_score=30,
                                    k_filter=None,
                                    scale1=False,
                                    scale2=False,
                                    n_components=50,
                                    n_features=200,
                                    alignments=[[[0], [1]]])
    corrected = integrator.integrate(key_correct='X_pca',
                                     row_normalize=True,
                                     n_components=30,
                                     k_weight=100,
                                     sd=1,
                                     alignments=[[[0], [1]]])

    corrected = np.concatenate(corrected)


    cells = sum([a.shape[0] for a in alist])
    features = alist[0].shape[1]

    try_avr = anndata.AnnData(X=csr_matrix(([], ([], [])),
                                               shape=(cells, features)),
                                  obs=pd.concat([a.obs for a in alist]),
                                  var=alist[0].var)

    try_avr.obsm['X_pca_integrate'] = corrected

    categorical_keys = []
    continuous_keys = []
    for i in alist[0].obs.columns.tolist():
        #print(type(total_subscan.obs[i].tolist()[0]))
        if type(alist[0].obs[i].tolist()[0]) == str:
            categorical_keys.append(i)
        else:
            continuous_keys.append(i)
    
    copy_try = try_avr.copy()
    transfer_results = integrator.label_transfer(
        ref=[0],
        qry=[1],
        categorical_key=categorical_keys,
        continuous_key=continuous_keys,
        key_dist='X_pca'
    )

    integrator.save_transfer_results_to_adata(copy_try, transfer_results)

    replicate = alist[1].copy()
    for i in categorical_keys:
        replicate.obs[i + '_transfer'] = copy_try[copy_try.obs.index.isin(set(replicate.obs.index)), :].obs[i+'_transfer'].tolist()
    
    for i in continuous_keys:
        replicate.obs[i + '_transfer'] = transfer_results[i]
        
    return copy_try, replicate

def harmony_on_merfish(try_avr):
    ho = run_harmony(data_mat=try_avr.obsm['X_pca_integrate'],
                     meta_data=try_avr.obs,
                     nclust=50,
                     vars_use=['batch'], 
                     max_iter_harmony=50)
    try_avr.obsm['X_harmony'] = ho.Z_corr.T
    try_avr.obsm['X_pca'] = try_avr.obsm['X_harmony']

    sc.pp.neighbors(try_avr)  
    return try_avr

def plot_umap_merfish_rna(try_avr, experiment_name):
    min_dist = max(0.1, 1 - try_avr.shape[0] / 60000)
    sc.tl.umap(try_avr, min_dist = min_dist, spread = 1)
    sc.pl.umap(try_avr, color = 'batch', title = experiment_name + ' scRNA Joint Embedding - Batch', save = experiment_name+'_batch.png')
    sc.pl.umap(try_avr, color = 'seurat_clusters_transfer', title = experiment_name +' scRNA Joint Embedding - Clusters', palette = 'tab20', show=False)
    try_avr.uns['seurat_clusters_transfer_colors'][8] = 'k'
    try_avr.uns['seurat_clusters_transfer_colors'][0] = 'r'
    try_avr.uns['seurat_clusters_transfer_colors'][3] = 'b'
    try_avr.uns['seurat_clusters_transfer_colors'][1] = 'tab:orange'
    try_avr.uns['seurat_clusters_transfer_colors'][6] = 'tab:purple'
    try_avr.uns['seurat_clusters_transfer_colors'][10] = 'y'
    sc.pl.umap(try_avr, color = 'seurat_clusters_transfer', title = experiment_name+' scRNA Joint Embedding - Clusters', save = experiment_name+'_cluster.png')
    return try_avr

def move_figure_path(lfile, experiment_name):
    source = r'C:\Users\amonell\Desktop\Scripts_Alex\MERSCOPE_PIPELINE\figures'

    destination = os.path.dirname(lfile)
    try:
        os.mkdir(destination + os.path.sep + 'cluster_images')
    except:
        None
    allfiles = ['umap'+experiment_name+'_cluster.png', 'umap'+experiment_name+'_batch.png']
    for f in allfiles:
        src_path = os.path.join(source, f)
        dst_path = os.path.join(destination, f)
        shutil.move(src_path, dst_path)
    return destination
        
def plot_clusters(a1, try_avr, experiment_name, lfile, destination, save = True):
    a1.obs['seurat_clusters'] = try_avr[try_avr.obs['experiment'] == experiment_name, :].obs['seurat_clusters_transfer'].tolist()
    sc.set_figure_params(dpi=75, dpi_save=300, facecolor='white')
    cluster_types = np.unique(a1.obs['seurat_clusters'])
    source = r'C:\Users\amonell\Desktop\Scripts_Alex\MERSCOPE_PIPELINE\figures'

    for i in [experiment_name]:
        for j in range(len(cluster_types)):
            embed_obj = a1[a1.obs['experiment'] == i, :]
            embed_obj.obs['color_new'] = [1 if k == cluster_types[j] else 0 for k in embed_obj.obs['seurat_clusters']]
            if save == True:
                sc.pl.embedding(embed_obj, basis='spatial', color= 'color_new', cmap = 'Reds', s = 5, title = i.capitalize() + ' Cluster ' + cluster_types[j], colorbar_loc=None, vmax = 1, vmin=-0.2, show=False, save = r'Cluster ' + cluster_types[j] + '-' +i.split('_')[0].capitalize())
                src_path = os.path.join(source, r'spatialCluster ' + cluster_types[j] + '-' +i.split('_')[0].capitalize()+'.pdf')
                dst_path = os.path.join(destination + os.path.sep + 'cluster_images', r'spatialCluster ' + cluster_types[j] + '-' +i.split('_')[0].capitalize()+'.pdf')
                shutil.move(src_path, dst_path)
            else:
                sc.pl.embedding(embed_obj, basis='spatial', color= 'color_new', cmap = 'Reds', show=True, s = 5, title = i.capitalize() + ' Cluster ' + cluster_types[j], colorbar_loc=None, vmax = 1, vmin=-0.2)
                
def plot_motifs(partial_ad, motif_cols, destination):
    source = r'C:\Users\amonell\Desktop\Scripts_Alex\MERSCOPE_PIPELINE\figures'
    try:
        os.mkdir(destination + os.path.sep + 'motif_images')
    except:
        None
    save = True
    for mot in motif_cols:
        if save == True:
            sc.pl.embedding(partial_ad, basis='spatial', color = mot, cmap = 'coolwarm', title = mot.split('_transfer')[0], show = False, save = mot.split('_transfer')[0].replace('/', '').replace('.', '').replace(':', '')+'.png')
            src_path = os.path.join(source, 'spatial'+mot.split('_transfer')[0].replace('/', '').replace('.', '').replace(':', '')+'.png')
            dst_path = os.path.join(destination + os.path.sep + 'motif_images', mot.split('_transfer')[0].replace('/', '').replace('.', '').replace(':', '')+'.png')
            shutil.move(src_path, dst_path)
        else:
            sc.pl.embedding(partial_ad, basis='spatial', color = mot, cmap = 'coolwarm', title = mot.split('_transfer')[0], show = True)

def probe_r_or_scanpy():
    messagebox.showinfo("Locate","Please locate the scanpy or h5 format Seurat object after clicking ok")
    should_load = messagebox.askyesno("Do you want to load the object at path C:/Users/amonell/Downloads/combined_filtered.h5")
    if should_load == True:
        adata_file = 'C:/Users/amonell/Downloads/combined_filtered.h5'
    else:
        root = Tk()
        root.withdraw()
        adata_file = askopenfilename()
    answer = simpledialog.askstring("Input", "Is the single cell object R or scanpy? (type R or scanpy)")
    if (answer == 'R') or (answer == 'r'):
        adata = diopy.input.read_h5(file = adata_file)
    else:
        adata = sc.read(adata_file)
    return adata