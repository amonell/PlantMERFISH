import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shutil, os
from tqdm.notebook import tqdm
import scanpy as sc

def read_in_segmentation_baysor(baysor_segmentation_folder, segmentation_name):
    segmentation = pd.read_csv(os.path.join(baysor_segmentation_folder, segmentation_name + '.csv'))
    segmentation_cell_stats = pd.read_csv(os.path.join(baysor_segmentation_folder, segmentation_name + '_cell_stats.csv'))
    xvals = segmentation_cell_stats['x'].tolist()
    yvals = segmentation_cell_stats['y'].tolist()
    color = segmentation_cell_stats['n_transcripts'].tolist()
    plt.scatter(xvals, yvals, c=color, s = 1)
    plt.savefig(os.path.join(baysor_segmentation_folder,'initial_seg_mosaic.png'))
    plt.close()
    return segmentation, segmentation_cell_stats

def create_baysor_cellxgene(segmentation, segmentation_cell_stats, baysor_segmentation_folder):    
    whole_gene_counts = []
    gene_list = np.unique(segmentation['gene'].tolist()).astype(str)
    segmentation_cell_stats_subset = segmentation_cell_stats[segmentation_cell_stats['n_transcripts'] > 50] 
    subset_segmentation = segmentation[segmentation['cell'].isin(segmentation_cell_stats_subset['cell'])]
    for i in tqdm(gene_list):
        temp = subset_segmentation[subset_segmentation['gene'] == i]
        temp_df = pd.DataFrame(np.array(np.unique(temp['cell'].tolist(), return_counts=True)).T)
        gene_counts = pd.merge(segmentation_cell_stats_subset, temp_df, left_on= 'cell', right_on=0, how = 'outer')[1]
        gene_counts = gene_counts.fillna(0)
        whole_gene_counts.append(gene_counts.astype(int).tolist())
    df_cell_by_gene = pd.DataFrame(np.array(whole_gene_counts))
    df_cell_by_gene.columns = segmentation_cell_stats_subset['cell'].tolist()
    df_cell_by_gene.index = gene_list
    df_cell_by_gene.T.to_csv(os.path.join(baysor_segmentation_folder, 'cell_by_gene.csv'), index=None)
    return df_cell_by_gene, gene_list, segmentation_cell_stats_subset

def create_baysor_adata(baysor_segmentation_folder, segmentation_cell_stats_subset, segmentation_cell_stats):
    adata = sc.read(os.path.join(baysor_segmentation_folder, 'cell_by_gene.csv'))
    adata.obs.index = segmentation_cell_stats_subset['cell'].tolist()
    adata = adata[adata.obs.index.isin(set(segmentation_cell_stats['cell'].tolist()))]
    adata.obs = pd.merge(adata.obs, segmentation_cell_stats, left_index=True, right_on='cell')
    adata.obs.index = adata.obs.cell
    adata.obs = adata.obs.drop(['cell'], axis = 1)
    xspatial = adata.obs['x'].tolist()
    yspatial = adata.obs['y'].tolist()
    adata.obsm['X_spatial'] = np.array((xspatial, yspatial)).T
    num_transcripts = list(np.sum(adata.X, axis=1))
    adata.obs['transcripts_per_cell'] = num_transcripts
    total_gene = list(np.sum(adata.X, axis=0))
    adata.var['total_gene_counts'] = total_gene
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    return adata, num_transcripts

def baysor_qc_plots(adata, baysor_segmentation_folder):
    ax = sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, multi_panel=False, show=False)
    ax.set_title('Different Genes per Cell, Average = '+str(np.round(np.mean(adata.obs['n_genes_by_counts']), 2)))
    ax.set_ylabel('Number of Genes')
    ax.set_xlabel('Cells')
    try: 
        os.mkdir('figures')
    except:
        None

    plt.savefig(os.path.join(baysor_segmentation_folder, 'violin_quality_metrics_cpg.png'))
    ax.set_xticks([])
    plt.close()
    ax = sc.pl.violin(adata, 'total_counts', jitter=0.4, multi_panel=False, show=False)
    ax.set_title('Total Gene Counts Per Cell, Average = '+str(np.round(np.mean(adata.obs['total_counts']), 2)))
    ax.set_ylabel('Number of Transcripts')
    ax.set_xlabel('Cells')
    ax.set_xticks([])
    plt.savefig(os.path.join(baysor_segmentation_folder,'violin_quality_metrics_tgc.png'))
    plt.close()

    ax = sc.pl.violin(adata, 'area', jitter=0.4, multi_panel=False, show=False)
    ax.set_title('Area Per Cell, Average = '+str(np.round(np.mean(adata.obs['area']), 2)))
    ax.set_ylabel('Area')
    ax.set_xlabel('Cells')
    ax.set_ylim(0,5000)
    ax.set_xticks([])
    plt.savefig(os.path.join(baysor_segmentation_folder,'violin_quality_metrics_tgc.png'))
    plt.close()
    
def baysor_postprocess_adata(adata, num_transcripts, baysor_segmentation_folder):
    for i in range(len(adata.X)):
        #by area
        #adata.X[i] /= total_unique[i]
        #by tpc
        adata.X[i] /= num_transcripts[i]
    adata = adata[adata.obs['transcripts_per_cell'] >= 50]
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    sc.tl.umap(adata)

    sc.pl.umap(adata, color='leiden', vmax=0.05, save = '_leiden.png', show = False)
    sc.pl.embedding(adata, basis='spatial', color='leiden', vmax=0.05, size=10, save = '_leiden.png', show = False)
    adata.write(baysor_segmentation_folder + os.path.sep + 'adata.h5ad')
    source = r'figures'
    destination = baysor_segmentation_folder
    allfiles = ['umap_leiden.png', 'spatial_leiden.png']
    for f in allfiles:
        src_path = os.path.join(source, f)
        dst_path = os.path.join(destination, f)
        shutil.move(src_path, dst_path)
    return adata