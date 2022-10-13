import pandas as pd
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import imageio as io
import scanpy as sc
import anndata as ad
import os
import logging
from cellpose import models, io
import glob
from PIL import Image
import cv2
from tkinter import Tk     # from tkinter import Tk for Python 3.x
from tkinter import filedialog
from tkinter import *
from tkinter import ttk
from tkinter import simpledialog
from tkinter import messagebox
import tkinter as tk
import os
import logging
from cellpose import models, io
import shutil

def get_image(d2):      
    dict_detected = {}

    xlist = d2['global_x'].astype(int).tolist()
    ylist = d2['global_y'].astype(int).tolist()
    print(len(xlist))
    for i in range(len(xlist)):
        dict_detected[xlist[i]] = []
    for i in range(len(ylist)):
        dict_detected[xlist[i]].append(ylist[i])

    #create image from transcript locations
    new_img = np.zeros((int(np.max(d2['global_x']+3000)), int(np.max(d2['global_y']+3000))))
    for i in dict_detected.keys():
        new_img[i][dict_detected.get(i)] = 1
        
    return dict_detected, new_img

def take_user_input():
    messagebox.showinfo("Option","Please locate the region0 folder for the experiment")
    root = Tk()
    root.withdraw()
    folder_selected1 = filedialog.askdirectory()

    messagebox.showinfo("Option","Next question appearing shortly")
    try:
        detected_tanscripts = pd.read_csv(folder_selected1+os.path.sep+'detected_transcripts.csv')
    except:
        print('detected_transcripts.csv not found in the folder')
    answer = simpledialog.askstring("Input", "Slices to run segmentation on (separated by commas)")
    answer = answer.replace(' ', '')
    answer = answer.split(',')
    answer = [float(i) for i in answer]   
    
    messagebox.showinfo("Option","Please select an output folder for the experiment analysis results")
    folder_selected = filedialog.askdirectory()
    return folder_selected1, detected_tanscripts, answer, folder_selected

def run_cellpose_custom(new_img, folder_selected):
    model = models.CellposeModel(gpu=True, pretrained_model='Models/CP_20220920_113001')
    # define CHANNELS to run segementation on
    # grayscale=0, R=1, G=2, B=3
    # channels = [cytoplasm, nucleus]
    # if NUCLEUS channel does not exist, set the second channel to 0
    channels = [0,0]
    masks = []
    flows = []
    styles = []
    diams = []
    ct = 0
    for i in range(2000, len(new_img[0]), 2000):
        for j in range(2000, len(new_img[1]), 2000):      
            ksize = (5, 5)

            # Using cv2.blur() method 
            image = cv2.blur(new_img[i-2000:i, j-2000:j], ksize) 
            image = image * 255

            try:
                assert len(np.unique(image))>1
                masks_, flows_, styles_= model.eval([image], channels=channels, diameter=22.92,flow_threshold=0.7, cellprob_threshold=-2)
                masks.extend(masks_)
                flows.extend(flows_)
                styles.extend(styles_)
                print(ct, end = ' ')
            except (AssertionError):
                try:
                    masks.append(np.zeros((len(image[0]), len(image[1]))))
                    flows.append([])
                    styles.append([])
                except:
                    masks.append([])
                    flows.append([])
                    styles.append([])   
            ct += 1
    temp_img_arr = []
    for i in range(2000, len(new_img[0]), 2000):
        for j in range(2000, len(new_img[1]), 2000):      
            ksize = (5, 5)

            # Using cv2.blur() method 
            image = cv2.blur(new_img[i-2000:i, j-2000:j], ksize) 
            image = image * 255
            temp_img_arr.append(image)
    io.save_masks(temp_img_arr, masks, flows, file_names=['stack_prestain_'+str(i).zfill(6)+'_cp_masks'+'.png' for i in range(len(temp_img_arr))], savedir = folder_selected + os.path.sep +'cellpose_predictions')
    try:
        os.makedirs(folder_selected + os.path.sep +'images_for_cellpose_prediction')
    except:
        None
    for i in range(len(temp_img_arr)):   
        io.imsave(folder_selected + os.path.sep +'images_for_cellpose_prediction'+os.path.sep+'image-'+str(i)+'.tiff', temp_img_arr[i])
    return masks
   
def create_total_image_from_masks(masks, new_img, folder_selected):
    from PIL import Image
    for i in range(len(masks)):
        masks[i] = masks[i]+(i*10000)
    h_list = []
    for j in range(len(new_img[1])//2000):
        h_list.append(np.hstack([masks[i] for i in range(j*(len(new_img[0])//2000),  j*(len(new_img[0])//2000)+ len(new_img[0])//2000)]))
    reconstruction = np.vstack(h_list)
    
    masks_read = []
    for i in glob.glob(folder_selected + os.path.sep +'cellpose_predictions'+ os.path.sep+ '*.png'):
        masks_read.append(np.array(Image.open(i)))

    h_list = []
    for j in range(len(new_img[1])//2000):
        h_list.append(np.hstack([masks_read[i] for i in range(j*(len(new_img[0])//2000),  j*(len(new_img[0])//2000) + len(new_img[0])//2000)]))
    vizualized_reconstruction = np.vstack(h_list)
    io.imsave(folder_selected + os.path.sep +'cellmask_reconstruction.tiff', vizualized_reconstruction)
    vizualized_reconstruction = np.clip(vizualized_reconstruction, 0, 1)
    return reconstruction, vizualized_reconstruction

def top3(folder_selected, detected_tanscripts, vizualized_reconstruction):
    adata = sc.read('DataPathogenPanel1/combined_filtered.h5')
    gene_to_id_table = pd.read_csv('DataPathogenPanel1/geneID_to_geneName_MERSCOPE_panel1.txt', sep='\t', index_col=0)
    should_explore = messagebox.askyesno('Save all top 3 gene images?', 'Should the top 3 DE gene per cluster transcript images be saved? This may take a while')
    if should_explore == True:
        try:
            os.mkdir(folder_selected + os.path.sep + 'Transcript_ClusterTop3_Images')
        except:
            None
        finding_cluster_markers = adata[:,adata.var.index.isin(set(gene_to_id_table['gene_name']))]
        sc.tl.rank_genes_groups(finding_cluster_markers, groupby='seurat_clusters')
        top3_DE = [i for i in finding_cluster_markers.uns['rank_genes_groups']['names']][:3]
        top3_list = []
        for i in range(len(top3_DE[0])):
            top3_list.append([top3_DE[0][i], top3_DE[1][i], top3_DE[2][i]])

        def display_gene_spatial(cell_image, gene_name, detected_transcript_df, cluster):
            toplot = gene_to_id_table[gene_to_id_table['gene_name'] == gene_name]['gene_id']
            xandy = detected_transcript_df[detected_transcript_df['gene'] == toplot.tolist()[0]]
            x = xandy['global_x'].tolist()
            y = xandy['global_y'].tolist()
            plt.figure(figsize=(10, 10), dpi=100)
            plt.ylim(len(cell_image.T),0)
            plt.imshow(cell_image.T, vmax = 2.3, cmap = 'Greys_r')
            plt.scatter(x, y, s = 0.15, color = 'cyan')
            plt.title('Experiment '+ os.path.basename(folder_selected)+': Gene - '+gene_name+', Marker for Cluster '+str(cluster))
            plt.savefig(folder_selected + os.path.sep + 'Transcript_ClusterTop3_Images'+os.path.sep + os.path.basename(folder_selected)+' - Gene - '+gene_name.split('.')[0].replace('/', '')+', Marker for Cluster '+str(cluster))
            plt.close()

        for i in range(len(top3_list)):
            for j in range(len(top3_list[0])):
                print(top3_list[i][j])
                display_gene_spatial(vizualized_reconstruction, top3_list[i][j], detected_tanscripts, i) 

def create_cell_by_gene_(detected_tanscripts, reconstruction, folder_selected):
    gene_list = np.unique(detected_tanscripts['gene'])
    cell_list = np.unique(reconstruction, return_counts=True)[0]
    cell_by_gene = pd.DataFrame(columns=gene_list, index=cell_list.astype(int))
    cell_by_gene = cell_by_gene.fillna(0)
    postions_not_in = []
    for i in tqdm(gene_list):
        hold_part = detected_tanscripts[detected_tanscripts['gene'] == i]
        x_part = hold_part['global_x'].astype(int).tolist()
        y_part = hold_part['global_y'].astype(int).tolist()   
        for k in range(len(x_part)):
            try:
                value_at = reconstruction[x_part[k]][y_part[k]]
                cell_by_gene.loc[value_at][i] += 1 
            except:
                postions_not_in.append((x_part,y_part, i))
    cell_by_gene.to_csv(folder_selected + os.path.sep + 'cell_by_gene.csv', index=False)
    return cell_by_gene, postions_not_in

def create_adata_and_add_observations(folder_selected, cell_by_gene, masks, reconstruction):
    adata = sc.read(folder_selected + os.path.sep + 'cell_by_gene.csv')
    adata.obs.index = cell_by_gene.index
    total_unique = []
    for i in range(len(masks)):
    #     print(np.unique(masks[i], return_counts=True)[0])
        [total_unique.append(j) for j in np.unique(masks[i], return_counts=True)[1]]
    adata.obs['area'] = total_unique
    unique_cells_index = np.unique(reconstruction, return_index=True)[1]
    spatial = []
    for i in range(len(unique_cells_index)):
        spatial.append([unique_cells_index[i] % len(reconstruction[0]), unique_cells_index[i] // len(reconstruction[0])])
    adata.obsm['X_spatial'] = np.array(spatial)
    adata = adata[adata.obs.index.astype(int) % 10000 != 0]
    num_transcripts = list(np.sum(adata.X, axis=1))
    adata.obs['transcripts_per_cell'] = num_transcripts
    total_gene = list(np.sum(adata.X, axis=0))
    adata.var['total_gene_counts'] = total_gene
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    return adata, num_transcripts

def plot_qc(adata):
    ax = sc.pl.violin(adata, 'n_genes_by_counts', jitter=0.4, multi_panel=False, show=False)
    ax.set_title('Counts per Gene, Average = '+str(np.round(np.mean(adata.obs['n_genes_by_counts']), 2)))
    ax.set_ylabel('Number of Transcripts')
    ax.set_xlabel('Gene')
    try: 
        os.mkdir('figures')
    except:
        None
    plt.savefig('figures/violin_quality_metrics_cpg.png')
    plt.close()
    ax = sc.pl.violin(adata, 'total_counts', jitter=0.4, multi_panel=False, show=False)
    ax.set_title('Total Gene Counts Per Cell, Average = '+str(np.round(np.mean(adata.obs['total_counts']), 2)))
    ax.set_ylabel('Number of Transcripts')
    ax.set_xlabel('Cell')
    plt.savefig('figures/violin_quality_metrics_tgc.png')
    plt.close()
def postprocess_adata(adata, num_transcripts, folder_selected):

    #normalization
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
    adata.write(folder_selected + os.path.sep + 'adata.h5ad')
    source = r'figures'
    destination = folder_selected
    allfiles = ['umap_leiden.png', 'spatial_leiden.png', 'violin_quality_metrics_cpg.png', 'violin_quality_metrics_tgc.png']
    for f in allfiles:
        src_path = os.path.join(source, f)
        dst_path = os.path.join(destination, f)
        shutil.move(src_path, dst_path)
