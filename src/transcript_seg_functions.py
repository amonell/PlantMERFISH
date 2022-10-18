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
from rtree import index
import math
import random

class MERSCOPE_analyzer:
    def __init__(self, detected_transcripts, output_folder, merscope_region_folder, transcripts_interest):
        self.detected_transcripts = detected_transcripts
        self.output_folder = output_folder
        self.merscope_region_folder = merscope_region_folder
        self.transcripts_interest = transcripts_interest
        self.baysor = False
        
    def get_image(self):      
        dict_detected = {}

        xlist = self.transcripts_interest['global_x'].astype(int).tolist()
        ylist = self.transcripts_interest['global_y'].astype(int).tolist()

        for i in range(len(xlist)):
            dict_detected[xlist[i]] = []
        for i in range(len(ylist)):
            dict_detected[xlist[i]].append(ylist[i])

        #create image from transcript locations
        new_img = np.zeros((int(np.max(self.transcripts_interest['global_x']+3000)), int(np.max(self.transcripts_interest['global_y']+3000))))
        for i in dict_detected.keys():
            new_img[i][dict_detected.get(i)] = 1

        self.dict_detected = dict_detected
        self.new_img = new_img
        
    def run_cellpose_custom(self):
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
        for i in range(2000, len(self.new_img[0]), 2000):
            for j in range(2000, len(self.new_img[1]), 2000):      
                ksize = (5, 5)

                # Using cv2.blur() method 
                image = cv2.blur(self.new_img[i-2000:i, j-2000:j], ksize) 
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
        for i in range(2000, len(self.new_img[0]), 2000):
            for j in range(2000, len(self.new_img[1]), 2000):      
                ksize = (5, 5)

                # Using cv2.blur() method 
                image = cv2.blur(self.new_img[i-2000:i, j-2000:j], ksize) 
                image = image * 255
                temp_img_arr.append(image)
        io.save_masks(temp_img_arr, masks, flows, file_names=['stack_prestain_'+str(i).zfill(6)+'_cp_masks'+'.png' for i in range(len(temp_img_arr))], savedir = self.output_folder + os.path.sep +'cellpose_predictions')
        try:
            os.makedirs(self.output_folder + os.path.sep +'images_for_cellpose_prediction')
        except:
            None
        for i in range(len(temp_img_arr)):   
            io.imsave(self.output_folder + os.path.sep +'images_for_cellpose_prediction'+os.path.sep+'image-'+str(i)+'.tiff', temp_img_arr[i])
        self.masks = masks 
        
    def create_total_image_from_masks(self):
        from PIL import Image
        for i in range(len(self.masks)):
            self.masks[i] = self.masks[i]+(i*10000)
        h_list = []
        for j in range(len(self.new_img[1])//2000):
            h_list.append(np.hstack([self.masks[i] for i in range(j*(len(self.new_img[0])//2000),  j*(len(self.new_img[0])//2000)+ len(self.new_img[0])//2000)]))
        self.reconstruction = np.vstack(h_list)

        masks_read = []
        for i in glob.glob(self.output_folder + os.path.sep +'cellpose_predictions'+ os.path.sep+ '*.png'):
            masks_read.append(np.array(Image.open(i)))

        h_list = []
        for j in range(len(self.new_img[1])//2000):
            h_list.append(np.hstack([masks_read[i] for i in range(j*(len(self.new_img[0])//2000),  j*(len(self.new_img[0])//2000) + len(self.new_img[0])//2000)]))
        self.vizualized_reconstruction = np.vstack(h_list)
        io.imsave(self.output_folder + os.path.sep +'cellmask_reconstruction.tiff', self.vizualized_reconstruction)
        self.vizualized_reconstruction = np.clip(self.vizualized_reconstruction, 0, 1)

    def top3(self):
        adata = sc.read('DataPathogenPanel1/combined_filtered.h5')
        gene_to_id_table = pd.read_csv('DataPathogenPanel1/geneID_to_geneName_MERSCOPE_panel1.txt', sep='\t', index_col=0)
        should_explore = messagebox.askyesno('Save all top 3 gene images?', 'Should the top 3 DE gene per cluster transcript images be saved? This may take a while')
        if should_explore == True:
            try:
                os.mkdir(self.output_folder + os.path.sep + 'Transcript_ClusterTop3_Images')
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
                plt.title('Experiment '+ os.path.basename(self.output_folder)+': Gene - '+gene_name+', Marker for Cluster '+str(cluster))
                plt.savefig(self.output_folder + os.path.sep + 'Transcript_ClusterTop3_Images'+os.path.sep + os.path.basename(self.output_folder)+' - Gene - '+gene_name.split('.')[0].replace('/', '')+', Marker for Cluster '+str(cluster))
                plt.close()

            for i in range(len(top3_list)):
                for j in range(len(top3_list[0])):
                    print(top3_list[i][j])
                    display_gene_spatial(self.vizualized_reconstruction, top3_list[i][j], self.detected_transcripts, i) 
                    
    def get_baysor_commands(self):
        do_bays = messagebox.askyesno('Run Baysor', 'Do you want to run Baysor, or only the Cellpose segmentation?')
        if do_bays == True:
            num_cells_total = len(np.unique(self.reconstruction))
            print('/mnt variable set to D:/Tatsuya/merscope', end='\n')
            try:
                os.mkdir('D:/Tatsuya/merscope/Baysor_segmentations/'+os.path.basename(self.output_folder))
            except:
                None            
            print('baysor run -s 250 -x global_x -y global_y -z global_z -o segmentation_'+os.path.basename(self.output_folder)+'.csv -g gene --num-cells-init '+str(int(num_cells_total*1.1))+' --n-clusters 1 --force-2d -i 1 -c /mnt/Baysor_segmentations/baysor_avr_fullrun_config.toml /mnt/'+self.merscope_region_folder.replace('D:/Tatsuya/merscope/', '')+'/detected_transcripts.csv', end = '\n')
            print('cp *'+os.path.basename(self.output_folder)+'* /mnt/Baysor_segmentations/'+os.path.basename(self.output_folder))
            self.baysor = True
        else:
            return

    def create_cell_by_gene_(self):
        gene_list = np.unique(self.detected_transcripts['gene'])
        cell_list = np.unique(self.reconstruction, return_counts=True)[0]
        cell_by_gene = pd.DataFrame(columns=gene_list, index=cell_list.astype(int))
        cell_by_gene = cell_by_gene.fillna(0)
        postions_not_in = []
        for i in tqdm(gene_list):
            hold_part = self.detected_tanscripts[self.detected_transcripts['gene'] == i]
            x_part = hold_part['global_x'].astype(int).tolist()
            y_part = hold_part['global_y'].astype(int).tolist()   
            for k in range(len(x_part)):
                try:
                    value_at = self.reconstruction[x_part[k]][y_part[k]]
                    self.cell_by_gene.loc[value_at][i] += 1 
                except:
                    postions_not_in.append((x_part,y_part, i))
        self.cell_by_gene.to_csv(self.output_folder + os.path.sep + 'cell_by_gene.csv', index=False)

    def create_adata_and_add_observations(self):
        self.adata = sc.read(self.output_folder + os.path.sep + 'cell_by_gene.csv')
        self.adata.obs.index = self.cell_by_gene.index
        total_unique = []
        for i in range(len(self.masks)):
        #     print(np.unique(masks[i], return_counts=True)[0])
            [total_unique.append(j) for j in np.unique(self.masks[i], return_counts=True)[1]]
        self.adata.obs['area'] = total_unique
        unique_cells_index = np.unique(self.reconstruction, return_index=True)[1]
        spatial = []
        for i in range(len(unique_cells_index)):
            spatial.append([unique_cells_index[i] % len(self.reconstruction[0]), unique_cells_index[i] // len(self.reconstruction[0])])
        self.adata.obsm['X_spatial'] = np.array(spatial)
        self.adata = self.adata[self.adata.obs.index.astype(int) % 10000 != 0]
        self.num_transcripts = list(np.sum(self.adata.X, axis=1))
        self.self.adata.obs['transcripts_per_cell'] = self.num_transcripts
        total_gene = list(np.sum(self.adata.X, axis=0))
        self.adata.var['total_gene_counts'] = total_gene
        sc.pp.calculate_qc_metrics(self.adata, inplace=True)
    
    def plot_qc(self):
        ax = sc.pl.violin(self.adata, 'n_genes_by_counts', jitter=0.4, multi_panel=False, show=False)
        ax.set_title('Counts per Gene, Average = '+str(np.round(np.mean(self.adata.obs['n_genes_by_counts']), 2)))
        ax.set_ylabel('Number of Transcripts')
        ax.set_xlabel('Gene')
        try: 
            os.mkdir('figures')
        except:
            None
        plt.savefig('figures/violin_quality_metrics_cpg.png')
        plt.close()
        ax = sc.pl.violin(self.adata, 'total_counts', jitter=0.4, multi_panel=False, show=False)
        ax.set_title('Total Gene Counts Per Cell, Average = '+str(np.round(np.mean(self.adata.obs['total_counts']), 2)))
        ax.set_ylabel('Number of Transcripts')
        ax.set_xlabel('Cell')
        plt.savefig('figures/violin_quality_metrics_tgc.png')
        plt.close()

    def postprocess_adata(self):

        #normalization
        for i in range(len(self.adata.X)):
            #by area
            #adata.X[i] /= total_unique[i]
            #by tpc
            self.adata.X[i] /= self.num_transcripts[i]
        self.adata = self.adata[self.adata.obs['transcripts_per_cell'] >= 50]
        sc.pp.log1p(self.adata)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata)
        sc.tl.leiden(self.adata)
        sc.tl.umap(self.adata)

        sc.pl.umap(self.adata, color='leiden', vmax=0.05, save = '_leiden.png', show = False)
        sc.pl.embedding(self.adata, basis='spatial', color='leiden', vmax=0.05, size=10, save = '_leiden.png', show = False)
        self.adata.write(self.output_folder + os.path.sep + 'adata.h5ad')
        source = r'figures'
        destination = self.output_folder
        allfiles = ['umap_leiden.png', 'spatial_leiden.png', 'violin_quality_metrics_cpg.png', 'violin_quality_metrics_tgc.png']
        for f in allfiles:
            src_path = os.path.join(source, f)
            dst_path = os.path.join(destination, f)
            shutil.move(src_path, dst_path)

    def read_in_segmentation_baysor(self):
        self.baysor_segmentation_folder = r'D:\Tatsuya\merscope\Baysor_segmentations'+os.path.sep+os.path.basename(self.output_folder)
        self.segmentation_name = 'segmentation_'+os.path.basename(self.output_folder)
        self.segmentation = pd.read_csv(os.path.join(self.baysor_segmentation_folder, self.segmentation_name + '.csv'))
        self.segmentation_cell_stats = pd.read_csv(os.path.join(self.baysor_segmentation_folder, self.segmentation_name + '_cell_stats.csv'))
        xvals = self.segmentation_cell_stats['x'].tolist()
        yvals = self.segmentation_cell_stats['y'].tolist()
        color = self.segmentation_cell_stats['n_transcripts'].tolist()
        plt.scatter(xvals, yvals, c=color, s = 1)
        plt.savefig(os.path.join(self.baysor_segmentation_folder,'initial_seg_mosaic.png'))
        plt.close()

    def create_baysor_cellxgene(self):    
        whole_gene_counts = []
        gene_list = np.unique(self.segmentation['gene'].tolist()).astype(str)
        self.segmentation_cell_stats_subset = self.segmentation_cell_stats[self.segmentation_cell_stats['n_transcripts'] > 50] 
        subset_segmentation = self.segmentation[self.segmentation['cell'].isin(self.segmentation_cell_stats_subset['cell'])]
        for i in tqdm(gene_list):
            temp = subset_segmentation[subset_segmentation['gene'] == i]
            temp_df = pd.DataFrame(np.array(np.unique(temp['cell'].tolist(), return_counts=True)).T)
            gene_counts = pd.merge(self.segmentation_cell_stats_subset, temp_df, left_on= 'cell', right_on=0, how = 'outer')[1]
            gene_counts = gene_counts.fillna(0)
            whole_gene_counts.append(gene_counts.astype(int).tolist())
        self.df_cell_by_gene = pd.DataFrame(np.array(whole_gene_counts))
        self.df_cell_by_gene.columns = self.segmentation_cell_stats_subset['cell'].tolist()
        self.df_cell_by_gene.index = gene_list
        self.df_cell_by_gene.T.to_csv(os.path.join(self.baysor_segmentation_folder, 'cell_by_gene.csv'), index=None)

    def create_baysor_adata(self):
        self.adata = sc.read(os.path.join(self.baysor_segmentation_folder, 'cell_by_gene.csv'))
        self.adata.obs.index = self.segmentation_cell_stats_subset['cell'].tolist()
        self.adata = self.adata[self.adata.obs.index.isin(set(self.segmentation_cell_stats['cell'].tolist()))]
        self.adata.obs = pd.merge(self.adata.obs, self.segmentation_cell_stats, left_index=True, right_on='cell')
        self.adata.obs.index = self.adata.obs.cell
        self.adata.obs = self.adata.obs.drop(['cell'], axis = 1)
        xspatial = self.adata.obs['x'].tolist()
        yspatial = self.adata.obs['y'].tolist()
        self.adata.obsm['X_spatial'] = np.array((xspatial, yspatial)).T
        self.num_transcripts = list(np.sum(self.adata.X, axis=1))
        self.adata.obs['transcripts_per_cell'] = self.num_transcripts
        total_gene = list(np.sum(self.adata.X, axis=0))
        self.adata.var['total_gene_counts'] = total_gene
        sc.pp.calculate_qc_metrics(self.adata, inplace=True)

    def baysor_qc_plots(self):
        ax = sc.pl.violin(self.adata, 'n_genes_by_counts', jitter=0.4, multi_panel=False, show=False)
        ax.set_title('Different Genes per Cell, Average = '+str(np.round(np.mean(self.adata.obs['n_genes_by_counts']), 2)))
        ax.set_ylabel('Number of Genes')
        ax.set_xlabel('Cells')
        try: 
            os.mkdir('figures')
        except:
            None

        plt.savefig(os.path.join(self.baysor_segmentation_folder, 'violin_quality_metrics_cpg.png'))
        ax.set_xticks([])
        plt.close()
        ax = sc.pl.violin(self.adata, 'total_counts', jitter=0.4, multi_panel=False, show=False)
        ax.set_title('Total Gene Counts Per Cell, Average = '+str(np.round(np.mean(self.adata.obs['total_counts']), 2)))
        ax.set_ylabel('Number of Transcripts')
        ax.set_xlabel('Cells')
        ax.set_xticks([])
        plt.savefig(os.path.join(self.baysor_segmentation_folder,'violin_quality_metrics_tgc.png'))
        plt.close()

        ax = sc.pl.violin(self.adata, 'area', jitter=0.4, multi_panel=False, show=False)
        ax.set_title('Area Per Cell, Average = '+str(np.round(np.mean(self.adata.obs['area']), 2)))
        ax.set_ylabel('Area')
        ax.set_xlabel('Cells')
        ax.set_ylim(0,5000)
        ax.set_xticks([])
        plt.savefig(os.path.join(self.baysor_segmentation_folder,'violin_quality_metrics_area.png'))
        plt.close()

    def baysor_postprocess_adata(self):
        for i in range(len(self.adata.X)):
            #by area
            #adata.X[i] /= total_unique[i]
            #by tpc
            self.adata.X[i] /= self.num_transcripts[i]
        self.adata = self.adata[self.adata.obs['transcripts_per_cell'] >= 50]
        sc.pp.log1p(self.adata)
        sc.tl.pca(self.adata)
        sc.pp.neighbors(self.adata)
        sc.tl.leiden(self.adata)
        sc.tl.umap(self.adata)

        sc.pl.umap(self.adata, color='leiden', vmax=0.05, save = '_leiden.png', show = False)
        sc.pl.embedding(self.adata, basis='spatial', color='leiden', vmax=0.05, size=10, save = '_leiden.png', show = False)
        self.adata.write(self.baysor_segmentation_folder + os.path.sep + 'adata.h5ad')
        source = r'figures'
        destination = self.baysor_segmentation_folder
        allfiles = ['umap_leiden.png', 'spatial_leiden.png']
        for f in allfiles:
            src_path = os.path.join(source, f)
            dst_path = os.path.join(destination, f)
            shutil.move(src_path, dst_path)
        self.adata.write(os.path.join(self.baysor_segmentation_folder, 'adata.h5ad'))

def take_user_input():
    messagebox.showinfo("Option","Please locate the region0 folder for the experiment")
    root = Tk()
    root.withdraw()
    folder_selected1 = filedialog.askdirectory()
    messagebox.showinfo("Option","Please select an output folder for the experiment analysis results")
    folder_selected = filedialog.askdirectory()
    question = simpledialog.askstring("Input", "Input segmentation slices or auto-detect? ('y' for input slices)")
    try:
        detected_tanscripts = pd.read_csv(folder_selected1+os.path.sep+'detected_transcripts.csv')
    except:
        print('detected_transcripts.csv not found in the folder')
    if question == 'y':
        while True:
            try:
                answer = simpledialog.askstring("Input", "Slices to run segmentation on (separated by commas)")
                answer = answer.replace(' ', '')
                answer = answer.split(',')
                answer = [float(i) for i in answer]   
                break
            except:
                print('bad input, try again')
                continue
    else:
        answer = [float(p) for p in choose_cellpose_segmentation_slices(detected_tanscripts)]
    

    return folder_selected1, detected_tanscripts, answer, folder_selected

def choose_cellpose_segmentation_slices(detected_transcripts):
    combinations_list = [{3}, {0}, {6}, {3, 4}, {5,6}, {0,1}, {3,4,5}, {0,1,2}, {3,4,5,6}, {1,2,3,4,5}, {1,2,3,4,5,6}, {0,1,2,3,4,5,6}]
    best_value = None
    best_comb = None
    for comb in tqdm(combinations_list):
        if len(comb) > 0:
            subsetted_obj = detected_transcripts[detected_transcripts['global_z'].astype(int).isin(comb)]
            xvalues = np.array(subsetted_obj['global_x'].tolist())
            yvalues = np.array(subsetted_obj['global_y'].tolist())
            lower_size = [i for i in range(0, int(len(xvalues)/200))]
            [lower_size.append(i) for i in range(int(len(xvalues)/2)-int(len(xvalues)/200), int(len(xvalues)/2))]
            [lower_size.append(i) for i in range(int(len(xvalues))-int(len(xvalues)/200), int(len(xvalues)))]

            idx = index.Index()
            for i in np.array(lower_size):
                idx.insert(i, np.array([xvalues[i], yvalues[i]]))
            #stochastically determine density
            randomlist = random.sample(range(0, len(xvalues)), min(int(len(lower_size)/2), 5000))
            avg_distances = []
            for p in range(len(randomlist)):
                all_neighbors = []
                nn_obj = idx.nearest([xvalues[p], yvalues[p]], num_results = 25, objects=True)
                for i in nn_obj:
                    all_neighbors.append([math.sqrt(math.pow((xvalues[p] - i.bounds[0]),2) + math.pow((yvalues[p] - i.bounds[3]),2))])
                avg_distances.append(np.mean(all_neighbors[1:]))
            avg_distances = np.array(avg_distances) 
            total_average = np.median(avg_distances)

            #4 because that was the density the cellpose model was trained on
            if best_value==None or (abs(4 - total_average) < best_value):
                best_value = abs(4 - total_average)
                best_comb = comb
    return list(best_comb)