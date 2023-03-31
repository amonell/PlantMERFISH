# PlantMERFISH

This repository contains code for analyzing Vizgen MERSCOPE data obtained for the publication:
"Time-resolved single-cell and spatial gene regulatory atlas of plants under pathogen attack", Nobori et al.

# Usage
- The Transcript_Segmentation.ipynb notebook takes the data folders output from the MERSCOPE as input, and then prompts the user for inputs to help them get cell segmentation boundaries from only transcripts. If they already have a cell segmentation that they are comfortable with, and have already assigned transcripts to cells to create a Scanpy Anndata object for their experiment, then they do not need this step.

- The Merging_MERFISH_Objects.ipynb notebook takes a dataset containing scRNA-seq data, in this case, a 10x multiome dataset, and integrates it with a MERFISH dataset. This means that a joint UMAP is calculated, and both categorical and continuous observations from the reference dataset are transferred to the MERFISH Anndata.

- The Psuedotime_Multiome folder contains the code used to perform pseudotime with Palantir on the 10x Multiome Dataset. Instructions on how to create the necessary environment are supplied in the folder as well.

- The smFISH_Quant folder contains the notebook used to identify bacterial colony locations and quantify genes images with smFISH in the MERSCOPE. Instructions on how to create the necessary environment to run this code are supplied in the folder as well.

- load_obj.R can be used to convert a Seurat experiment object to a Scanpy Anndata object.

- Get_Transcript_Images.ipynb can be used to get images of different RNA species across a tissue section

# Requirements
This code was developed using Python version 3.10.4. You can find the required Python packages for the general MERFISH data processing in the requirements.txt file.
