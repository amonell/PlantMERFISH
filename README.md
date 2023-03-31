# PlantMERFISH

This repository contains code for analyzing Vizgen MERSCOPE data obtained for the publication:
"Time-resolved single-cell and spatial gene regulatory atlas of plants under pathogen attack", Nobori et al.

# Usage
- The Transcript_Segmentation.ipynb notebook takes the data folders output from the MERSCOPE as input, and then prompts the user for inputs to help them get cell segmentation boundaries from only transcripts. If they already have a cell segmentation that they are comfortable with, and have already assigned transcripts to cells to create a Scanpy Anndata object for their experiment, then they do not need this step.

- The Merging_MERFISH_Objects.ipynb notebook takes a dataset containing scRNA-seq data, in this case, a 10x multiome dataset, and integrates it with a MERFISH dataset. This means that a joint UMAP is calculated, and both categorical and continuous observations from the reference dataset are transferred to the MERFISH Anndata.

- 
# Requirements
This code was developed using Python version 3.10.4. You can find the required Python packages in the requirements.txt file.
