import subprocess
import argparse
import os
import re
import pandas as pd
import scanpy as sc
import scvi
import anndata
import scipy.io as sio
import copy

def load_anndata_from_folder_tag(folder):
	#######
	os.chdir(folder)
	files_in_folder = os.listdir(folder)
	file_pattern_2 = r'__.+expression_matrix.mtx'
	#######
	matching_elements_2 = [element for element in files_in_folder if re.search(file_pattern_2, element)]
	#######
	new_string = re.sub(r'expression_matrix.mtx', "", matching_elements_2[0])
	#######
	return new_string


def load_anndata_from_folder(folder):
	os.chdir(folder)
	files_in_folder = os.listdir(folder)
	#######
	file_pattern_1 = r'__.+cell_metadata.csv'
	file_pattern_2 = r'__.+expression_matrix.mtx'
	file_pattern_3 = r'__.+gene_metadata.csv'
	#######
	matching_elements_1 = [element for element in files_in_folder if re.search(file_pattern_1, element)]
	matching_elements_2 = [element for element in files_in_folder if re.search(file_pattern_2, element)]
	matching_elements_3 = [element for element in files_in_folder if re.search(file_pattern_3, element)]
	#######
	expression_matrix = anndata.read_mtx(matching_elements_2[0])
	cell_metadata = pd.read_csv(matching_elements_1[0], index_col=0)
	gene_metadata = pd.read_csv(matching_elements_3[0], index_col=0)
	#######
	anndata_obj = copy.deepcopy(expression_matrix)
	anndata_obj.obs = cell_metadata
	anndata_obj.var['gene_name'] = gene_metadata['x'].tolist()
	######
	return anndata_obj


def Process_anndatasets(anndata_obj,anndata_obj_name):
	### Get vae and lave ###
	anndata_obj.layers["counts"] = anndata_obj.X
	###
	sc.pp.highly_variable_genes(
        anndata_obj,
        flavor="seurat_v3",
        n_top_genes=3000,
        layer="counts",
        batch_key="batch",
        subset=True
	)
	###
	scvi.model.SCVI.setup_anndata(anndata_obj, layer="counts", batch_key="batch")
	vae = scvi.model.SCVI(anndata_obj, n_layers=2, n_latent=30, gene_likelihood="nb")
	vae.train()
	###
	lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=anndata_obj,
        labels_key="celltype",
        unlabeled_category="Unknown",
	)
	###
	lvae.train(max_epochs=20, n_samples_per_label=100)
	###
	anndata_obj.obsm["X_scANVI"] = lvae.get_latent_representation(anndata_obj)
	###
	### next perform cluster analysis ###
	###
	from scvi.model.utils import mde
	anndata_obj.obsm["X_mde_scanvi"] = mde(anndata_obj.obsm["X_scANVI"])
	###
	### output the celltype,batch and cluster plot !!! #######
	###
	FN1 = anndata_obj_name + '_batch.png'
	sc.pl.embedding(
    	anndata_obj,
    	basis="X_mde_scanvi",
    	color="batch",
    	frameon=False,
    	save=FN1
	)
	###
	FN2 = anndata_obj_name + '_celltype.png'
	sc.pl.embedding(
    	anndata_obj,
    	basis="X_mde_scanvi",
    	color="celltype",
    	frameon=False,
    	save=FN2
	)
	###
	sc.pp.neighbors(anndata_obj, use_rep="X_scANVI")
	sc.tl.leiden(anndata_obj)
	###
	FN3 = anndata_obj_name + '_cluster.png'
	sc.pl.embedding(
    	anndata_obj,
    	basis="X_mde_scanvi",
    	color="leiden",
    	frameon=False,
    	save=FN3
	)
	### Next output the meta data and the corrected matrix ###
	###
	corrected_data_lvae = lvae.get_normalized_expression()
	###
	### output corrected matrix #####
	### mtx,gene and cells ##########
	### type(corrected_data_lvae) ###
	corrected_data_lvae_array = corrected_data_lvae.values
	###
	import numpy as np
	from scipy.sparse import csr_matrix
	from scipy.io import mmwrite
	matrix = csr_matrix(corrected_data_lvae_array)
	###
	FN4 = anndata_obj_name + '_backmergeCorrect.mtx'
	mmwrite(FN4, matrix)
	###
	FN5 = anndata_obj_name + '_backmergeCorrect_Cells.csv'
	###
	### Cells = corrected_data_lvae.index.values.tolist()
	###
	Cells_meta = anndata_obj.obs
	###
	#Cells2 = Cells_meta.index.values.tolist()
	###
	#Cells == Cells2
	####
	Cells_meta.to_csv(FN5,index=False)
	####
	FN6 = anndata_obj_name + '_backmergeCorrect_Genes.csv'
	####
	Gene = anndata_obj.var['gene_name'].to_frame()
	####
	Gene.to_csv(FN6)
	####
	#### Next save tbe anndatasets ##########
	####
	import pickle
	FN7 = anndata_obj_name + '_backmergeCorrect.pkl'
	with open(FN7, "wb") as f:
    	pickle.dump(anndata_obj, f)
	#####
	print('Done!')





