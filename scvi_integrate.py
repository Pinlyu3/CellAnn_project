#!/usr/bin/python

#####################
# import the module #
#####################
import subprocess
import argparse
import os
import re
import pandas as pd
import scanpy as sc
import scvi
import scib
from scvi.model.utils import mde
import pymde

#####################
# get the parameter #
#####################

p = argparse.ArgumentParser(usage = "Get parameter", description = "run scvi piplines")

p.add_argument('--folder')
p.add_argument('--file')

args = p.parse_args()

folder = args.folder
file = args.file

print "folder = " + folder
print "Genome = " + file

#######
####### input is a folder and a txt data frame ##################
####### let us to test ######
####### 
####### folder = "/home/lp123/Desktop/CellAnn_integrate_folder/H5ad"
####### file = "/home/lp123/Desktop/CellAnn_integrate_folder/reference_prepare_table1_202211.txt"
#######

import pandas as pd

####### change the file folder #######

os.chdir(folder)

####### read the txt file ############
####### 

processing_table = pd.read_csv(file,sep='\t')

file_lists = list(processing_table['reference_background_h5ad'])

#######
####### get the output file names ########
#######

file_lists_out = []
for i in range(len(file_lists)):
    #print(file_lists[i])
    new_file = file_lists[i].replace('.h5ad','.tsv')
    #print(new_file)
    file_lists_out.append(new_file)

####### 


##### let us define the main function #####
def Mainscvi(file_input,file_output):
    global_seurat = sc.read(file_input)
    ###
    global_seurat.layers["counts"] = global_seurat.X
    ###
    sc.pp.highly_variable_genes(
        global_seurat,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="batch",
        subset=True
    )
    ####
    scvi.model.SCVI.setup_anndata(global_seurat, layer="counts", batch_key="batch")
    vae = scvi.model.SCVI(global_seurat, n_layers=2, n_latent=30, gene_likelihood="nb")
    vae.train()
    ####
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=global_seurat,
        labels_key="celltype",
        unlabeled_category="Unknown",
    )
    ####
    lvae.train(max_epochs=20, n_samples_per_label=100)
    global_seurat.obsm["X_scANVI"] = lvae.get_latent_representation(global_seurat)
    #### ####
    global_seurat_scANVI = pd.DataFrame(global_seurat.obsm['X_scANVI'])
    global_seurat_scANVI.add_prefix('SCANVI_')
    Cellnames = list(pd.DataFrame(global_seurat.obs['_scvi_batch']).index)
    global_seurat_scANVI['Cell_ID'] = Cellnames
    ####
    global_seurat_scANVI.to_csv(file_output, sep="\t",index=False)



#### OK!!!! we will output each samples ####


#### Next we will try to process each file #####
for i in range(len(file_lists)):
    print("input =", file_lists[i])
    print("output =", file_lists_out[i])
    Mainscvi(file_lists[i],file_lists_out[i])



print("All Done !!!!! ")

#### OK!!!!! #####



