##### #####
##### #####


##### we will perform the Final Figure5: ###########
#####

##### first we will see the clusters of the reference datasets: #######
#####
ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R


pmid31075224_eye_test_CellAnn_Step1_input.txt

setwd("/zp1/data/plyu3/UCSC_old/pmid31075224")

library(Seurat)

seurat_obj <- readRDS("pmid31075224_seurat_obj_new")
seurat_obj <- FindNeighbors(seurat_obj,reduction='dims',dims=1:2)
seurat_obj <- FindClusters(seurat_obj,reduction=0.3)

devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/CellAnn/main/prepare_CellAnn.R")

prepare_CellAnn(seurat_obj,folder="/zp1/data/plyu3/UCSC_old/pmid31075224",sample_name='pmid31075224_eye_test',matrix_name='RNA',dims='dims',cluster='seurat_clusters')

library(ggplot2)

######
mkdir /zp1/data/plyu3/CellAnn_figure5

#######
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure51A.png")
png(png_file,height=4000,width=6000,res=72*12)
print(DimPlot(seurat_obj, reduction = "dims",group.by=c('seurat_clusters'),label = TRUE, label.size = 2.5, repel = TRUE,raster = FALSE))
dev.off()

#######
Multi-species single-cell transcriptomic analysis of ocular compartment regulons

34584087_scRNA_Human_Ocular_Retina

folder = "/home/plyu3/cell_ann_part_datasets1/34584087_scRNA_Human Ocular_Retina"

files = list.files()

setwd(folder)
seurat_obj_2 = readRDS("pmid_34584087_scRNA_Human_Ocular_Retina_seurat_obj_new")

library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure52_1.png")
png(png_file,height=4000,width=6000,res=72*12)
print(DimPlot(seurat_obj_2, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE))
dev.off()

#######
####### sample sample #######
#######


#######
####### we will set the UMAPs of this study ########
#######
Integration of eQTL and a Single-Cell Atlas in the Human Eye Identifies Causal Genes for Age-Related Macular Degeneration
31995762_total

folder = "/home/plyu3/cell_ann_part_datasets1/31995762_total"
setwd(folder)

seurat_obj_2 = readRDS("31995762_total_seurat_obj_new")

library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure52_2.png")
png(png_file,height=4000,width=6000,res=72*12)
print(DimPlot(seurat_obj_2, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE))
dev.off()



#######
Construction of a human cell landscape at single-cell level

32214235_Human_cell_landscape_eye

setwd("/zp1/data/plyu3/Altas_add_CellAnn/Human_cell_landscape/eye_Human_cell_landscape")

seurat_obj_2 = readRDS("eye_seurat_obj_new")

library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure52_3.png")
png(png_file,height=4000,width=8000,res=72*12)
print(DimPlot(seurat_obj_2, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE))
dev.off()

#######
####### OK!!! Next #######
#######

folder = "/zp1/data/plyu3/UCSC_old/pmid31075224"
TAG = "31075224_Human_eye"

setwd(folder)

seurat_obj_2 = readRDS("pmid31075224_seurat_obj_new")


library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure53_1.png")
png(png_file,height=4000,width=7000,res=72*12)
print(DimPlot(seurat_obj_2, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE,cols = cols))
dev.off()

col_list = c("#EF9000","#F6BA00","#9EA220","#86BF38","#EC6F64","#D11536","#936DAD","#A74997","#61BFB9","#0097BB","#026AB1")

cols = c("amacrine cell" = "#EF9000","bipolar neuron" = "#936DAD","endothelial cell" = "#61BFB9","foveal cone photoreceptor" = "#9EA220","glial cell"="#D11536","microglial cell"="#EC6F64","pericyte cell"="#A74997","peripheral cone photoreceptor"="#9EA220","retina horizontal cell"="grey","retinal ganglion cell"="#026AB1","retinal rod cell"="#86BF38","Unknown"="black")

"amacrine cell"
"bipolar neuron"
"endothelial cell"
"foveal cone photoreceptor"
"glial cell"
"microglial cell"
"pericyte cell"
"peripheral cone photoreceptor"
"retina horizontal cell"
"retinal ganglion cell"
"retinal rod cell"
"Unknown"



alignment_results <- read.table("/zp1/data/plyu3/CellAnn_figure5/Alignment_results_T8Fu0y4Jqv.txt",sep="\t",header=T)

alignment_results$Custom.label[9] <- "Bipolar"
alignment_results$Custom.label[10] <- "MG"
alignment_results$Custom.label[11] <- "Rod"
alignment_results$Custom.label[12] <- "Rod"
alignment_results$Custom.label[13] <- "Bipolar"
alignment_results$Custom.label[14] <- "Bipolar"
alignment_results$Custom.label[15] <- "MG"
alignment_results$Custom.label[16] <- "Bipolar"
alignment_results$Custom.label[17] <- "Bipolar"
alignment_results$Custom.label[18] <- "MG"
alignment_results$Custom.label[19] <- "MG"
alignment_results$Custom.label[20] <- "RGC"
alignment_results$Custom.label[21] <- "MG"


alignment_results$Custom.label[22] <- "MG"
alignment_results$Custom.label[23] <- "RGC"
alignment_results$Custom.label[24] <- "Cone"
alignment_results$Custom.label[25] <- "HC"
alignment_results$Custom.label[26] <- "endothelial"
alignment_results$Custom.label[27] <- "Bipolar"
alignment_results$Custom.label[28] <- "Bipolar"
alignment_results$Custom.label[29] <- "Rod"
alignment_results$Custom.label[30] <- "Pericytes"
alignment_results$Custom.label[31] <- "MG"
alignment_results$Custom.label[32] <- "endothelial"
alignment_results$Custom.label[33] <- "Microglia"
alignment_results$Custom.label[34] <- "AC"


alignment_results$Custom.label[35] <- "RGC"
alignment_results$Custom.label[36] <- "Rod"
alignment_results$Custom.label[37] <- "Rod"
alignment_results$Custom.label[38] <- "Cone"
alignment_results$Custom.label[39] <- "Cone"
alignment_results$Custom.label[40] <- "Microglia"
alignment_results$Custom.label[41] <- "Astrocyte"
alignment_results$Custom.label[42] <- "Rod"
alignment_results$Custom.label[43] <- "MG"

table(alignment_results$Custom.label)

input4 <- read.table("pmid31075224_eye_test_CellAnn_Step4_input.txt",sep='\t',header=T)

m = match(input4$cluster,alignment_results$cluster)
input4$celltype = alignment_results$Custom.label[m]

m1 = match(colnames(seurat_obj_3),input4$cell)
seurat_obj_3$celltype = input4$celltype[m1]

library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure53_1.png")
png(png_file,height=4000,width=7000,res=72*12)
print(DimPlot(seurat_obj_3, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE,cols = cols))
dev.off()


cols = c("AC" = "#EF9000","Bipolar" = "#936DAD","endothelial" = "#61BFB9","Cone" = "#9EA220","MG"="#D11536","Microglia"="#EC6F64","Pericytes"="#A74997","HC"="grey","RGC"="#026AB1","Rod"="#86BF38","Astrocyte"="blue")


library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure53_2.png")
png(png_file,height=4000,width=6000,res=72*12)
print(DimPlot(seurat_obj_3, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE,cols = cols))
dev.off()



### then we add colors ######
### 






####### Mouse datasets ########
#######
Gene regulatory networks controlling vertebrate retinal regeneration

33004674_Mouse_injury_retina

folder = "/zp1/data/plyu3/Eye_add_CellAnn/33004674_Mouse_injury_retina"
setwd(folder)


seurat_obj_2 = readRDS("33004674_Mouse_injury_retina_seurat_obj_new")

library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_figure5")
png_file = paste0("Figure52_4.png")
png(png_file,height=4000,width=6000,res=72*12)
print(DimPlot(seurat_obj_2, reduction = "dims",group.by=c('celltype'),label = FALSE, label.size = 2.5, repel = TRUE,raster = FALSE))
dev.off()


#######
#######
#######
#######
#######








































