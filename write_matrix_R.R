
library(Seurat)
setwd("/zp1/data/plyu3/CellAnn_final_check")
check_list <- read.table("reference_table2_final_path.txt",sep="\t",header=T)

check_process_Step3 <- function(check_list,list_number=1:439){
	######
	######
	######
	for(i in list_number){
		print(i)
		tmp_folder = check_list$Folder[i]
		print(tmp_folder)
		######
		setwd(tmp_folder)
		files = list.files()
		index = grep("_seurat_obj_new_clean$",files)
		######
		files_need = files[index]
		files_need_tag = gsub("_seurat_obj_new_clean","",files_need)
		######
		tmp_seurat = readRDS(files_need)
		###### count matrix is the count matrix !!!! ###########
		library(Seurat)
		write.csv(Matrix::t(tmp_seurat[['RNA']]@counts), file = "py_expression_matrix.csv", row.names = TRUE)
		write.csv(tmp_seurat@meta.data, file = "py_cell_metadata.csv", row.names = TRUE)
		write.csv(rownames(tmp_seurat), file = "py_gene_metadata.csv")
		###### count matrix ####################################
		######
		#####
		###### new seurat ###
		#####
	}
	print('Done!!!')
}


check_process_Step3(check_list,list_number=1:439)

print("Done!")

