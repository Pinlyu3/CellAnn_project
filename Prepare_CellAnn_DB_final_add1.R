ssh plyu3@omb2.onc.jhmi.edu
U[9C20&&

conda activate seurat4
R
####

source("/zp1/data/plyu3/UCSC/Prepare_CellAnn_DB_final_source.R")
source("/zp1/data/plyu3/UCSC/Prepare_CellAnn_DB_final_add1.R")

#####
#####
#####


Read_Seurat_file_in_the_folder_UCSC <- function(input_folder){
	##### first we init a table ######
	checked_items = c('Seurat_files_exists:','Seurat_file_celltypes:','Seurat_file_dims:','Seurat_file_matrix_equal:','Seurat_file_matrix_max:','Seurat_file_matrix_log:','Seurat_file_matrix_log2:','Seurat_file_matrix_counts:','Seurat_file_Genes:','Seurat_file_matrix_counts_sum:')
	##### 
	checked_table = data.frame(checked_items=checked_items,results="Need to check!")
	##### loading library !!! #####
	library(Seurat)
	##### we first go into the folder ######
	setwd(input_folder)
	##### 'Seurat_files_exists:' #####
	files = list.files()
	#####
	##### searching seurat files ########
	##### 
	k = grep("_seurat_obj$",files)
	#####
	if(length(k) != 1){
		checked_table$results[1] = 'No file!'
		write.table(checked_items,file='checked_items.txt',sep="\t",quote=F,row.names=F)
		return(checked_items)
	}
	if(length(k) == 1){
		checked_table$results[1] = 'OK!'
		###### Next we will load the seurat objects !!!! ##########
		seurat_file = files[k]
		seurat_obj <- readRDS(seurat_file)
		###### Next we will check whether celltype exists: ########
		k1 = which(colnames(seurat_obj@meta.data) == "celltype")
		if(length(k1) != 1){
			checked_items$results[2] = 'No celltype!'
			write.table(checked_table,file='checked_items.txt',sep="\t",quote=F,row.names=F)
			return(checked_items)
		}
		if(length(k) == 1){
			checked_table$results[2] = 'OK!'
			###### Next we will check the dims !!!! ##############
			k2 = which(colnames(seurat_obj@meta.data) %in% c("dim1","dim2","dim3") == T)
			if(length(k2) < 2){
				checked_table$results[3] = 'No dims!'
				#write.table(checked_table,file='checked_items.txt',sep="\t",quote=F,row.names=F)
				#return(checked_items)
			}
			if(length(k2) > 1){
				checked_table$results[3] = 'OK!'
			}
			##### Next check whether matrix is equal!!! #######
			##### RNA@data matrix and RNA@counts matrix #######
			data_mat = seurat_obj[['RNA']]@data[1:100,1:100]
			count_mat = seurat_obj[['RNA']]@counts[1:100,1:100]
			#####
			res = all.equal(data_mat,count_mat)
			#####
			checked_table$results[4] = res
			#####
			k_NAN = which(is.nan(colSums(seurat_obj)) == T)
			#####
			if(length(k_NAN) > 0){
				seurat_obj = seurat_obj[,-k_NAN]
			}
			#####
			count_mat = seurat_obj[['RNA']]@counts[,1:100]
			max = max(as.vector(count_mat))
			#####
			if(max < 20){
				checked_table$results[5] = max
			}
			if(max > 50){
				checked_table$results[5] = max
			}
				####### top5 sums #####
				####
				sums = as.character(round(Matrix::colSums(exp(count_mat))[1:5],2))
				####
				sums_res = paste(sums,collapse='::')
				####
				checked_table$results[6] = sums_res
				####
				sums = as.character(round(Matrix::colSums(2^(count_mat))[1:5],2))
				sums_res = paste(sums,collapse='::')
				checked_table$results[7] = sums_res
				####
				sums = as.character(round(Matrix::colSums(count_mat)[1:5],2))
				sums_res = paste(sums,collapse='::')
				checked_table$results[10] = sums_res
				#### print the non-Zero number in the matrix #######
				####
				count_mat_v <- as.vector(count_mat)
				count_mat_v = count_mat_v[which(count_mat_v > 0)]
				count_mat_v = round(count_mat_v,3)
				sums_res = paste(count_mat_v[1:5],collapse='::')
				checked_table$results[8] = sums_res
				#### 
				######## Next we will check the genes !!! ##########
				####
				#### plot the first rownames of the seurat objects !!! #####
				####
				Genes = rownames(seurat_obj)[1:5]
				####
				sums_res = paste(Genes[1:5],collapse='::')
				checked_table$results[9] = sums_res
				####
				write.table(checked_table,file='checked_items.txt',sep="\t",quote=F,row.names=F)
				return(checked_table)
		}
	}
	#####
	#####
}



Prepare_Seurat_file_in_the_folder_UCSC <- function(input_folder,norm="log"){
	#########
	library(Seurat)
	##### we first go into the folder ######
	setwd(input_folder)
	##### 'Seurat_files_exists:' #####
	files = list.files()
	#####
	##### searching seurat files ########
	##### 
	k = grep("_seurat_obj$",files)
	#####
	#########
	seurat_file = files[k]
	seurat_obj <- readRDS(seurat_file)
	print(dim(seurat_obj))
	#########
	seurat_obj_colsum = Matrix::colSums(seurat_obj[['RNA']]@counts)
	k_NAN = which(is.nan(seurat_obj_colsum) == T)
	if(length(k_NAN) > 0){
		seurat_obj = seurat_obj[,-seurat_obj_colsum]
	}
	#########
	if(dim(seurat_obj)[2] > 250000){
		index = sample(1:dim(seurat_obj)[2],250000)
		seurat_obj = seurat_obj[,index]
	}
	##### first we convert matrix to counts ##########
	if(norm == "log"){
		######
		mat = seurat_obj[['RNA']]@counts
		mat_new = exp(mat)-1
		######
	}
	if(norm == "counts"){
		######
		mat = seurat_obj[['RNA']]@counts
		mat_new = mat
		######
	}
	if(norm == "log2"){
		######
		mat = seurat_obj[['RNA']]@counts
		mat_new = 2^(mat)-1
		######
	}
	#####
	k_NA_count = which(is.na(rowSums(mat_new)) == F)
	if(length(k_NA_count) > 0){
		mat_new = mat_new[k_NA_count,]
	}
	#####
	seurat_obj_new = CreateSeuratObject(mat_new)
	seurat_obj_new@meta.data = seurat_obj@meta.data
	#####
	##### next we add dims to this seurat_obj #######
	#####
	dims_index = c("dim1","dim2","dim3")
	k = match(dims_index,colnames(seurat_obj_new@meta.data))
	k = k[is.na(k) == F]
	if(length(k) > 0){
		k = k[is.na(k) == F]
		UMAPs_mat = as.matrix(seurat_obj@meta.data[,k])
		UMAPs_mat_new = Norm_UMAPs(UMAPs_mat)
	}
	if(length(k) == 0){
		print("add custom UMAPs !!!!")
		seurat_obj_new <- NormalizeData(seurat_obj_new)
		seurat_obj_new <- FindVariableFeatures(seurat_obj_new)
		seurat_obj_new <- ScaleData(seurat_obj_new)
		seurat_obj_new <- RunPCA(seurat_obj_new)
		seurat_obj_new <- RunUMAP(seurat_obj_new,dims=1:30)
		#######
		UMAPs_mat = as.matrix(seurat_obj_new[["umap"]]@cell.embeddings)
		UMAPs_mat_new = Norm_UMAPs(UMAPs_mat)
	}
	#####
	#####
	#####
	seurat_obj_new[['dims']] <- CreateDimReducObject(embeddings = UMAPs_mat_new, key = "dims_", assay = 'RNA')
	#####
	#png_file = paste0("dim_check.png")
	#png(png_file,height=5000,width=6000,res=72*12)
	#print(DimPlot(seurat_obj_new, reduction = "dims",group.by=c('celltype'),label = TRUE, label.size = 2.5, repel = TRUE,raster = FALSE))
	#dev.off()
	###### 
	###### Next we will create the files #######
	######
	seurat_file_new = paste(seurat_file,'_new',sep="")
	######
	saveRDS(seurat_obj_new,file=seurat_file_new)
}



#### OK!!! ######
#### then we write small functions ########
####

#### 这个是OK的 ！！！##### 
####
Norm_UMAPs <- function(UMAPs_mat){
	#### each rows and cols norms to the -50 50 dims #######
	####
	UMAPs_mat_new = UMAPs_mat
	####
	umap_dims = dim(UMAPs_mat)[2]
	####
	for(i in 1:dim(UMAPs_mat)[2]){
		tmp_vector = UMAPs_mat[,i]
		tmp_vector_mean = tmp_vector - mean(tmp_vector)
		tmp_vector_scale = 100 / (max(tmp_vector) - min(tmp_vector))
		tmp_vector_norm = tmp_vector_mean*tmp_vector_scale
		UMAPs_mat_new[,i] = round(tmp_vector_norm,3)
	}
	####
	return(UMAPs_mat_new)
}


runDEGs_Ref <- function(Seurat_Obj,method='COSG',idents='celltype',num_of_genes = 100){
	Idents(Seurat_Obj) = idents
	#######
	library(COSG)
	#######
	if(method == 'COSG'){
		marker_cosg <- cosg(
 				Seurat_Obj,
 				groups=c('all'),
 				assay='RNA',
 				slot='data',
 				mu=1,
 				n_genes_user=num_of_genes)
		res = marker_cosg$names
	}
	#######
	#######
	return(res)
}

#####
#####
##### input is the log transformed matrix #######
#####
#####

CellAnn_Avg_Mat <- function(data_mat,data_cluster,log='log',scale_factor=10000){
	######
	tag_cluster = unname(data_cluster)
	tag_cluster_level = levels(as.factor(tag_cluster))
	###### normalized back datasets ######
	if(log == 'log'){
		data_mat_exp = exp(data_mat)
		data_mat_exp = data_mat_exp-1
	}
	if(log == 'log2'){
		data_mat_exp = 2^(data_mat)
		data_mat_exp = data_mat_exp-1
	}
	print(paste('Sums:',head(colSums(data_mat_exp[,c(1:2)]))))
	###### data_mat_exp is 1e5 normalize #######
	merge_mat = c()
	for(i in 1:length(tag_cluster_level)){
		index = which(data_cluster %in% tag_cluster_level[i] == T)
		index_mat = data_mat_exp[,index]
		######
		index_sum = rowSums(index_mat)
		######
		merge_mat = c(merge_mat,index_sum)
	}
	###
	merge_mat = matrix(merge_mat,nrow=dim(data_mat)[1])
	###
	rownames(merge_mat) = rownames(data_mat)
	colnames(merge_mat) = tag_cluster_level
	### colSums(merge_mat)
	scale = colSums(merge_mat)/scale_factor
	merge_mat = sweep(merge_mat,2,scale,FUN='/')
	### default norm ####
	merge_mat = round(log(merge_mat+1),3)
	return(merge_mat)
	### return is a log transformed matrix ####
}

##############
############## this function subset each celltype annotation ######
############## keys = "celltype"
############## 变成 subset seurat ##################################
############## 

subset_each_ct <- function(data,resolution){
	#####
	matrix_list = list()
	##### ct ######
	ct = levels(as.factor(data$celltype))
	#####
	cells_list = list()
	#####
	for(j in 1:length(ct)){
		print(ct[j])
		print(paste("j = ",j))
		#######
		k = which(data$celltype == ct[j])
		#######
		if(length(k) > 5){
			print('large')
			#######
			sub_seurat = subset(data,subset = celltype == ct[j])
			#######
			dim_length = dim(sub_seurat@reductions$dims@cell.embeddings)[2]
			#######
			sub_seurat = FindNeighbors(sub_seurat,reduction="dims",dims=c(1:dim_length))
			sub_seurat = FindClusters(sub_seurat,resolution=resolution)
			sub_seurat_mat = sub_seurat[['RNA']]@data
			data_cluster = sub_seurat$seurat_clusters
			data_cluster = paste(ct[j],data_cluster,sep='@sub')
			sub_seurat_avg = CellAnn_Avg_Mat(data_mat=sub_seurat_mat,data_cluster)
			matrix_list = c(matrix_list,list(sub_seurat_avg))
			######
			cells_table = data.frame(cells=colnames(sub_seurat_mat),new_sub_cluster = data_cluster)
			cells_list = c(cells_list,list(cells_table))
		}
	}
	matrix_list = do.call('cbind',matrix_list)
	cells_list = do.call('rbind',cells_list)
	####
	return(list(matrix_list,cells_list))
}


######
######

Prepare_Seurat_file_in_the_folder_UCSC_STEP2 <- function(input_folder){
	#########
	library(Seurat)
	##### we first go into the folder ######
	setwd(input_folder)
	##### 'Seurat_files_exists:' #####
	files = list.files()
	#####
	##### searching seurat files ########
	##### 
	k = grep("_seurat_obj_new$",files)
	#####
	#########
	seurat_file = files[k]
	seurat_obj <- readRDS(seurat_file)
	######### Performing log-normalization ###########
	##### first we convert matrix to counts ##########
	#seurat_obj[['RNA']]@counts = seurat_obj[['RNA']]@counts*100
	seurat_obj <- NormalizeData(seurat_obj)
	#####
	#####
	subset_matrix_res = subset_each_ct(seurat_obj,resolution=0.5)
	##### OK! what is the next ? #####################
	#####
	subset_matrix = subset_matrix_res[[1]]
	subset_cell = subset_matrix_res[[2]]
	#####
	m = match(colnames(seurat_obj),subset_cell$cells)
	seurat_obj$cell_id = subset_cell$cells[m]
	seurat_obj$celltype_sub = subset_cell$new_sub_cluster[m]
	##### rm NA cells #####
	k_NA = which(is.na(seurat_obj$celltype_sub) == T)
	if(length(k_NA) > 0){
		seurat_obj = seurat_obj[,-k_NA]
	}
	#####
	##### OK!!! we finished annotation these cells, next to check to see the results: #########
	#####
	png_file = paste0("sub_ct_check.png")
	png(png_file,height=5000,width=26000,res=72*12)
	print(DimPlot(seurat_obj, reduction = "dims",group.by=c('celltype_sub'),label = FALSE,raster = FALSE))
	dev.off()
	###### OK!! Next we will see the results !!! ############
	###### first we generate the avg expression matrix !!! ##
	######
	celltype_avg = CellAnn_Avg_Mat(seurat_obj[['RNA']]@data,seurat_obj$celltype)
	sub_celltype_avg = CellAnn_Avg_Mat(seurat_obj[['RNA']]@data,seurat_obj$celltype_sub)
	###### then we will call markers ########################
	######
	celltype_marker = runDEGs_Ref(seurat_obj,method='COSG',idents='celltype',num_of_genes=100)
	sub_celltype_marker = runDEGs_Ref(seurat_obj,method='COSG',idents='celltype_sub',num_of_genes=100)
	###### then we will save these files ####################
	######
	saveRDS(seurat_obj,file=files[k])
	######
	saveRDS(celltype_avg,file=paste0(files[k],"_avg_mat"))
	saveRDS(sub_celltype_avg,file=paste0(files[k],"_avg_sub_mat"))
	saveRDS(celltype_marker,file=paste0(files[k],"_marker"))
	saveRDS(sub_celltype_marker,file=paste0(files[k],"_sub_marker"))
	######
	###### also need to provide dimplot file !! #######
	###### 
	dim_files = data.frame(seurat_obj@reductions$dims@cell.embeddings)
	dim_files_width = dim(dim_files)[2]
	colnames(dim_files) <- paste0("dim",1:dim_files_width)
	######
	dim_files$cluster = seurat_obj$celltype
	###### re-order dim files #######
	dim_files = dim_files[,c(dim_files_width+1,1:dim_files_width)]
	print(head(dim_files))
	saveRDS(dim_files,file=paste0(files[k],"_Dimplot"))
	######
	print("Finished !!!")
}

########## should convert to “pmid” + PMID number + “_Dimplot” tags ######
##########

Prepare_Seurat_file_in_the_folder_UCSC_STEP2_copy_folder <- function(From_folder,To_folder,TAG){
	############# we need to search files in the From_folder ##########
	setwd(From_folder)
	#############
	files = list.files()
	#############
	### first: "_avg_mat" ##
	k1 = grep("_avg_mat$",files)
	k2 = grep("_avg_sub_mat$",files)
	k3 = grep("_Dimplot$",files)
	k4 = grep("new_marker$",files)
	k5 = grep("new_sub_marker$",files)
	#############
	file_avg_mat_target = paste0('pmid',TAG,'_avg_mat')
	file_avg_sub_mat_target = paste0('pmid',TAG,'_avg_sub_mat')
	file_Dimplot_target = paste0('pmid',TAG,'_Dimplot')
	file_marker_target = paste0('pmid',TAG,'_marker')
	file_sub_marker_target = paste0('pmid',TAG,'_sub_marker')
	#############
	file_avg_mat_target2 = paste0(To_folder,file_avg_mat_target)
	file_avg_sub_mat_target2 = paste0(To_folder,file_avg_sub_mat_target)
	file_Dimplot_target2 = paste0(To_folder,file_Dimplot_target)
	file_marker_target2 = paste0(To_folder,file_marker_target)
	file_sub_marker_target2 = paste0(To_folder,file_sub_marker_target)
	############# we use a command to copy these files ####################
	#############
	command1 = paste('cp',files[k1],file_avg_mat_target2)
	command2 = paste('cp',files[k2],file_avg_sub_mat_target2)
	command3 = paste('cp',files[k3],file_Dimplot_target2)
	command4 = paste('cp',files[k4],file_marker_target2)
	command5 = paste('cp',files[k5],file_sub_marker_target2)
	print(command1)
	print(command2)
	print(command3)
	print(command4)
	print(command5)
	#############
	system(command1)
	system(command2)
	system(command3)
	system(command4)
	system(command5)
	##############
	Datasheet_table = data.frame(PMID=TAG,reference_index=paste0('pmid',TAG),reference_avg_mat=file_avg_mat_target,reference_avg_sub_mat=file_avg_sub_mat_target,reference_marker=file_marker_target,reference_sub_marker=file_sub_marker_target)
	##############
	setwd(To_folder)
	files = list.files()
	##############
	k6 = grep("reference_table2_add.txt",files)
	##############
	if(length(k6) == 0){
		##### then we write it as a new table ############
		write.table(Datasheet_table,file="reference_table2_add.txt",sep="\t",row.names=F,quote=F)
	}
	if(length(k6) > 0){
		print("Find!!!")
		##### then we write it as a new table ############
		Datasheet_table_Ori = read.table(files[k6],sep="\t",header=T)
		#####
		Datasheet_table = rbind(Datasheet_table_Ori,Datasheet_table)
		#####
		Datasheet_table = Datasheet_table[!duplicated(Datasheet_table$PMID),]
		#####
		write.table(Datasheet_table,file="reference_table2_add.txt",sep="\t",row.names=F,quote=F)
	}
}

Broad_matrix_process <- function(file,sep="\t"){
	Sys.setenv("VROOM_CONNECTION_SIZE"= 131072 * 500)
	library(readr)
	matrix = read_delim(file,delim=sep)
	#### to check the format of the matrix #####
	####
	mat_input = as.matrix(matrix[,-1])
	rownames(mat_input) = data.frame(matrix[,1])[,1]
	print(mat_input[1:5,1:5])
	#######
	return(mat_input)
}

Broad_cluster_process <- function(file,sep='\t'){
	library(readr)
	#####
	cluster = read_delim(file,delim=sep)
	cluster = cluster[-1,]
	#####
	print(head(cluster))
	#####
	cluster = data.frame(cluster)
	return(cluster)
}

Broad_create_seurat_obj <- function(mat_input,cluster,file){
	######
	library(Seurat)
	###### saveRDS #######
	seurat_obj = CreateSeuratObject(mat_input)
	######
	m = match(colnames(seurat_obj),cluster[,1])
	cluster = cluster[m,]
	######
	k = m[is.na(m) == F]
	print(paste0("ratio:",length(k)/length(m)))
	######
	head(cluster)
	head(colnames(seurat_obj),n=2)
	######
	seurat_obj@meta.data = cbind(seurat_obj@meta.data,cluster)
	######
	k = which(colnames(seurat_obj@meta.data) %in% c("dim1","dim2","dim3") == T)
	if(length(k) > 0){
		for(j in k){
			seurat_obj@meta.data[,j] = as.numeric(seurat_obj@meta.data[,j])
		}
	}
	######
    k_NA = which(is.na(seurat_obj$celltype) == T)
    if(length(k_NA) > 0){
    	print(k_NA)
    	seurat_obj = seurat_obj[,-k_NA]
    }
	######
	saveRDS(seurat_obj,file=file)
}


######## Next we will upload the files to the 
######## we will first update the datasbase_sheet #########
########

######## tissue names all ########

















