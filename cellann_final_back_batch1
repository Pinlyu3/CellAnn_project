import pandas as pd
import os

file_list_path = "/zp1/data/plyu3/CellAnn_final_check/check_list3_update1.xlsx"
df = pd.read_excel(file_list_path)

df_cl = df[df['Background'] == 'M_Tabula_Muris']

df_cl.columns

file_need_path_all = df_cl['Folder'].tolist()

####### next process each folder #########

import sys
custom_module_path = '/zp1/data/plyu3/CellAnn_final_check/module'
sys.path.append(custom_module_path)

import cellann_final
import os

for file_need_path in file_need_path_all:
	print(file_need_path)
	#### change the dir now !!!! ####
	os.chdir(file_need_path)
	anndata_obj_name = cellann_final.load_anndata_from_folder_tag(file_need_path)
	anndata_obj = cellann_final.load_anndata_from_folder(file_need_path)
	####
	cellann_final.Process_anndatasets(anndata_obj,anndata_obj_name)


#######
print("All Done!!!")
