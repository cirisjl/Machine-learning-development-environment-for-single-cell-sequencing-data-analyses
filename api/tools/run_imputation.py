import os
import subprocess
import sys
from tools.formating.formating import *
from tools.imputation.MAGIC import magic_impute
    

def run_imputation(dataset, input, output, methods, layer=None, genes=None, ncores=12, show_error=True):
    if methods is None:
        print("No imputation method is selected.")
        return None
    methods = [x.upper() for x in methods if isinstance(x,str)]
    adata, counts, csv_path = load_anndata_to_csv(input, output, layer, show_error)
         
    if adata is None:
        print("File format is not supported.")
        return None 
    
    if "MAGIC" in methods:
        if 'MAGIC_imputed' not in adata.layers.keys(): 
            try:
                data_magic = magic_impute(counts, genes)
                adata.layers['MAGIC_imputed'] = data_magic
                output = get_output_path(dataset, input, method='MAGIC_imputation')
                adata.write_h5ad(output, compression='gzip')
                print("AnnData object for MAGIC imputation is saved successfully")
            except Exception as e:
                print("MAGIC imputation is failed")
                if show_error: print(e)
        else: 
            print("'MAGIC_imputed' layer already exists.")

    if "scGNN" in methods:
        if 'scGNN_imputed' not in adata.layers.keys(): 
            try:
                output = get_output_path(dataset, input, method='scGNN_imputation')
                print("AnnData object for scGNN imputation is saved successfully")          
            except Exception as e:
                print("scGNN imputation is failed")
                if show_error: print(e)
        else: 
            print("'scGNN_imputed' layer already exists.") 
    
    if "SAVER" in methods:
        if 'SAVER_imputed' not in adata.layers.keys(): 
            try:
                output = get_output_path(dataset, input, method='SAVER_imputation')
                report_path = get_report_path(dataset, output, "SAVER")
                saver_path = os.path.abspath("imputation/SAVER.Rmd")
                s = subprocess.call(["R -e \"rmarkdown::render('" + saver_path + "', params=list(dataset='" + str(dataset) + "', input='" + csv_path + "', output='" + output + "', output_format='AnnData', ncores=" + str(ncores) + "), output_file='" + report_path + "')\""], shell = True)
                print(s)
            except Exception as e:
                print("SAVER imputation is failed")
                if show_error: print(e)
        else: 
            print("'SAVER_imputed' layer already exists.")

    return {'status': 'Success'}

    
    

        
            

