import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output


def run_normalization(dataset, input, userID, output, methods, species, default_assay='RNA', output_format='AnnData',idtype='ENSEMBL', show_umap = True, show_error = True):
    
    if methods is None:
        print("No normalization method is selected.")
        return None
    
    #Get the absolute path for the given input
    input = get_input_path(input, userID)
    #Get the absolute path for the given output
    output = get_output(output, userID)

    methods = [x.upper() for x in methods if isinstance(x,str)]
    output = get_output_path(dataset, output, method='normalization', format='Seurat')
    # methods = list_py_to_r(methods)
    methods = methods_list(methods)

    try:
        report_path = get_report_path(dataset, output, "normalization")

        # Get the absolute path of the current file
        current_file = os.path.abspath(__file__)

        # Construct the relative path to the desired file
        relative_path = os.path.join(os.path.dirname(current_file), 'normalization', 'normalization.Rmd')

        # Get the absolute path of the desired file
        rmd_path = os.path.abspath(relative_path)

        # rmd_path = os.path.abspath("normalization/normalization.Rmd")
        s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', output_format='" + output_format + "', methods='" + methods + "', default_assay='" + default_assay + "', species='" + str(species) + "', idtype='" + str(idtype) + "', show_umap='" + str(show_umap) + "', show_error='" + str(show_error) + "'), output_file='" + report_path + "')\""], shell = True)
        print(s)
    except Exception as e:
        print("Normalization is failed")
        if show_error: print(e)

    return {'status': 'Success'}