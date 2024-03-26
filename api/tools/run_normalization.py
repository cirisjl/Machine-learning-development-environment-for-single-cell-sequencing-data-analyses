import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output
from utils.redislogger import *


def run_normalization(task_id, dataset, input, userID, output, methods, species, default_assay='RNA', output_format='AnnData',idtype='ENSEMBL', show_umap = True, show_error = True):
    
    if methods is None:
        redislogger.warning(task_id, "No normalization method is selected.")
        return None
    
    #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    output = get_output(output, userID, task_id)

    # methods = [x.upper() for x in methods if isinstance(x,str)]
    output = get_output_path(dataset, output, method='normalization', format='Seurat')
    # methods = list_py_to_r(methods)
    methods = list_to_string(methods)

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
        redislogger.info(task_id, s)
    except Exception as e:
        redislogger.error(task_id, "Normalization is failed.")
        if show_error: redislogger.error(task_id, f"Normalization is failed: {e}")

    return {'status': 'Success'}