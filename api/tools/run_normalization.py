import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *


def run_normalization(dataset, input, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL', show_umap = True, show_error = True):
    
    if methods is None:
        print("No normalization method is selected.")
        return None
    output = get_output_path(dataset, output, method='normalization')
    methods = list_py_to_r(methods)

    try:
        report_path = get_report_path(dataset, output, "normalization")

        # Get the absolute path of the current file
        current_file = os.path.abspath(__file__)

        # Construct the relative path to the desired file
        relative_path = os.path.join(os.path.dirname(current_file), 'normalization', 'normalization.Rmd')

        # Get the absolute path of the desired file
        rmd_path = os.path.abspath(relative_path)

        # rmd_path = os.path.abspath("normalization/normalization.Rmd")
        s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', output_format='" + output_format + "', methods='" + methods + "', default_assay='" + default_assay + "', species=" + str(species) + "', idtype=" + str(idtype) + "', show_umap=" + str(show_umap) + ", show_error=" + str(show_error) + "), output_file='" + report_path + "')\""], shell = True)
        print(s)
    except Exception as e:
        print("Normalization is failed")
        if show_error: print(e)

    return {'status': 'Success'}