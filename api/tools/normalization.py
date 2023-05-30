import os
import subprocess
# import sys
# sys.path.append('..')
from tools.formating.formating import *


def normalize(dataset, input, output, methods, default_assay='RNA', output_format='AnnData', species=None, idtype='ENSEMBL', show_umap = True, show_error = True):
    
    if methods is None:
        print("No normalization method is selected.")
        return None
    output = output_path_check(dataset, output)
    methods = list_py_to_r(methods)

    try:
        report_path = get_report_path(dataset, output, "normalization")
        rmd_path = os.path.abspath("normalization.Rmd")
        s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(dataset='" + str(dataset) + "', input='" + input + "', output='" + output + "', output_format='" + output_format + "', methods='" + methods + "', default_assay='" + default_assay + "', species=" + str(species) + "', idtype=" + str(idtype) + "', show_umap=" + str(show_umap) + ", show_error=" + str(show_error) + "), output_file='" + report_path + "')\""], shell = True)
        print(s)
    except Exception as e:
        print("Normalization is failed")
        if show_error: print(e)

    return {'status': 'Success'}