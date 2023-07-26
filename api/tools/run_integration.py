import os
import subprocess
from tools.formating.formating import *


def run_integration(datasets, inputs, methods, default_assay='RNA', genes=None, reference=12, show_error=True):
    if methods is None:
        print("No integration method is selected.")
        return None
    output = get_output_path(datasets, input, method='integration')
    methods = [x.upper() for x in methods if isinstance(x,str)]
    adata, counts, csv_path = load_anndata_to_csv(input, output, layer, show_error)

    methods = list_py_to_r(methods)

    try:
        report_path = get_report_path(datasets, output, "normalization")
        rmd_path = os.path.abspath("integration/integration.Rmd")
        s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(dataset='" + str(datasets) + "', input='" + input + "', output='" + output + "', output_format='" + output_format + "', methods='" + methods + "', default_assay='" + default_assay + "', species=" + str(species) + "', idtype=" + str(idtype) + "', show_umap=" + str(show_umap) + ", show_error=" + str(show_error) + "), output_file='" + report_path + "')\""], shell = True)
        print(s)
    except Exception as e:
        print("Normalization is failed")
        if show_error: print(e)

    return {'status': 'Success'}
