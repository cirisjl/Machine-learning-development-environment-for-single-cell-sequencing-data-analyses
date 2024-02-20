import os
import subprocess
from livelogs.logger import liveLogger
from tools.formating.formating import *
from config.celery_utils import get_input_path, get_output


def run_integration(task_id, datasets, inputs,userID,output, methods, species, default_assay='RNA', output_format='Seurat', genes=None, reference=12, show_error=True):
    if methods is None:
        liveLogger.push("No integration method is selected.")
        return None
    # output = get_output_path(datasets, input, method='integration')
    # methods = [x.upper() for x in methods if isinstance(x,str)]
    # adata, counts, csv_path = LoadAnndata_to_csv(input, output, layer, show_error)

    # methods = list_py_to_r(methods)
    abs_inputList = []

    if inputs is not None:
        for input in inputs:
            if input is not None:
                abs_inputList.append(get_input_path(input, userID))


    if datasets is not None:
        dataset = datasets[0]
    datasets = list_to_string(datasets)
    methods = list_to_string(methods)
    input = list_to_string_default(abs_inputList)
    

    # #Get the absolute path for the given input
    # input = get_input_path(input, userID)
    #Get the absolute path for the given output
    output = get_output(output, userID, task_id)

    output = get_output_path(dataset, output, method='Integration', format='Seurat')
    liveLogger.push("Calling RMD1")


    try:
        report_path = get_report_path(dataset, output, "integration")
        liveLogger.push("Calling RMD2")
        # Get the absolute path of the current file
        current_file = os.path.abspath(__file__)
        # Construct the relative path to the desired file
        relative_path = os.path.join(os.path.dirname(current_file), 'integration', 'integration.Rmd')
        # Get the absolute path of the desired file
        rmd_path = os.path.abspath(relative_path)
        liveLogger.push("Calling RMD")
        liveLogger.push(datasets)
        liveLogger.push(methods)
        liveLogger.push(input)
        liveLogger.push("R -e \"rmarkdown::render('" + rmd_path + "', params=list(datasets='" + str(datasets) + "', inputs='" + str(input) + "', output_folder='" + output + "', output_format='" + output_format + "', methods='" + str(methods) + "', default_assay='" + default_assay + "', reference='" + str(reference) + "', show_error='" + str(show_error) + "'), output_file='" + report_path + "')\"")
        liveLogger.push("hello")
        s = subprocess.call(["R -e \"rmarkdown::render('" + rmd_path + "', params=list(datasets='" + str(datasets) + "', inputs='" + str(input) + "', output_folder='" + output + "', output_format='" + output_format + "', methods='" + str(methods) + "', default_assay='" + default_assay + "', reference='" + str(reference) + "', show_error='" + str(show_error) + "'), output_file='" + report_path + "')\""], shell = True)
        liveLogger.push(s)
    except Exception as e:
        liveLogger.push("Integration  is failed")
        if show_error: liveLogger.push(e)

    return {'status': 'Success'}
