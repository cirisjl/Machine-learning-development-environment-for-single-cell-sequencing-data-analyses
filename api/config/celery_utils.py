from celery import current_app as current_celery_app
from celery.result import AsyncResult

from .celery_config import settings
from constants.declarations import USER_STORAGE
import os


def create_celery():
    celery_app = current_celery_app
    celery_app.config_from_object(settings, namespace='CELERY')
    celery_app.conf.update(task_track_started=True)
    celery_app.conf.update(task_serializer='json')
    celery_app.conf.update(result_serializer='json')
    celery_app.conf.update(accept_content=['json'])
    celery_app.conf.update(result_persistent=True)
    celery_app.conf.update(worker_send_task_events=False)
    celery_app.conf.update(worker_prefetch_multiplier=1)

    return celery_app


# def get_task_info(task_id):
#     """
#     return task info for the given task_id
#     """
#     task_result = AsyncResult(task_id)
#     summary = "Processing"
#     if task_result.ready():
#         if task_result.result is not None:
#             summary = task_result.result
#         else:
#             summary = task_result.traceback
#     result = {
#         "taskId": task_id,
#         "task_status": task_result.status,
#         "task_result": summary
#     }
#     print(result)
#     return result

#     from celery.result import AsyncResult

def get_task_info(task_id):
    """
    Return task information for the given task_id
    """
    task_result = AsyncResult(task_id)
    if task_result.ready():
        if task_result.successful():
            # Task completed successfully
            summary = task_result.get()  # This retrieves the result returned by the task
        else:
            # Task failed
            if task_result.failed():
                summary = str(task_result.result)  # Getting the exception raised in the task
                traceback = task_result.traceback  # To get the traceback if you need detailed debug info
                summary += "\n" + traceback
    else:
        summary = "Processing"

    result = {
        "taskId": task_id,
        "task_status": task_result.status,
        "task_result": summary
    }
    return result



def get_input_path(input, userID):
    """
    return the absolute input path for a given input
    """
    input_path = None
    if input is not None and userID is not None:
        input_path = USER_STORAGE + userID + input
    return input_path


def get_output(output, userID, task_id):
    """
    return the absolute input path for a given input
    """
    output_path = None
    if output is not None and userID is not None:
        output_path = output + "/" + task_id
    
    return output_path
    

def benchmarks_output_path(input):
    """
    return the absolute output path for a given input
    """
    output_path = None
    if input is not None:
        # Check if the input path is a directory
        if os.path.isdir(input):
            # If it's a directory, append '/QC' to it
            output_path = os.path.join(input, 'QC')
        else:
            # If it's a file, get the parent directory and append '/QC'
            parent_dir = os.path.dirname(input)
            output_path = os.path.join(parent_dir, 'QC')
    
    return output_path
