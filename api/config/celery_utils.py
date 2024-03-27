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


def get_task_info(task_id):
    """
    return task info for the given task_id
    """
    task_result = AsyncResult(task_id)
    result = {
        "task_id": task_id,
        "task_status": task_result.status,
        "task_result": task_result.result
    }
    return result


def get_input_path(input, userID):
    """
    return the absolute input path for a given input
    """
    if input is not None and userID is not None:
        input_path = USER_STORAGE + userID + input
        return input_path


def get_output(output, userID, task_id):
    """
    return the absolute input path for a given input
    """
    if output is not None and userID is not None:
        output_path = output + "/" + task_id
        return output_path

def benchmarks_output_path(input, task_id):
    """
    return the absolute output path for a given input
    """
    if input is not None and task_id is not None:
        # Check if the input path is a directory
        if os.path.isdir(input_path):
            # If it's a directory, append '/Results' to it
            output_path = os.path.join(input_path, 'Results', task_id)
        else:
            # If it's a file, get the parent directory and append '/Results'
            parent_dir = os.path.dirname(input_path)
            output_path = os.path.join(parent_dir, 'Results', task_id)
    
        return output_path
