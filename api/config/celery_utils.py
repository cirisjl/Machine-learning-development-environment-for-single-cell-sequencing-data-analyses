from celery import current_app as current_celery_app
from celery.result import AsyncResult
from datetime import datetime, timezone

from .celery_config import settings
from constants.declarations import USER_STORAGE
from utils.mongodb import count_jobs_by_status, get_job_from_db, get_jobs_by_status
import os


CELERY_JOB_QUEUES = ("celery", "tools", "benchmarks", "workflows", "web")
RUNNING_JOB_STATUSES = ("Processing", "Running", "STARTED")
QUEUED_JOB_STATUSES = ("Queued", "PENDING")


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


# def get_task_info(job_id):
#     """
#     return task info for the given job_id
#     """
#     task_result = AsyncResult(job_id)
#     summary = "Processing"
#     if task_result.ready():
#         if task_result.result is not None:
#             summary = task_result.result
#         else:
#             summary = task_result.traceback
#     result = {
#         "job_id": job_id,
#         "task_status": task_result.status,
#         "task_result": summary
#     }
#     print(result)
#     return result

#     from celery.result import AsyncResult


def get_task_info(job_id):
    """
    Return task information for the given job_id, checking Celery first,
    then falling back to the database for missing/purged tasks.
    """
    task_result = AsyncResult(job_id)
    
    # 1. Initialize safe defaults to avoid UnboundLocalError
    task_status = task_result.status  # e.g., 'PENDING', 'STARTED', 'SUCCESS'
    summary = "Processing"

    # 2. Handle terminal states (SUCCESS, FAILURE, REVOKED)
    if task_result.ready():
        if task_result.successful():
            summary = task_result.get()
        elif task_result.failed():
            # Safely format the exception and traceback
            traceback_info = task_result.traceback or "No traceback available"
            summary = f"{task_result.result}\n{traceback_info}"
        else:
            # Catch-all for other ready states like 'REVOKED'
            summary = str(task_result.result)

    # 3. Handle PENDING states (Celery defaults to PENDING for unknown IDs)
    elif task_status == 'PENDING':
        backend_key = task_result.backend.get_key_for_task(job_id)
        
        # If the task isn't actually in the Celery backend, check the DB
        if not task_result.backend.get(backend_key):
            results = get_job_from_db(job_id)
            
            # Safely check if 'results' is a valid dictionary before accessing keys
            if results:
                print(f"DB results for job_id {job_id}: {results.get('Status')}")
                # print(f"DB results for job_id {job_id}: {results.get('results')}")
                
                task_status = results.get('Status', task_status)
                summary = results.get('results', "No results available")
            else:
                task_status = "UNKNOWN"
                summary = "Task not found in Celery backend or Database"

    # 4. Return the standard dictionary
    return {
        "job_id": job_id,
        "task_status": task_status,
        "task_result": summary
    }


def _get_broker_queue_counts(celery_app, queue_names=CELERY_JOB_QUEUES):
    counts = {}

    for queue_name in queue_names:
        try:
            with celery_app.connection_for_read() as connection:
                channel = connection.channel()
                try:
                    queue_state = channel.queue_declare(queue=queue_name, passive=True)
                    counts[queue_name] = queue_state.message_count
                finally:
                    channel.close()
        except Exception:
            counts[queue_name] = None

    return counts


def get_celery_job_snapshot(timeout=1.0):
    celery_app = current_celery_app
    queue_counts = _get_broker_queue_counts(celery_app)

    running_count = count_jobs_by_status(RUNNING_JOB_STATUSES)
    db_queued_count = count_jobs_by_status(QUEUED_JOB_STATUSES)
    broker_queued_count = sum(count for count in queue_counts.values() if isinstance(count, int))
    has_broker_counts = any(isinstance(count, int) for count in queue_counts.values())
    queued_count = broker_queued_count if has_broker_counts else db_queued_count

    return {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "summary": {
            "running": running_count,
            "queued": queued_count,
            "broker_queued": broker_queued_count,
            "db_queued": db_queued_count,
        },
        "queues": queue_counts,
        "jobs": {
            "running": get_jobs_by_status(RUNNING_JOB_STATUSES),
            "queued": get_jobs_by_status(QUEUED_JOB_STATUSES),
        },
    }


def get_input_path(input, userID):
    """
    return the absolute input path for a given input
    """
    input_path = None
    if input is not None and userID is not None:
        input_path = USER_STORAGE + userID + input
    return input_path


def get_output(output, userID, job_id):
    """
    return the absolute input path for a given input
    """
    output_path = None
    if output is not None and userID is not None:
        output_path = output + "/" + job_id
    
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
