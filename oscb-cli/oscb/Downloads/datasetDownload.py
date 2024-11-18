import os
import time
import requests
from tqdm import tqdm

def downloadDataset(dataset_id, destination_path, process_type="Normalization", method=""):
    # Define the base URL for the API
    base_url = "http://130.127.133.108:5005/api"

    # Step 1: Submit the task
    submit_url = f"{base_url}/oscb-cli/downloadDataset"
    payload = {
        "dataset_id": dataset_id
    }

    response = requests.post(submit_url, json=payload)

    if response.status_code != 200:
        print(f"Error submitting the task: {response.status_code}, {response.text}")
        return

    task_info = response.json()
    job_id = task_info.get("job_id")
    if not job_id:
        print("Task submission failed: Missing job ID.")
        return

    print(f"Task submitted successfully. Job ID: {job_id}")

    # Step 2: Poll for task status
    task_status_url = f"{base_url}/job/downloadDataset/{job_id}"
    filename = None  # Initialize task_result

    while True:
        status_response = requests.get(task_status_url)
 
        # Check if the file is ready for download
        if status_response.headers.get("content-disposition"):
            content_disposition = status_response.headers.get('content-disposition')
            filename = content_disposition.split('filename=')[1].strip('"')
            print("Task completed. File ready for download.")
            break # Move to download step

        # Parse the JSON response if it's not empty
        if status_response.text.strip():
            task_result = status_response.json()
            print("Task results")
            print(task_result)
            print(task_result["job_status"])
            if "status" in task_result and task_result["status"] == "Task result does not contain a file path" and task_result["job_status"] == "PROCESSING":
                print("Task completed, but no file to download.")
                return

        print("Task is still processing. Waiting...")
        time.sleep(2)  # Wait before polling again

    # Step 3: Download the file
    download_dir = os.path.abspath(destination_path)

    # Make sure the destination directory exists
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    download_path = os.path.join(download_dir, filename)

    # Now perform the download
    file_response = requests.get(task_status_url, stream=True)

    if file_response.status_code == 200:
        total_size = int(file_response.headers.get('content-length', 0))
        progress_bar = tqdm(total=total_size, unit="B", unit_scale=True, desc=filename)
        with open(download_path, 'wb') as file:
            for chunk in file_response.iter_content(chunk_size=1024):
                file.write(chunk)
                progress_bar.update(len(chunk))

        print(f"\nFile downloaded successfully: {download_path}")
    else:
        print(f"Error downloading the file: {file_response.status_code}, {file_response.text}")

if __name__ == "__main__":
    dataset_id = "U-h-Heart-Wang-2024@kbcfh"
    destination_path = "datasets"
    downloadDataset(dataset_id, destination_path)
