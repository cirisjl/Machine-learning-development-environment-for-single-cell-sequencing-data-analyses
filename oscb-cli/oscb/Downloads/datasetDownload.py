import os
import time
import requests
from tqdm import tqdm
import platform
import hashlib
from pathlib import Path
import websockets
import asyncio
import re



def downloadDataset(dataset_id, destination_path, process_type, method):
    # Define the base URL for the API
    base_url = "http://130.127.133.171:5005/api"

    ws_base_url = "ws://130.127.133.171:5005/wsapi"  # WebSocket base URL


    user_id = get_persistent_machine_id()
    # Step 1: Submit the task
    submit_url = f"{base_url}/oscb-cli/downloadDataset"
    payload = {
        "dataset_id": dataset_id,
        "user_id": user_id,
        "process_type":process_type,
        "method":method
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

    # Start a coroutine to fetch logs in real time
    asyncio.run(fetch_logs_from_websocket(ws_base_url, job_id))


    # Step 2: Poll for task status
    task_status_url = f"{base_url}/job/downloadDataset/{job_id}"
    filename = None  # Initialize task_result

    while True:
        try:
            # Make a GET request to check task status
            response = requests.get(task_status_url)

            # Check if the file is ready for download
            if response.headers.get("content-disposition"):
                content_disposition = response.headers.get("content-disposition")
                filename = content_disposition.split("filename=")[1].strip('"')
                print("Task completed. File ready for download.")
                break  # Move to download step
            
            # Parse the JSON response
            task_result = response.json() if response.text.strip() else None

            if not task_result:
                print("No valid task response from the server. Exiting!")
                return

            # Handle task statuses
            job_status = task_result.get("job_status")
            status = task_result.get("status")

            if job_status == "SUCCESS":
                if status == "Task result does not contain a file path":
                    print("Task completed, but no file to download.")
                if status == "Invalid Response for the Task Submitted":
                    print("Invalid Response for the Task Submitted")
                return
            elif job_status == "FAILURE":
                error_message = task_result.get("error_message", "Unknown error occurred.")
                print(f"Task failed: {error_message}")
                return
            else:
                print("Task is still processing. Waiting...")
                time.sleep(5)  # Wait before polling again

        except requests.RequestException as e:
            print(f"Error occurred while checking task status: {e}")
            return

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

def get_persistent_machine_id():
    # Generate a new unique ID based on system properties
    system_properties = f"{platform.node()}-{platform.system()}-{platform.processor()}-{platform.machine()}"
    hashed_id = hashlib.sha256(system_properties.encode()).hexdigest()
    return hashed_id

async def fetch_logs_from_websocket(base_url, job_id):
    """
    Connect to WebSocket and display logs for the given job ID in real-time.
    """
    ws_url = f"{base_url}/log/{job_id}"
    try:
        async with websockets.connect(ws_url) as websocket:
            print(f"Connected to WebSocket for job ID: {job_id}. Receiving logs...\n")
            async for message in websocket:
                # Remove HTML tags from the message
                clean_message = re.sub(r'<.*?>', '', message)
                
                # Apply color and print each message on a new line
                if "ERROR" in clean_message:
                    print(f'\033[91m{clean_message}\033[0m')  # Red for ERROR
                elif "WARNING" in clean_message:
                    print(f'\033[93m{clean_message}\033[0m')  # Yellow for WARNING
                elif "SUCCESS" in clean_message:
                    print(f'\033[92m{clean_message}\033[0m')  # Green for SUCCESS
                else:
                    print(f'{clean_message}\n')  # Ensure each message is printed on a new line
    except websockets.exceptions.ConnectionClosed as e:
        print(f"WebSocket connection closed unexpectedly: {e}")
    except Exception as e:
        print(f"Error occurred while fetching logs from WebSocket: {e}")

if __name__ == "__main__":
    dataset_id = "U-h-Heart-Wang-2024@kbcfh"
    destination_path = "datasets"
    process_type="quality_control"
    method = "scanpy"

    downloadDataset(dataset_id, destination_path, process_type, method)

    # # Example usage
    # machine_id = get_persistent_machine_id()
    # print(f"Machine ID: {machine_id}")
