import os
import requests
from urllib.parse import unquote
from tqdm import tqdm
import re

def download_file_by_id(id):
    try:
        # Define the base URL where your server is hosted
        base_url = 'http://node-0.sindhuja1.biomizzou-pg0.wisc.cloudlab.us:3001'

        # Define the directory where you want to save the downloaded files
        download_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'datasets')

        # Make sure the 'datasets' directory exists, if not, create it
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

        # Make a HEAD request to check if the file exists
        response = requests.head(f'{base_url}/download2/{id}')

        # Check if the server responded with a success status code
        if response.status_code == 200:
            # Extract the content disposition from the response headers
            content_disposition = response.headers.get('content-disposition')

            # Extract the filename from the content disposition header, if available
            filename = id
            if content_disposition:
                filename = unquote(content_disposition.split("filename=")[1])

            # Remove invalid characters from the filename
            filename = re.sub(r'[<>:"/\\|?*]', '', filename)

            # Define the path where you want to save the downloaded file
            download_path = os.path.join(download_dir, filename)

            # Check if the file already exists in the download directory
            if os.path.exists(download_path):
                print('File already exists:', download_path)
                return

            # Make a GET request to download the file
            with requests.get(f'{base_url}/download2/{id}', stream=True) as response:
                # Initialize tqdm to show download progress
                total_size = int(response.headers.get('content-length', 0))
                progress_bar = tqdm(total=total_size, unit='B', unit_scale=True, desc=filename)
                with open(download_path, 'wb') as file:
                    for data in response.iter_content(chunk_size=1024):
                        file.write(data)
                        progress_bar.update(len(data))

            print('\nFile downloaded successfully:', download_path)
        else:
            print('Error downloading file:', response.reason)
    except Exception as e:
        print('Error downloading file:', e)