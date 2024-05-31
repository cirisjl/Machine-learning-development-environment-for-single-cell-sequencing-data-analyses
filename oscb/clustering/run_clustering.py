import subprocess

# Define the path to your .h5ad file
file_path = r"C:\Users\sindh\oscb\oscb\download_dataset\datasets\hca_heart_neuronal_raw (2).h5ad"

# Call the clustering script with the file path
subprocess.run(['python', 'clustering_analysis.py', file_path])
