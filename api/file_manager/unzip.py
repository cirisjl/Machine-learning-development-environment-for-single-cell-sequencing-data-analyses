import os
import sys
import zipfile
import rarfile
import tarfile
from tarfile import is_tarfile

# zip files
def extract_zip_all(zip_file_path):
    extract_path = os.path.dirname(zip_file_path)
    try:
        with zipfile.ZipFile(zip_file_path,'r') as zip_ref:
            file_list = zip_ref.namelist()
            zip_ref.extractall(extract_path)
            if len(file_list) == 1:
                extract_path = os.path.dirname(zip_file_path) + "/" + file_list[0]
                
            os.remove(zip_file_path) # Delete the zip file
    except Exception as e:
        print('An error occurred:', str(e))
        return None
        
    return extract_path

# rar files
def extract_rar_all(rar_file_path):
    extract_path = os.path.dirname(rar_file_path)
    try:
        with rarfile.RarFile(rar_file_path,'r') as rar_ref:
            file_list = rar_ref.namelist()
            rar_ref.extractall(extract_path)
            if len(file_list) == 1:
                extract_path = os.path.dirname(rar_file_path) + "/" + file_list[0]
                
            os.remove(rar_file_path) # Delete the zip file
    except Exception as e:
        print('An error occurred:', str(e))
        return None
        
    return extract_path
    
# tar, tar.gz, bz2, xz files
def extract_tar_all(gz_file_path):
    extract_path = os.path.dirname(gz_file_path)
    try:
        with tarfile.open(gz_file_path,'r') as tar_ref:
            file_list = tar_ref.getnames()
            tar_ref.extractall(extract_path)
            extract_path = os.path.dirname(gz_file_path) + "/" + file_list[0]
            subfolder = traversal_subfolder(extract_path)
            if subfolder is not None or len(subfolder) != 0:
                extract_path += '/' + subfolder[0]
                
            os.remove(gz_file_path) # Delete the zip file
    except Exception as e:
        print('An error occurred:', str(e))
        return None
        
    return extract_path

# Find the deepest subfolder
def traversal_subfolder(path):
    list = []
    if (os.path.exists(path)):
        files = os.listdir(path)
        for file in files:
            m = os.path.join(path, file)
            if (os.path.isdir(m)):
                h = os.path.split(m)
                list.append(h[1])
        return list

# Unzip file if it's compressed, toherwise return its original path
def unzip_file_if_compressed(file_path):
    extract_path = file_path
    if os.path.exists(file_path):
        if file_path.endswith(".zip"):
            extract_path = extract_zip_all(file_path)
        elif file_path.endswith(".rar"):
            extract_path = extract_rar_all(file_path)
        elif is_tarfile(file_path) or file_path.endswith(".gz"):
            extract_path = extract_tar_all(file_path)
        else:
            print("No compressed file is found in: " + file_path)

    return extract_path