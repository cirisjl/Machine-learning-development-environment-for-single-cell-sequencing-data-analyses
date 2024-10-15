import os
from tools.formating.formating import load_anndata, convert_seurat_sce_to_anndata, load_anndata, change_file_extension, get_cell_metadata
from utils.unzip import unzip_file_if_compressed
from exceptions.custom_exceptions import CeleryTaskException


def load_metadata(job_id, file_dict):
    fileDetails = file_dict['fileDetails']
    userID = file_dict['userID']
    assay_name = file_dict['assay_name']
    results = {}
    adata = None
    adata_path = None
    seurat_path = None

    # Check if only one file is provided
    if len(fileDetails) == 1:
        file = unzip_file_if_compressed(job_id, fileDetails[0])
        format = None
        try:
            results['inputfile'] = fileDetails
            adata = load_anndata(file)
            if file.endswith('.h5ad'):
                adata_path = file
                format = "h5ad"
            else:
                adata_path = change_file_extension(fileDetails[0], 'h5ad')
                adata.write_h5ad(adata_path, compression='gzip')
            results["adata_path"] = adata_path
    
            if file.endswith(('.h5Seurat', 'h5seurat')):
                seurat_path = file
                format = "h5seurat"
                results['seurat_path'] = seurat_path

        except Exception as e:
            # Handle the exception and return an error response
            detail = f"Error during loading dataset: {str(e)}"
            raise CeleryTaskException(detail)

    elif len(fileDetails) > 1:
        parent_directory = os.path.dirname(unzip_file_if_compressed(job_id, fileDetails[0]))
        
        # Optionally, verify that all files are in the same directory
        if not all(os.path.dirname(unzip_file_if_compressed(job_id, file)) == parent_directory for file in fileDetails):
            detail="Not all files are in the same directory."
            raise CeleryTaskException(detail)

        try:
            # Now, use the parent directory to load the dataset
            adata = load_anndata(parent_directory)
            adata_path = os.path.join(parent_directory, ".h5ad")
            adata.write_h5ad(adata_path, compression='gzip')
            results["inputfile"] = fileDetails
            results["adata_path"] = adata_path

        except Exception as e:
            detail=str(e)
            raise CeleryTaskException(detail)
    
    cell_metadata, cell_metadata_head, obs_names, nCells, nGenes, layers, info, adata_size, embeddings = get_cell_metadata(adata, adata_path=adata_path)
    adata = None
    results['cell_metadata'] = cell_metadata
    results['cell_metadata_head'] = cell_metadata_head
    results['obs_names'] = obs_names
    results['nCells'] = nCells
    results['nGenes'] = nGenes
    results['layers'] = layers
    results['info'] = info
    results['adata_size'] = adata_size
    results['embeddings'] = embeddings

    return results