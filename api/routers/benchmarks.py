from starlette.responses import JSONResponse
from fastapi import HTTPException, Body, APIRouter, status
from schemas.schemas import ConversionRequest, ConvertRequest, ConversionResponse, InputFilesRequest, CombinedQCResult, AnndataMetadata, DataSplitRequest, SubsetDataRequest, BenchmarksRequest, QualityControlRequest, UMAPRequest
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata, change_file_extension, get_metadata_from_anndata, get_metadata_from_seurat, get_md5
from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from tools.utils.datasplit import sc_train_val_test_split, subset_by_obskey
from typing import List
from services.clustering import clustering_task
from tools.visualization.plot import plot_table
from utils.unzip import unzip_file_if_compressed
from pathlib import Path
import shutil
import tempfile
from fastapi.encoders import jsonable_encoder
import os


router = APIRouter(prefix='/convert', tags=['benchmarks'], responses={404: {"description": "API Not found"}})


# @router.post('/api/convert_to_anndata', response_model=ConversionResponse)
# async def convert_to_annData(request_data: ConversionRequest):
#     """
#     Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object or the list of assay names of Seurat object
#     """
#     try:
#         file_path = unzip_file_if_compressed(request_data.path)
#         adata_path, assay_names, default_assay = convert_seurat_sce_to_anndata(file_path)

#         if assay_names is None:
#             assay_names = []
#         if adata_path is None:
#             adata_path = "Not available"

#         return ConversionResponse(assay_names=assay_names, adata_path=adata_path, message="OK")

#     except Exception as e:
#         raise HTTPException(status_code=500, detail=str(e))

# @router.post('/api/convert_sce_to_annData', response_model=dict)
# async def receive_data(data: List[dict]):
#     response_data = []

#     for entry in data:
#         path = unzip_file_if_compressed(entry.get('fileDetails'))
#         assay = entry.get('assayName')

#         if path and assay:
#             adata_path, assay_names, default_assay = convert_seurat_sce_to_anndata(path, assay)

#             if adata_path and adata_path != None:
#                 adata_path = adata_path.lstrip('[1] ').rstrip('\n')

#             response_data.append({
#                 'path': path,
#                 'assay': assay,
#                 'adata_path': adata_path
#             })

#     return {'data': response_data, 'message': 'Data processed successfully'}


@router.post('/publishDatasets/validation')
async def process_input_files_validation(request: InputFilesRequest):
    input_files = request.inputFiles
    result = []

    for input in input_files:
        file = unzip_file_if_compressed(input.fileDetails)
        assay = input.assay

        try:
            if file.endswith('.h5Seurat') or file.endswith('.h5seurat') or file.endswith('.rds') or file.endswith(".Robj"):
                # It's an H5Seurat or RDS file, call runQCSeurat method
                info, default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap = get_metadata_from_seurat(file)

                if assay_names is None:
                    assay_names = []
                
                result.append({
                        "inputfile": file,
                        "info": info,
                        "format": "h5seurat",
                        "default_assay": default_assay,
                        "assay_names": assay_names
                    })
        except Exception as e:
            # Handle the exception and return an error response
            print(str(e))
            raise HTTPException(status_code=500, detail=str(e))        

    return result


@router.post("/publishDatasets/run/quality_control")
async def run_quality_control(file_mappings: QualityControlRequest):
    try:
        result = []
        print(file_mappings)
        assay = file_mappings.assay
        doublet_rate = file_mappings.doublet_rate
        input_path = file_mappings.fileDetails
        min_genes = file_mappings.min_genes
        max_genes = file_mappings.max_genes
        min_cells = file_mappings.min_cells
        target_sum = file_mappings.target_sum
        n_top_genes = file_mappings.n_top_genes
        n_neighbors = file_mappings.n_neighbors
        n_pcs = file_mappings.n_pcs
        resolution = file_mappings.resolution
        regress_cell_cycle = file_mappings.regress_cell_cycle
        use_default = file_mappings.use_default

        md5 = get_md5(input_path)
        pp_stage = "Raw"
        method_id = None
        parameters = {
            "min_genes": min_genes, # default: 200, step: 25, range: [0, 20000], scale: 200(default), 1000, 5000, 10000, 15000, 20000(No limit)
            "max_genes": max_genes, # default: 20000(=No limit, default), step: 25, range: [0, 20000], scale: 200, 1000, 5000, 10000, 15000, 20000(=No limit, default)
            "min_cells": min_cells, # default: 2, step:1, range: [1, 200], scale: 2(default), 10, 50, 100, 200
            # "max_cells": max_cells,
            "target_sum": target_sum, # default: 0(None), step: 1e4, range:[0, 1e6], scale: 0(None, default), 1e4, 1e5, 1e6
            "n_top_genes": n_top_genes, # highly variable genes:  default: 2000, step:25, range: [100, 10000], scale: 500, 1000, 2000(default), 5000, 10000
            "n_neighbors": n_neighbors, # default: 15, step:1, range: [2, 100], scale: 1, 5, 10, 15(default), 20, 50, 100
            "n_pcs": n_pcs, # default: 0(None), step:1, range: [0, 200], scale: 0(None, default), 5, 10, 20, 40, 50, 125, 200 
            "resolution": resolution, # default: 1, step:0.05, range: [0, 5], scale: 0, 0.1, 0.25, 0.5, 1(default), 2.5, 5
            "doublet_rate": doublet_rate, # default: 0.08, step:0.001, range: [0, 0.5], scale: 0, 0.8%, 2.3%, 3.8%, 4.6%, 6.1%, 8%(default), 12.5%, 20_, 50% Please show the scale in the form of percentage
            "regress_cell_cycle": regress_cell_cycle, # default: false, values: [true, false]
            "use_default": use_default # default: true, values: [true, false]
        }

        if input_path.endswith('.h5Seurat') or input_path.endswith('.h5seurat') or input_path.endswith('.rds') or input_path.endswith(".Robj"):
            # It's an H5Seurat or RDS file, call runQCSeurat method
            default_assay, assay_names, adata_path, adata, output, ddl_assay_names = run_seurat_qc(input_path, assay=assay, min_genes=200, max_genes=0, min_UMI_count=2, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=0.5, dims=10, doublet_rate=0.075, regress_cell_cycle=False)
            # default_assay, assay_names, adata_path, adata, output= run_seurat_qc(input_path, assay=assay, min_genes=min_genes, max_genes=max_genes, min_UMI_count=min_cells, max_UMI_count=0, percent_mt_max=5, percent_rb_min=0, resolution=resolution, dims=n_neighbors, regress_cell_cycle=regress_cell_cycle)
            
            if ddl_assay_names:
                result.append({
                    "inputfile": input_path,
                    "format": "h5seurat",
                    "default_assay": default_assay,
                    "assay_names": assay_names,
                    "ddl_assay_names": ddl_assay_names
                })

                return result
            
            info, layers, cell_metadata_obs, gene_metadata, nCells, nGenes, genes, cells, embeddings, umap_plot, umap_plot_3d, violin_plot, scatter_plot, highest_expr_genes_plot = get_metadata_from_anndata(adata)

            if(use_default):
                method_id = "seurat_qc"
            else:
                method_id = "seurat_qc" + "-" + min_genes + "-" + max_genes + "-" + min_cells + "-" + target_sum + "-" + n_top_genes + "-" + n_neighbors + "-" + n_pcs + "-" + resolution + "-" + regress_cell_cycle

            pp_results = {
                "stage": pp_stage,
                "task": "QC",
                "method": "seurat",
                "method_id": method_id,
                "parameters": parameters,
                "files": adata_path
            }

            # Return metadata in the API response
            metadata =  {
                "layers": layers,
                "cell_metadata_obs": cell_metadata_obs.to_dict(),
                "gene_metadata": gene_metadata.to_dict(),
                "nCells": nCells,
                "nGenes": nGenes,
                "genes": genes,
                "cells": cells,
                "embeddings": embeddings
            }
            if assay_names is None:
                assay_names = []
            
            result.append({
                    "inputfile": input_path,
                    "info": info,
                    "format": "h5seurat",
                    "default_assay": default_assay,
                    "assay_names": assay_names,
                    "adata_path": adata_path,
                    "output": output,
                    "umap_plot": umap_plot,
                    "umap_plot_3d": umap_plot_3d,
                    "violin_plot": violin_plot,
                    "scatter_plot": scatter_plot,
                    "highest_expr_genes_plot": highest_expr_genes_plot,
                    "md5": md5,
                    "metadata": metadata,
                    "pp_results": pp_results,
                    "message": "Quality control completed successfully."
                })
        else:
            # Load the annData object
            adata = load_anndata(input_path)

            # Run Scanpy QC
            try:
                scanpy_results = run_scanpy_qc(adata, min_genes=200, max_genes=None, min_cells=2, target_sum=1e4, n_top_genes=None, n_neighbors=15, n_pcs=None, resolution=1, regress_cell_cycle=False)
                # scanpy_results = run_scanpy_qc(adata, min_genes=min_genes, max_genes=max_genes, min_cells=min_cells, target_sum=target_sum, n_top_genes=n_top_genes, n_neighbors=n_neighbors, n_pcs=n_pcs, resolution=resolution, regress_cell_cycle=regress_cell_cycle)
                info, layers, cell_metadata_obs, gene_metadata, nCells, nGenes, genes, cells, embeddings, umap_plot, umap_plot_3d, violin_plot, scatter_plot, highest_expr_genes_plot = get_metadata_from_anndata(scanpy_results)

                adata_path = change_file_extension(input_path, 'h5ad')
                scanpy_results.write_h5ad(adata_path)

                if(use_default):
                    method_id = "scanpy_qc"
                else:
                    method_id = "scanpy_qc" + "-" + min_genes + "-" + max_genes + "-" + min_cells + "-" + target_sum + "-" + n_top_genes + "-" + n_neighbors + "-" + n_pcs + "-" + resolution + "-" + regress_cell_cycle

                pp_results = {
                    "stage": pp_stage,
                    "task": "QC",
                    "method": "scanpy",
                    "method_id": method_id,
                    "parameters": parameters,
                    "files": adata_path
                }

                # Return metadata in the API response
                metadata =  {
                    "layers": layers,
                    "cell_metadata_obs": cell_metadata_obs.to_dict(),
                    "gene_metadata": gene_metadata.to_dict(),
                    "nCells": nCells,
                    "nGenes": nGenes,
                    "genes": genes,
                    "cells": cells,
                    "embeddings": embeddings
                }
                
                result.append({
                    "inputfile": input_path,
                    "info": info,
                    "format": "h5ad",
                    "adata_path": adata_path,
                    "umap_plot": umap_plot,
                    "umap_plot_3d": umap_plot_3d,
                    "violin_plot": violin_plot,
                    "scatter_plot": scatter_plot,
                    "highest_expr_genes_plot": highest_expr_genes_plot,
                    "md5": md5,
                    "metadata": metadata,
                    "pp_results": pp_results,
                    "message": "Quality control completed successfully."
                    
                })
            except Exception as e:
                # logger.exception("Error during Scanpy QC")
                print(e)
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail=f"Error during Scanpy QC: {str(e)}"
                )
        return result

    except Exception as error:
        print(f"Error during quality control: {error}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error during quality control: {str(error)}"
        )


router.post("/api/umap")
async def data_split(user_request: UMAPRequest):
    from tools.visualization.plot import plot_UMAP
    try:
        # Access user data
        adata_path = user_request.adata_path
        layer = user_request.layer
        clustering_plot_type = user_request.clustering_plot_type
        selected_cell_intersection = user_request.selected_cell_intersection
        n_dim = user_request.n_dim

        adata = load_anndata(adata_path)
        umap_json = plot_UMAP(adata, layer=layer, clustering_plot_type=clustering_plot_type, selected_cell_intersection=selected_cell_intersection, n_dim=n_dim)
    
        return umap_json
    except Exception as e:
        # Handle any errors
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/api/data-split")
async def data_split(user_data: DataSplitRequest):
    try:
        # Access user data
        data_filepath = user_data.data
        train_fraction = user_data.train_fraction
        validation_fraction = user_data.validation_fraction
        test_fraction = user_data.test_fraction

        adata = load_anndata(data_filepath)

        train, validation, test = sc_train_val_test_split(adata, train_fraction, validation_fraction, test_fraction)
       
       # Extract directory and filename from the data filepath
        data_directory = Path(data_filepath).parent
        data_filename = Path(data_filepath).stem

        # Define a temporary directory to store the files
        temp_dir = tempfile.TemporaryDirectory(dir=data_directory)

        # Write AnnData objects to files with unique filenames in the temporary directory
        train.write(Path(temp_dir.name) / f"{data_filename}_train.h5ad")
        validation.write(Path(temp_dir.name) / f"{data_filename}_validation.h5ad")
        test.write(Path(temp_dir.name) / f"{data_filename}_test.h5ad")

        # Compress files into a single archive in the same directory
        shutil.make_archive(data_directory / f"{data_filename}_data_split", 'zip', temp_dir.name)

        # Return the path to the compressed archive
        archive_path = data_directory / f"{data_filename}_data_split.zip"
        return {"result": "Data split successfully.", "archive_path": archive_path}
    except Exception as e:
        # Handle any errors
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/api/subset-data")
async def subset_data(user_data: SubsetDataRequest):
    try:
        # Access user data
        data_filepath = user_data.data
        obskey = user_data.obskey
        values = user_data.values

        adata = load_anndata(data_filepath)

        adata_sub = subset_by_obskey(adata, obskey, values)
       
       # Extract directory and filename from the data filepath
        data_directory = Path(data_filepath).parent
        data_filename = Path(data_filepath).stem

        # Return the path to the compressed archive
        archive_path = data_directory / f"{data_filename}_sub.h5ad"
        adata_sub.write(archive_path)

        return {"result": "AnnData is subset successfully.", "archive_path": archive_path}
    except Exception as e:
        # Handle any errors
        raise HTTPException(status_code=500, detail=str(e))
    

@router.post("/publishDatasets/benchmarks")
async def process_task_data(data: BenchmarksRequest):
    try:
        # Access the data received
        task_type = data.task_type
        items = data.data
        results = []
        
        for item in items:
            if task_type.lower() == 'clustering':  # Check if task_type is 'clustering'
                adata_path = item.adata_path
                task_label = item.task_label
                datasetId = item.datasetId
                clustering_results = clustering_task(adata_path, task_label, datasetId, task_type)
                results.append(clustering_results)
       
        return results

    except Exception as e:
        # Handle exceptions as needed
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")
    

@router.post("/api/getTablePlot")
async def process_files(file_paths: List[str]):
    
    try:
        results = []
        for file_path in file_paths:
            adata = load_anndata(file_path)
            if adata is not None:
                # Convert AnnData to DataFrame
                dataframe = adata.to_df()
                tablePlot = plot_table(dataframe)
                results.append(tablePlot)
        return results
    except Exception as e:
        # Handle exceptions as needed
        raise HTTPException(status_code=500, detail=f"Internal Server Error: {str(e)}")


@router.post('/api/to_adata_or_srat')
async def convert_files_to_adata_or_srat(request_body: ConvertRequest):

    fileDetails = request_body.fileDetails
    assay_name = request_body.assay_name
    results = []

    # Check if only one file is provided
    if len(fileDetails) == 1:
        file = unzip_file_if_compressed(fileDetails[0])
        try:
            # Check if the file is of a specific type
            if file.endswith(('.h5Seurat', 'h5seurat', 'rds', '.Robj')):
                # Process as specified for these file types
                adata_path, assay_names, default_assay = convert_seurat_sce_to_anndata(file, assay=assay_name)
                results.append({
                    'adata_path': adata_path,
                    'assay_names': assay_names,
                    'default_assay': default_assay,
                    "inputfile": fileDetails,
                    "format": "h5seurat"
                })
            else:
                adata_path = change_file_extension(file, 'h5ad')
                adata = load_anndata(file)
                adata.write_h5ad(adata_path)
                results.append({"inputfile": fileDetails, "adata_path": adata_path, "format": "h5ad"})
        except Exception as e:
            # Handle the exception and return an error response
            raise HTTPException(status_code=500, detail=str(e))

    elif len(fileDetails) > 1:
        parent_directory = os.path.dirname(unzip_file_if_compressed(fileDetails[0]))
        
        # Optionally, verify that all files are in the same directory
        if not all(os.path.dirname(unzip_file_if_compressed(file)) == parent_directory for file in fileDetails):
            raise HTTPException(status_code=400, detail="Not all files are in the same directory.")

        try:
            # Now, use the parent directory to load the dataset
            adata = load_anndata(parent_directory)
            adata_path = os.path.join(parent_directory, "combined_dataset.h5ad")
            adata.write_h5ad(adata_path)
            results.append({"inputfile": fileDetails, "adata_path": adata_path, "format": "h5ad"})
        except Exception as e:
            raise HTTPException(status_code=500, detail=str(e))
 
    return results