from starlette.responses import JSONResponse
from fastapi import HTTPException, Body, APIRouter, status
from schemas.schemas import ConversionRequest, ConversionResponse, InputFilesRequest, CombinedQCResult, AnndataMetadata, DataSplitRequest,BenchmarksRequest
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata, change_file_extension, get_metadata_from_anndata
from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from tools.utils.utils import sc_train_val_test_split
from typing import List
from services.clustering import clustering_task
from tools.visualization.plot import plot_table
import logging
from pathlib import Path
import shutil
import zipfile
import tempfile
from fastapi.encoders import jsonable_encoder
import os

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

router = APIRouter(prefix='/convert', tags=['conversion'], responses={404: {"description": "API Not found"}})


@router.post('/api/convert_to_anndata', response_model=ConversionResponse)
async def convert_to_annData(request_data: ConversionRequest):
    """
    Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object or the list of assay names of Seurat object
    """
    try:
        adata_path, assay_names = convert_seurat_sce_to_anndata(request_data.path)

        if assay_names is None:
            assay_names = []
        if adata_path is None:
            adata_path = "Not available"

        return ConversionResponse(assay_names=assay_names, adata_path=adata_path, message="OK")

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post('/api/convert_sce_to_annData', response_model=dict)
async def receive_data(data: List[dict]):
    response_data = []

    for entry in data:
        path = entry.get('fileDetails')
        assay = entry.get('assayName')

        if path and assay:
            adata_path, assay_names = convert_seurat_sce_to_anndata(path, assay)

            if adata_path and adata_path != None:
                adata_path = adata_path.lstrip('[1] ').rstrip('\n')

            response_data.append({
                'path': path,
                'assay': assay,
                'adata_path': adata_path
            })

    return {'data': response_data, 'message': 'Data processed successfully'}


@router.post('/publishDatasets/validation')
async def process_input_files_validation(request: InputFilesRequest):
    input_files = request.inputFiles
    result = []

    for input in input_files:
        file = input.fileDetails
        assay = input.assay

        try:
            if file.endswith('.h5Seurat') or file.endswith('.h5seurat') or file.endswith('.rds') or file.endswith(".Robj"):
                # It's an H5Seurat or RDS file, call runQCSeurat method
                default_assay, assay_names, adata_path, adata , output= run_seurat_qc(file, assay=assay)
                layers, cell_metadata_obs, gene_metadata, nCells, nGenes, genes, cells, embeddings, umap_plot, violin_plot, scatter_plot, highest_expr_genes_plot = get_metadata_from_anndata(adata)

                 # Return metadata in the API response
                metadata =  {
                    "layers": layers,
                    "cell_metadata_obs": cell_metadata_obs.to_dict(),
                    "gene_metadata": gene_metadata.to_dict(),
                    "nCells": nCells,
                    "nGenes": nGenes,
                    "genes": genes,
                    "cells": cells,
                    "embeddings": embeddings,
                    "message": "Quality control completed successfully"
                }
                if assay_names is None:
                    assay_names = []
                
                result.append({
                        "inputfile": file,
                        "format": "h5seurat",
                        "default_assay": default_assay,
                        "assay_names": assay_names,
                        "adata_path": adata_path,
                        "output": output,
                        "umap_plot": umap_plot,
                        "violin_plot": violin_plot,
                        "scatter_plot": scatter_plot,
                        "highest_expr_genes_plot": highest_expr_genes_plot,
                        "metadata": metadata
                    })
            else:
                # It's a different file, call load_annData method
                adata = load_anndata(file)
                if (os.path.isdir(file)):
                    file = os.path.join(file, "anndata_object.h5ad")

                adata_path = change_file_extension(file, 'h5ad')
                adata.write_h5ad(adata_path)
                result.append({"inputfile": file, "format": "h5ad", "adata_path": adata_path})
        
        except Exception as e:
            # Handle the exception and return an error response
            raise HTTPException(status_code=500, detail=str(e))

    return result


@router.post("/publishDatasets/run/quality_control")
async def run_quality_control(file_mappings: List[dict]):
    try:
        result = []
        for mapping in file_mappings:
            format = mapping.get("format")
            input_path = mapping.get("fileDetails")
            path = mapping.get("adata_path")

            if format == "h5ad":
                # Load the annData object
                adata = load_anndata(path)

                # Run Scanpy QC
                try:
                    scanpy_results = run_scanpy_qc(adata)
                    layers, cell_metadata_obs, gene_metadata, nCells, nGenes, genes, cells, embeddings, umap_plot, violin_plot, scatter_plot, highest_expr_genes_plot = get_metadata_from_anndata(scanpy_results)

                    # Return metadata in the API response
                    metadata =  {
                        "layers": layers,
                        "cell_metadata_obs": cell_metadata_obs.to_dict(),
                        "gene_metadata": gene_metadata.to_dict(),
                        "nCells": nCells,
                        "nGenes": nGenes,
                        "genes": genes,
                        "cells": cells,
                        "embeddings": embeddings,
                        "message": "Quality control completed successfully"
                    }
                    
                    result.append({
                        "inputfile": input_path,
                        "format": format,
                        "adata_path": path,
                        "umap_plot": umap_plot,
                        "violin_plot": violin_plot,
                        "scatter_plot": scatter_plot,
                        "highest_expr_genes_plot": highest_expr_genes_plot,
                        "metadata": metadata
                    })
                except Exception as e:
                    logger.exception("Error during Scanpy QC")
                    raise HTTPException(
                        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail=f"Error during Scanpy QC: {str(e)}"
                    )
        return result

    except Exception as error:
        logger.exception(f"Error during quality control: {error}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error during quality control: {str(error)}"
        )

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
        return {"result": "Data split successful", "archive_path": archive_path}
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
