from starlette.responses import JSONResponse
from fastapi import HTTPException, Body, APIRouter, status
from schemas.schemas import ConversionRequest, ConversionResponse, InputFilesRequest, CombinedQCResult, AnndataMetadata
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata, change_file_extension, get_metadata_from_anndata
from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from typing import List
import logging

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
                default_assay, assay_names, adata_path, annData = run_seurat_qc(file, assay=assay)
                print("API ressss")
                print(annData)
                if assay_names is None:
                    assay_names = []
                result.append({
                        "inputfile": file,
                        "format": "h5seurat",
                        "default_assay": default_assay,
                        "assay_names": assay_names,
                        "adata_path": adata_path
                    })
            else:
                # It's a different file, call load_annData method
                adata = load_anndata(file)
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
                    layers, cell_metadata_obs, umap_coords, gene_metadata, nCells, nGenes, genes, cells, embeddings, traces = get_metadata_from_anndata(scanpy_results)

                    # Return metadata in the API response
                    metadata =  {
                        "layers": layers,
                        "cell_metadata_obs": cell_metadata_obs.to_dict(),
                        "umap_coords": umap_coords.to_dict(),
                        "gene_metadata": gene_metadata.to_dict(),
                        "nCells": nCells,
                        "nGenes": nGenes,
                        "genes": genes,
                        "cells": cells,
                        "embeddings": embeddings,
                        "message": "Quality control completed successfully"
                    }
                    # Return UMAP traces and other metadata in the API response
                    return JSONResponse(content={"traces": traces, "metadata": metadata, "message": "UMAP traces and metadata generated successfully"})

                except Exception as e:
                    logger.exception("Error during Scanpy QC")
                    raise HTTPException(
                        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                        detail=f"Error during Scanpy QC: {str(e)}"
                    )

    except Exception as error:
        logger.exception(f"Error during quality control: {error}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Error during quality control: {str(error)}"
        )

    # If the function reaches this point, it means the quality control process failed
    logger.error("Quality control process failed for unknown reasons")
    raise HTTPException(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        detail="An error occurred during quality control"
    )