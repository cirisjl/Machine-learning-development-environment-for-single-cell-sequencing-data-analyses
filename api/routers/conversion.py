from starlette.responses import JSONResponse
from fastapi import HTTPException, Body, APIRouter
from schemas.schemas import ConversionRequest, ConversionResponse, InputFilesRequest, CombinedQCResult, AnndataMetadata
from tools.formating.formating import convert_seurat_sce_to_anndata, load_anndata, change_file_extension, get_metadata_from_anndata
from tools.qc.scanpy_qc import run_scanpy_qc
from tools.qc.dropkick_qc import run_dropkick_qc
from tools.qc.seurat_qc import run_seurat_qc
from typing import List

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
                default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap, adata_path = run_seurat_qc(file, assay=assay)
                if assay_names is None:
                    assay_names = []
                result.append({
                        "inputfile": file,
                        "format": "h5seurat",
                        "default_assay": default_assay,
                        "assay_names": assay_names,
                        "metadata": metadata,
                        "nCells": nCells,
                        "nGenes": nGenes,
                        "genes": genes,
                        "cells": cells,
                        "HVGsID": HVGsID,
                        "pca": pca,
                        "tsne": tsne,
                        "umap": umap,
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
    qc_results = []

    try:
        for mapping in file_mappings:
            format = mapping.get("format")
            input_path = mapping.get("fileDetails")
            path = mapping.get("adata_path")

            if format == "seurat":
                
                default_assay, assay_names, metadata, nCells, nGenes, genes, cells, HVGsID, pca, tsne, umap, adata_path = run_seurat_qc(input_path, assay)
                qc_results.append({
                    "inputfile": input_path,
                    "format": "h5seurat",
                    "default_assay": default_assay,
                    "assay_names": assay_names,
                    "metadata": metadata,
                    "nCells": nCells,
                    "nGenes": nGenes,
                    "genes": genes,
                    "cells": cells,
                    "HVGsID": HVGsID,
                    "pca": pca,
                    "tsne": tsne,
                    "umap": umap,
                    "adata_path": adata_path
                })

            elif format == "h5ad":

                # Load the annData object
                adata = load_anndata(path)

                # Run Scanpy QC
                try:
                    scanpy_results = run_scanpy_qc(adata)
                    layers, cell_metadata_obs, cell_metadata_obsm, gene_metadata, nCells, nGenes, genes, cells, embeddings = get_metadata_from_anndata(scanpy_results)
                    scanpy_metadata = AnndataMetadata(
                        layers=layers,
                        cell_metadata_obs=cell_metadata_obs.to_dict(),
                        cell_metadata_obsm=cell_metadata_obsm.to_dict(),
                        gene_metadata=gene_metadata.to_dict(),
                        nCells=nCells,
                        nGenes=nGenes,
                        genes=genes,
                        cells=cells,
                        embeddings=embeddings
                    )
                except Exception as e:
                    print("Scanpy QC failed")
                    print(e)

                # Run Dropkick QC
                # try:
                #     dropkick_results = run_dropkick_qc(adata)
                #     layers, cell_metadata_obs, cell_metadata_obsm, gene_metadata, nCells, nGenes, genes, cells, embeddings = get_metadata_from_anndata(dropkick_results)
                #     dropkick_metadata = AnndataMetadata(
                #         layers=layers,
                #         cell_metadata_obs=cell_metadata_obs.to_dict(),
                #         cell_metadata_obsm=cell_metadata_obsm.to_dict(),
                #         gene_metadata=gene_metadata.to_dict(),
                #         nCells=nCells,
                #         nGenes=nGenes,
                #         genes=genes,
                #         cells=cells,
                #         embeddings=embeddings
                #     )
                # except Exception as e:
                #     print("DropKick QC failed")
                #     print(e)

            # # Append combined metadata to qc_results
            #     qc_results.append({
            #         "inputfile": input_path,
            #         "format": "h5ad",
            #         "combined_results": CombinedQCResult(
            #             scanpy_results=scanpy_metadata,
            #             # dropkick_results=dropkick_metadata
            #         ).dict()
            #     })

    except Exception as error:
        raise HTTPException(status_code=500, detail=f"An error occurred during quality control: {str(error)}")

    return JSONResponse(content=scanpy_metadata, status_code=200)
