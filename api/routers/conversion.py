from starlette.responses import JSONResponse
from fastapi import HTTPException, Body, APIRouter
from schemas.schemas import ConversionRequest, ConversionResponse, InputFilesRequest
from tools.formating.formating import ConvertSeuratSCEtoAnndata, run_seurat_qc, LoadAnndata, change_file_extension
from typing import List

router = APIRouter(prefix='/convert', tags=['conversion'], responses={404: {"description": "API Not found"}})


@router.post('/api/convert_to_anndata', response_model=ConversionResponse)
async def convert_to_annData(request_data: ConversionRequest):
    """
    Convert Seurat/Single-Cell Experiment object to Anndata object and return the path of Anndata object or the list of assay names of Seurat object
    """
    try:
        adata_path, assay_names = ConvertSeuratSCEtoAnndata(request_data.path)

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
            adata_path, assay_names = ConvertSeuratSCEtoAnndata(path, assay)

            if adata_path and adata_path != None:
                adata_path = adata_path.lstrip('[1] ').rstrip('\n')

            response_data.append({
                'path': path,
                'assay': assay,
                'adata_path': adata_path
            })

    return {'data': response_data, 'message': 'Data processed successfully'}


@router.post('/publishDatasets/validation')
async def process_input_files(request: InputFilesRequest):
    input_files = request.inputFiles
    result = []

    for file in input_files:
        try:
            if file.endswith('.h5Seurat') or file.endswith('.h5seurat') or file.endswith('.rds') or file.endswith(".Robj"):
                # It's an H5Seurat or RDS file, call runQCSeurat method
                default_assay, assay_names = run_seurat_qc(file)
                result.append({
                        "file": file,
                        "format": "H5Seurat or RDS",
                        "default_assay": default_assay,
                        "assay_names": assay_names
                    })
            else:
                # It's a different file, call load_annData method
                adata = LoadAnndata(file)
                adata_path = change_file_extension(file)
                adata.write_h5ad(adata_path)
                result.append({"file": file, "format": "h5ad", "result": adata_path})
        
        except Exception as e:
            # Handle the exception and return an error response
            raise HTTPException(status_code=500, detail=str(e))

    return result
