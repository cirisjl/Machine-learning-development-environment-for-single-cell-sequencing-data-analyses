from typing import List, Optional
from pydantic import BaseModel



class Dataset(BaseModel):
    dataset: str
    input: str
    output: str
    userID: str
    output_format: str
    methods: List[str] = None
    default_assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    path_of_scrublet_calls: Optional[str] = './scrublet_calls.tsv'
    species: Optional[str] = None
    idtype: Optional[str] = None 
    genes: Optional[List[str]] = None
    ncores: Optional[int] = 12
    colour_by: Optional[str] = None
    shape_by_1: Optional[str] = None
    shape_by_2: Optional[str] = None
    show_umap: Optional[bool] = True
    show_error: Optional[bool] = True
    

class IntegrationDataset(BaseModel):
    dataset: List[str] = None
    input: List[str] = None
    output: str
    userID: str
    output_format: str
    methods: List[str] = None
    default_assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    path_of_scrublet_calls: Optional[str] = './scrublet_calls.tsv'
    species: Optional[str] = None
    idtype: Optional[str] = None 
    genes: Optional[List[str]] = None
    reference: Optional[int] = 12
    colour_by: Optional[str] = None
    shape_by_1: Optional[str] = None
    shape_by_2: Optional[str] = None
    show_umap: Optional[bool] = True
    show_error: Optional[bool] = True
   
class PathRequest(BaseModel):
    path: str

# Define data models using Pydantic for request and response bodies
class ConversionRequest(BaseModel):
    path: str

class ConversionResponse(BaseModel):
    assay_names: list
    adata_path: str
    message: str

class InputFile(BaseModel):
    fileDetails: str
    assay: Optional[str] = None

class InputFilesRequest(BaseModel):
    inputFiles: List[InputFile]

class AnndataMetadata(BaseModel):
    layers: list
    cell_metadata: dict
    gene_metadata: dict
    nCells: int
    nGenes: int
    genes: list
    cells: list
    embeddings: list

class CombinedQCResult(BaseModel):
    scanpy_results: AnndataMetadata  # Assuming you have the AnndataMetadata model defined
    # dropkick_results: AnndataMetadata

class DataSplitRequest(BaseModel):
    data: str
    train_fraction: float
    validation_fraction: float
    test_fraction: float

class TaskDataRequest(BaseModel):
    adata_path: str
    task_label: str
    datasetId: str

class BenchmarksRequest(BaseModel):
    task_type: str
    data: List[TaskDataRequest]

class ConvertRequest(BaseModel):
    fileDetails: List[str]
    assay_name: Optional[str] = None