from typing import List, Optional
from pydantic import BaseModel, Field



# class QCParameters(BaseModel):
#     assay: Optional[str] =  Field(default='RNA')
#     doublet_rate: float = Field(default=0)
#     min_genes: int = Field(default=200)
#     max_genes: int = Field(default=0)
#     min_cells: int = Field(default=2)
#     target_sum: int = Field(default=1e4)
#     n_top_genes: int = Field(default=2000)
#     n_neighbors: int = Field(default=15)
#     n_pcs: int = Field(default=1)
#     resolution: float = Field(default=1)
#     regress_cell_cycle: Optional[bool] = False
#     use_default: Optional[bool] = True



class QCParameters(BaseModel):
    min_genes: int = 200
    max_genes: int = 0
    min_cells: int = 2
    target_sum: int =1e4
    n_top_genes: int = 2000
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 1
    doublet_rate: Optional[float] = 0
    regress_cell_cycle: Optional[bool] = False
    use_default: Optional[bool] = True
    clustering_only: Optional[bool] = False
    dimension_reduction_only: Optional[bool] = False

 
    # Bioconductor
    colour_by: Optional[str] = 'NULL' # Color by for plots
    shape_by_1: Optional[str] = 'NULL'  # Shape by 1 for plots
    shape_by_2: Optional[str] = 'NULL'  # Shape by 2 for plots



class imputationParameters(BaseModel):
    genes: Optional[List[str]] = None
    ncores: Optional[int] = 12



class normalizationParameters(BaseModel):
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 1
    use_default: Optional[bool] = True



class Dataset(BaseModel):
    dataset: Optional[str] = None # Tittle of datasets
    input: str
    output: Optional[str]= None
    userID: Optional[str]= None
    task_id: Optional[str]= None
    output_format: Optional[str] = 'AnnData'
    methods: Optional[List[str]]= None
    assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str]= None
    species: Optional[str]= 'human' # c("human", "mouse") Species of the database for annotation. Allowed input is human or mouse.
    idtype: Optional[str]= 'SYMBOL' # idtype should be one of "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
    qc_params: QCParameters = Field(default_factory=QCParameters)
    imputation_params: imputationParameters = Field(default_factory=imputationParameters)
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



class UMAPRequest(BaseModel):
    adata_path: str
    layer: str
    clustering_plot_type: str
    selected_cell_intersection: list
    n_dim: int



class CombinedQCResult(BaseModel):
    scanpy_results: AnndataMetadata  # Assuming you have the AnndataMetadata model defined
    # dropkick_results: AnndataMetadata



class DataSplitRequest(BaseModel):
    data: str
    train_fraction: float
    validation_fraction: float
    test_fraction: float



class SubsetDataRequest(BaseModel):
    data: str
    obskey: str
    values: list



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



class QualityControlRequest(BaseModel):
    fileDetails: str
    assay: Optional[str] = None
    min_genes: int
    max_genes: int
    min_cells: int
    target_sum: float
    n_top_genes: int
    n_neighbors: int
    n_pcs: int
    resolution: float
    regress_cell_cycle: bool
    use_default: bool
    doublet_rate: float
    unique_id: str