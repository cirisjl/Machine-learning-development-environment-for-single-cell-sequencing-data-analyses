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
#     resolution: float = Field(default=0.5)
#     regress_cell_cycle: Optional[bool] = False
#     use_default: Optional[bool] = True



class QCParameters(BaseModel):
    methods: Optional[List[str]]= None
    assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    min_genes: int = 200
    max_genes: int = 0
    min_cells: int = 2
    target_sum: int =1e4
    n_top_genes: int = 2000
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 0.5
    doublet_rate: Optional[float] = 0
    regress_cell_cycle: Optional[bool] = False
    use_default: Optional[bool] = True 
    # Bioconductor
    colour_by: Optional[str] = 'NULL' # Color by for plots
    shape_by_1: Optional[str] = 'NULL'  # Shape by 1 for plots
    shape_by_2: Optional[str] = 'NULL'  # Shape by 2 for plots



class imputationParameters(BaseModel):
    methods: Optional[List[str]]= None
    assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    genes: Optional[List[str]] = None
    ncores: Optional[int] = 12
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 0.5



class normalizationParameters(BaseModel):
    methods: Optional[List[str]]= None
    assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 0.5
    use_default: Optional[bool] = True



class reductionParameters(BaseModel):
    assay: Optional[str] = 'RNA' # Required for Seurat
    layer: Optional[str] = None
    n_neighbors: int = 15
    n_pcs: int = 1 # Scanpy
    resolution: float = 0.5
    use_default: Optional[bool] = True



class Dataset(BaseModel):
    dataset: Optional[str] = None # Title of datasets
    input: str
    output: Optional[str] = None
    userID: Optional[str] = None
    description: Optional[str] = None,
    dataset_id: Optional[str] = None
    method: Optional[str] = None,
    process: Optional[str] = None,
    output_format: Optional[str] = 'AnnData'
    species: Optional[str] = 'human' # c("human", "mouse") Species of the database for annotation. Allowed input is human or mouse.
    idtype: Optional[str] = 'SYMBOL' # idtype should be one of "SYMBOL", "ENSEMBL", "ENTREZID" or "REFSEQ".
    cluster_label: Optional[str] = None
    qc_params: QCParameters = Field(default_factory=QCParameters)
    imputation_params: imputationParameters = Field(default_factory=imputationParameters)
    normalization_params: normalizationParameters = Field(default_factory=normalizationParameters)
    reduction_params: reductionParameters = Field(default_factory=reductionParameters)
    show_umap: Optional[bool] = True
    show_error: Optional[bool] = True   
    


class integrationParameters(BaseModel):
    # description: Optional[str] = None
    # dataset_id: Optional[str] = None
    # method: Optional[str] = None
    dims: Optional[int] = 30
    npcs: Optional[int] = 30
    default_assay: Optional[str] = 'RNA' # Required for Seurat
    reference: Optional[str] = None
    # show_umap: Optional[bool] = True
    # show_error: Optional[bool] = True



class IntegrationDataset(BaseModel):
    process: Optional[str] = None
    description: Optional[str] = None
    datasetIds: List[str] = None
    dataset: List[str] = None
    input: List[str] = None
    output: str
    userID: str
    # output_format: str
    methods: List[str] = None
    method: Optional[str] = None,
    params: integrationParameters = Field(default_factory=integrationParameters)
   


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
    benchmarksId: str
    datasetId: str
    userID: Optional[str] = None
    adata_path: str
    train_fraction: float
    validation_fraction: float
    test_fraction: float



class SubsetDataRequest(BaseModel):
    benchmarksId: str
    datasetId: str
    userID: Optional[str] = None
    adata_path: str
    obskey: str
    values: list



class TaskDataRequest(BaseModel):
    adata_path: str
    task_label: str
    datasetId: str



class BenchmarksRequest(BaseModel):
    benchmarksId: str
    datasetId: str
    userID: Optional[str] = None
    task_type: str
    adata_path: str
    label: str
    # data: List[TaskDataRequest]



class UploadRequest(BaseModel):
    fileDetails: List[str]
    assay_name: Optional[str] = None
    userID: str



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


class ProcessResultsRequest(BaseModel):
    process_ids: List[str]
    clustering_plot_type: Optional[str] = 'leiden'
    annotation: Optional[str] = None
    record_type: Optional[str] = None
    selected_cell_intersection: Optional[List[str]]= None