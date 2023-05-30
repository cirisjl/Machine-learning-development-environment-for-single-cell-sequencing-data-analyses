from typing import List, Optional
from pydantic import BaseModel



class Dataset(BaseModel):
    dataset: str
    input: str
    output: str
    output_format: str
    methods: List[str] = None
    default_assay: Optional[str] = 'RNA'
    layer: Optional[str] = None
    path_of_scrublet_calls: Optional[str] = './scrublet_calls.tsv'
    species: Optional[str] = None
    idtype: Optional[str] = None 
    genes: Optional[List[str]] = None
    ncores: Optional[int] = 12
    show_umap: Optional[bool] = True
    show_error: Optional[bool] = True
