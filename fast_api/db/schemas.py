from datetime import date
from pydantic import BaseModel

class Spatial_Datasets(BaseModel):
    dataset_id:str
    publication_id:str
    species:str
    tissue:str
    organ_part:str
    n_spot:str

    class Config:
        orm_mode = True

class Dataset_Factor(BaseModel):
    id:int
    publication_id:str
    dataset_id:str
    factor_key:str
    factor_value:str

    class Config:
        orm_mode = True

class Clustering_Publication(BaseModel):
    publication_id:int
    species:str	
    pmid:int
    library_prep_protocol:str	
    title:str	
    journal:str	
    authors:str	
    citation:str	
    abstract:str	
    doi:str	
    cells:int	
    publish_date:date	
    load_date:date	
    n_samples:int	

    class Config:
        orm_mode = True	

class Clustering_Dataset_Factor(BaseModel):
    id	:int
    publication_id:int
    dataset_id:int
    factor_key:str
    factor_value:str
    unit	:str		

    class Config:
        orm_mode = True
