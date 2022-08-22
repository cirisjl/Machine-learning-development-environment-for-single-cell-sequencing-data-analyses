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