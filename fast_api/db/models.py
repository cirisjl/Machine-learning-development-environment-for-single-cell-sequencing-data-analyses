from sqlalchemy import Boolean, Column, ForeignKey, Integer, String
from sqlalchemy.orm import relationship

from .database import Base

class Spatial_Dataset(Base):
    __tablename__ = "spatial_dataset"

    dataset_id=Column(String,primary_key=True,index=True)
    publication_id=Column(String)
    species=Column(String)
    tissue=Column(String)
    organ_part=Column(String)
    n_spot=Column(String)

class Dataset_Factor(Base):
    __tablename__ = "dataset_factor"

    id=Column(Integer,primary_key=True,index=True)
    publication_id=Column(String)
    dataset_id=Column(String)
    factor_key=Column(String)
    factor_value=Column(String)