from sqlite3 import Date
from sqlalchemy import Boolean, Column, ForeignKey, Integer, String, DateTime
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

class Clustering_Publication(Base):
    __tablename__ = "clustering_publication"

    publication_id=Column(Integer)
    species=Column(String)	
    pmid=Column(Integer,primary_key=True,index=True)
    library_prep_protocol=Column(String)
    title=Column(String)
    journal=Column(String)
    authors=Column(String)
    citation=Column(String)
    abstract=Column(String)
    doi=Column(String)
    cells=Column(Integer)
    publish_date=Column(DateTime)
    load_date=Column(DateTime)
    n_samples=Column(Integer)		

class Clustering_Dataset_factor(Base):
    __tablename__ = "clustering_dataset_factor"
    id=Column(Integer)	
    publication_id=Column(Integer,primary_key=True,index=True)
    dataset_id=Column(Integer)
    factor_key=Column(String,primary_key=True,index=True)
    factor_value=Column(String)
    unit=Column(String)			

