from sqlalchemy.orm import Session

from . import models, schemas

def get_spatial_datasets(db:Session):
    return db.query(models.Spatial_Dataset).all()

def get_spatial_dataset_by_id(db:Session, id:str):
    return db.query(models.Spatial_Dataset).filter(models.Spatial_Dataset.dataset_id==id).first()

def get_dataset_factors(db:Session):
    return db.query(models.Dataset_Factor).all()

def get_clustering_publications(db:Session):
    return db.query(models.Clustering_Publication).all()

def get_clustering_dataset_factor(db:Session):
    return db.query(models.Clustering_Dataset_factor).all()