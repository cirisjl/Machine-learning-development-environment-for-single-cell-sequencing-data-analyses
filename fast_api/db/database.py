from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from auth.credentials import db_credentials

user=db_credentials["user"]
database=db_credentials["db"]
port=db_credentials["port"]
password=db_credentials["password"]
host=db_credentials["host"]

SQLALCHEMY_DATABASE_URL = "mysql+pymysql://{}:{}@{}:{}/{}".format(user,password,host,port,database)

engine = create_engine(SQLALCHEMY_DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
