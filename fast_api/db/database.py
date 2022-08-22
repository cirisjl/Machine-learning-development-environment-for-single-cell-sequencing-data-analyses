from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

user="root"
password="Akhil721"
host="localhost"
port=3306
database="aiready"

SQLALCHEMY_DATABASE_URL = "mysql+pymysql://{}:{}@{}:{}/{}".format(user,password,host,port,database)

engine = create_engine(SQLALCHEMY_DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()
