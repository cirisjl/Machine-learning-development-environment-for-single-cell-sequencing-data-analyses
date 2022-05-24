# 
FROM python:3.9

# 
WORKDIR /code

# 
COPY ./requirements.txt /code/requirements.txt

# 
RUN pip install --no-cache-dir --upgrade -r /code/requirements.txt

# Set environment variables
ENV place=docker

# 
COPY . /code/

EXPOSE 5000
# 
CMD ["gunicorn", "main:app", "--bind", "0.0.0.0:5000",\
    "--worker-class", "uvicorn.workers.UvicornWorker"]
