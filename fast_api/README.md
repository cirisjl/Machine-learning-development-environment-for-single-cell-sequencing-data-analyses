# This README shows the proccess on FastAPI
## create a screen session and start the server
```bash
screen -S scree_name
``` 
## Installation

create a virtual env and activate it. 

```bash
cd fast-api
python3 -m venv env
source env/bin/activate
```

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install fledge-db.

```bash
pip3 install -r requirements.txt
```


Create a `config.yaml` file in the root of the project with the following contents:

```yaml
database:
  host: host
  db: db
  password: password
  port: port
  user: user

admin:
  user: admin
  hashed_password: hashed_password

mode: mode
```

Generate shh keys: 
```bash
cd auth
ssh-keygen -t rsa -b 4096 -E SHA512 -m PEM -P "" -f RS512.key
```


add path to the env file to the PYTHONPATH environment variable.
```bash
export PYTHONPATH=/<path_to>/fast-api 
```
## run the server manually

```bash
gunicorn main:app --workers 4 --worker-class uvicorn.workers.UvicornWorker --bind 0.0.0.0:5000
```

## run the server in a docker container

```bash
docker-compose up --build
```

## just run docker
```bash
docker build -t dev-fast-api-bio -f Dockerfile.dev .
docker run -p 5000:5000 dev-fast-api-bio
```


## License
[MIT](https://choosealicense.com/licenses/mit/)