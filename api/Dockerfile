# FROM u15ec104/ai-single-cell:celery-backend-base
FROM stormgod/oscb-api

# Install python packages
RUN pip install --upgrade --no-cache \
    dash~=2.9.1 pymongo dash-bootstrap-components dash_bootstrap_templates python-multipart numpyro==0.12.1
    # louvain plotly psutil nvidia-ml-py pymde
    # gputil

RUN mkdir /app
WORKDIR /app
COPY . .

# CMD ["bash"]
EXPOSE 5000
EXPOSE 5555