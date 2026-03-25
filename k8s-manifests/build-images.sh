#!/bin/bash
set -e

# Point terminal to Minikube's internal Docker daemon
echo "Configuring Docker environment for Minikube..."
eval $(minikube docker-env)

echo "Building local images directly inside Minikube..."

# Go back to the project root directory where docker-compose context is defined
cd ..

echo "Building Celery API / Workers (oscb-api:latest)..."
docker build -t oscb-api:latest ./api/

# Note: for MySQL, kubernetes is using the standard mysql:8.0 image with a config map for initialization, 
# but if a local build is absolutely essential we can build it. The original compose file did build it.
echo "Building MySQL (oscb-mysql:latest)..."
docker build -t oscb-mysql:latest ./oscb-ui/oscb-mysql/

echo "Building Node Backend (oscb-node:latest)..."
docker build -t oscb-node:latest ./oscb-ui/oscb-node/

echo "Building React Frontend (oscb-react:latest)..."
docker build -t oscb-react:latest ./oscb-ui/oscb-react/

echo "All images built successfully in Minikube's Docker registry."
