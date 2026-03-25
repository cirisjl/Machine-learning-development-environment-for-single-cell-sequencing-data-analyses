#!/bin/bash
set -e

# Download Minikube
echo "Downloading Minikube..."
curl -LO https://storage.googleapis.com/minikube/releases/latest/minikube-linux-amd64
sudo install minikube-linux-amd64 /usr/local/bin/minikube
rm minikube-linux-amd64

# Download kubectl
echo "Downloading kubectl..."
curl -LO "https://dl.k8s.io/release/$(curl -L -s https://dl.k8s.io/release/stable.txt)/bin/linux/amd64/kubectl"
sudo install -o root -g root -m 0755 kubectl /usr/local/bin/kubectl
rm kubectl

# Start Minikube with Docker driver
echo "Starting Minikube..."
# Adjust CPU/Memory defaults depending on your machine
minikube start --driver=docker --cpus="max" --memory="max"

echo "Minikube is running!"
minikube status
