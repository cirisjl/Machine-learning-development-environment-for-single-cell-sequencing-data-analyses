#!/bin/bash
set -e

cd ..

echo "Creating the namespace and configs..."
kubectl apply -f k8s-manifests/00-namespace.yaml
kubectl apply -f k8s-manifests/01-configmaps.yaml
kubectl apply -f k8s-manifests/02-secrets.yaml
kubectl apply -f k8s-manifests/03-pvc.yaml

echo "Creating config maps from local init files..."
# MongoDB init script
kubectl create configmap initmongo-js --from-file=./initMongo.js -n oscb --dry-run=client -o yaml | kubectl apply -f -
# MySQL initialization schema
kubectl create configmap oscb-schema-sql --from-file=./oscb-ui/oscb-mysql/oscb_schema_latest.sql -n oscb --dry-run=client -o yaml | kubectl apply -f -

echo "Deploying databases and message brokers..."
kubectl apply -f k8s-manifests/04-mongodb.yaml
kubectl apply -f k8s-manifests/05-redis.yaml
kubectl apply -f k8s-manifests/06-rabbitmq.yaml
kubectl apply -f k8s-manifests/07-mysql.yaml

echo "Waiting for databases to initialize (sleeping for 10 seconds)..."
sleep 10

echo "Deploying application services..."
kubectl apply -f k8s-manifests/08-celery-api.yaml
kubectl apply -f k8s-manifests/09-workers.yaml
kubectl apply -f k8s-manifests/10-dashboard.yaml
kubectl apply -f k8s-manifests/11-directus.yaml
kubectl apply -f k8s-manifests/12-oscb-ui.yaml
kubectl apply -f k8s-manifests/13-jupyter.yaml

echo "Deployment submitted! Check pod status with: kubectl get pods -n oscb"
echo "To get the URLs for NodePort services (like the React UI), use: minikube service list -n oscb"
