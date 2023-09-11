FROM u15ec104/ai-single-cell:celery-backend-base

WORKDIR /usr/src/app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .

# Set the FLASK_APP environment variable
ENV FLASK_APP=main.py

# Expose port 5003 for the Flask application to listen on
EXPOSE 5003

CMD [ "python3", "-m" , "flask", "run", "--host=0.0.0.0","--port=5003"]
