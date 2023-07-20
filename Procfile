web: gunicorn app:server
worker: celery -A celery_worker worker --loglevel=INFO --concurrency=2