web: gunicorn app:server
worker: celery -A app.celery_app worker --loglevel=INFO --concurrency=2