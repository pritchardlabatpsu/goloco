import os
from celery import Celery
from dash import DiskcacheManager, CeleryManager
import logging
import app

if 'REDIS_URL' in os.environ:
    redis_url = os.environ['REDIS_URL']
else:
    redis_url = 'redis://localhost:6379/0'

celery_app = Celery(__name__, 
                    BROKER_URL=redis_url, 
                    CELERY_RESULT_BACKEND=redis_url, 
                    BROKER_POOL_LIMIT=0, 
                    include=['app']
                    )