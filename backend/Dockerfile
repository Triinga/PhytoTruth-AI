# Build context is the REPOSITORY ROOT, not backend/.
# The application resolves experiments/ and artifacts/ relative to the repo
# root (prediction_service.py and main.py use Path(__file__).parents[2]), so
# the image must mirror the repository layout: backend/ lives under /app.
FROM python:3.12-slim

WORKDIR /app

COPY backend/requirements.txt ./backend/requirements.txt
RUN pip install --no-cache-dir -r backend/requirements.txt

COPY backend/app ./backend/app
COPY backend/data ./backend/data

# Executed artifacts the API reads at runtime.
COPY experiments/ml_activity_chembl335 ./experiments/ml_activity_chembl335
COPY experiments/evidence_intelligence/results.json ./experiments/evidence_intelligence/
COPY artifacts/npass_screen_100000 ./artifacts/npass_screen_100000

WORKDIR /app/backend

EXPOSE 8000
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
