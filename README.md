# Solar Chart Engine (True-Sky 13-Sign)

FastAPI + Skyfield + Astropy. Deployed on Railway.

## Endpoints
- `POST /chart` — returns placements JSON and a wheel SVG URL
- `GET /docs` — interactive API docs

## Run locally
```
pip install -r requirements.txt
uvicorn app.main:app --reload
```
