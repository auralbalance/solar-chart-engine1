# app/main.py
import os, uuid
from fastapi import FastAPI, HTTPException
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from astro.calc import compute_chart
from astro.wheel import make_svg_wheel

# --- ensure runtime folders exist BEFORE mounting static ---
os.makedirs("static", exist_ok=True)
os.makedirs("skyfield_data", exist_ok=True)

class ChartRequest(BaseModel):
    name: str
    date: str            # YYYY-MM-DD
    time: str | None = None  # HH:MM (24h) or None
    place: str

app = FastAPI()
app.mount("/static", StaticFiles(directory="static"), name="static")

@app.post("/chart")
async def chart_endpoint(req: ChartRequest):
    try:
        chart_data, time_estimated = compute_chart(req.name, req.date, req.time, req.place)
        wheel_id = uuid.uuid4().hex[:10]
        svg_path = os.path.join("static", f"{wheel_id}.svg")
        make_svg_wheel(chart_data, svg_path)
        resp = dict(chart_data)
        resp["wheel_svg_url"] = f"/static/{wheel_id}.svg"
        if time_estimated:
            resp["time_estimated"] = True
        return JSONResponse(content=resp)
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.get("/")
def root():
    return {"ok": True, "try": "/docs"}

