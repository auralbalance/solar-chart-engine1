# app/main.py
from fastapi import FastAPI, Query
from typing import Optional
from astro.calc import compute_chart, get_ephemeris
import os

app = FastAPI(title="Solar Chart API — True Sky", version="1.0.0")

@app.get("/health")
def health():
    _, _, eph_name = get_ephemeris()
    return {"status": "ok", "ephemeris": eph_name, "houses": 13, "profile": "IAU true-sky"}

@app.get("/chart")
def chart(
    name: str = Query(..., description="Display name"),
    date: str = Query(..., description="YYYY-MM-DD"),
    time: Optional[str] = Query(None, description="HH:MM (24h)"),
    place: str = Query(..., description="City, Country"),
):
    payload, _ = compute_chart(name=name, date=date, time_str=time, place=place)
    return payload
