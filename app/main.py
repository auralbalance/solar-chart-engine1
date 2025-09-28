# app/main.py
from fastapi import FastAPI, Query
from typing import Optional
from astro.calc import compute_chart, get_ephemeris
import os

app = FastAPI(title="Solar Chart API", version="1.0.0")

@app.get("/health")
def health():
    # Touch ephemeris once, report the active name
    _, _, eph_name = get_ephemeris()
    return {
        "status": "ok",
        "ephemeris": eph_name,
        "house_mode_default": os.getenv("HOUSE_MODE", "13"),
        "zodiac_profile_default": os.getenv("ZODIAC_PROFILE", "CUSTOM_13"),
    }

@app.get("/chart")
def chart(
    name: str = Query(..., description="Display name"),
    date: str = Query(..., description="YYYY-MM-DD"),
    time: Optional[str] = Query(None, description="HH:MM (24h)"),
    place: str = Query(..., description="City, Country"),
    house_mode: Optional[str] = Query(None, description='"12" or "13"'),
    zodiac_profile: Optional[str] = Query(None, description='"IAU_13" or "CUSTOM_13"'),
):
    payload, _ = compute_chart(
        name=name,
        date=date,
        time_str=time,
        place=place,
        house_mode=house_mode,
        zodiac_profile=zodiac_profile,
    )
    return payload
