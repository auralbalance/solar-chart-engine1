# app/main.py
from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional
from astro.calc import compute_chart, get_ephemeris, get_ecliptic_bands, sign_from_ecliptic_lon
import os

ALLOWED_ORIGINS = os.getenv("ALLOWED_ORIGINS", "*").split(",")

app = FastAPI(title="Solar Chart API — True Sky", version="1.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=[o.strip() for o in ALLOWED_ORIGINS],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

METHOD_TEXT = (
    "We compute true-sky placements using IAU constellation boundaries. "
    "Positions come from JPL ephemerides via Skyfield. "
    "We transform each body to the ecliptic J2000 frame and classify its ecliptic longitude along the ecliptic line. "
    "No tropical signs. No custom bands. "
    "Houses use a 13-house whole-sign model. House 1 starts at the Ascendant’s sign. One house per constellation. "
    "The Midheaven floats. Chiron uses JPL Horizons. "
    "Time and place resolve by IANA timezone and OpenStreetMap geocoding with caching."
)

@app.get("/health")
def health():
    _, _, eph_name = get_ephemeris()
    return {"status": "ok", "ephemeris": eph_name, "houses": 13, "profile": "IAU true-sky"}

@app.get("/method")
def method():
    return {"method": METHOD_TEXT}

@app.get("/chart")
def chart(
    name: str = Query(..., description="Display name"),
    date: str = Query(..., description="YYYY-MM-DD"),
    time: Optional[str] = Query(None, description="HH:MM"),
    place: str = Query(..., description="City, Country"),
):
    payload, _ = compute_chart(name=name, date=date, time_str=time, place=place)
    return payload

@app.get("/debug/bands")
def debug_bands():
    return {"bands": get_ecliptic_bands()}

@app.get("/debug/locate")
def debug_locate(lon: float = Query(..., ge=0.0, le=360.0)):
    bands = get_ecliptic_bands()
    sign = sign_from_ecliptic_lon(lon)
    prev_band = None
    next_band = None
    current = None
    for i, b in enumerate(bands):
        s, e = b["start"], b["end"]
        if (s <= e and s <= lon < e) or (s > e and (lon >= s or lon < e)):
            current = b
            prev_band = bands[i - 1] if i > 0 else bands[-1]
            next_band = bands[(i + 1) % len(bands)]
            break
    return {"lon": lon, "sign": sign, "current_band": current, "prev_band": prev_band, "next_band": next_band}
