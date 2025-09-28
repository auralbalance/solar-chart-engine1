# astro/calc.py

from __future__ import annotations
import os
from datetime import datetime, time as dtime
from typing import Tuple, Dict, Any, Optional, List

import numpy as np

# Geo + TZ
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder
from pytz import timezone as tzlib, utc

# Astronomy
from skyfield.api import Loader, Topos, Star
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import (
    SkyCoord,
    GeocentricTrueEcliptic,
    get_constellation,
    ICRS,
)

# Optional: Chiron via JPL Horizons. We fail gracefully if not available.
try:
    from astroquery.jplhorizons import Horizons  # type: ignore
    HAS_HORIZONS = True
except Exception:
    HAS_HORIZONS = False


# ---------------------------------------------------------------------
# EPHEMERIS LOADER (use de440s: includes Pluto center and modern accuracy)
# ---------------------------------------------------------------------
def get_ephemeris():
    """
    Load JPL DE440s ephemeris, cached under /tmp/skyfield-data (Railway-friendly).
    """
    data_dir = "/tmp/skyfield-data"
    os.makedirs(data_dir, exist_ok=True)
    load = Loader(data_dir)
    return load("de440s.bsp"), load


# ---------------------------------------------------------------------
# GEO / TIME HELPERS
# ---------------------------------------------------------------------
def geocode_place(place_text: str) -> Tuple[float, float, str]:
    """
    Resolve 'place' to (lat, lon, tz). Uses OpenStreetMap Nominatim + TimezoneFinder.
    """
    geocoder = Nominatim(user_agent="solar-chart-api", timeout=15)
    loc = geocoder.geocode(place_text, addressdetails=False, language="en")
    if not loc:
        raise ValueError(f"Place '{place_text}' not found.")
    lat = float(loc.latitude)
    lon = float(loc.longitude)

    tf = TimezoneFinder()
    tz = tf.timezone_at(lng=lon, lat=lat)
    if not tz:
        raise ValueError("Timezone not found for the given location.")
    return lat, lon, tz


def to_utc(date_str: str, time_str: Optional[str], tz_str: str) -> Tuple[datetime, bool]:
    """
    Convert local date/time + IANA tz → timezone-aware UTC datetime.
    If time is None, default to 12:00 and set estimated flag True.
    """
    est = False
    if not time_str:
        clock = dtime(12, 0)
        est = True
    else:
        hh, mm = map(int, time_str.split(":"))
        clock = dtime(hh, mm)
    dt_local = datetime.strptime(date_str, "%Y-%m-%d").replace(
        hour=clock.hour, minute=clock.minute
    )
    tz = tzlib(tz_str)
    dt_loc = tz.localize(dt_local)
    return dt_loc.astimezone(utc), est


# ---------------------------------------------------------------------
# 13-SIGN MAPPING
# We build ecliptic bands by sampling β=0° at J2000, mapping IAU constellations
# to 13-sign labels. Everything else (planets, ASC, MC, Nodes, Chiron) is
# computed and converted to ecliptic J2000 to match these bands exactly.
# ---------------------------------------------------------------------

# IAU → 13-sign label map
# Libra→Libra, Scorpius→Scorpius, Ophiuchus included, Capricornus, etc.
CONSTELLATION_TO_13 = {
    "Aries": "Aries",
    "Taurus": "Taurus",
    "Gemini": "Gemini",
    "Cancer": "Cancer",
    "Leo": "Leo",
    "Virgo": "Virgo",
    "Libra": "Libra",
    "Scorpius": "Scorpius",
    "Ophiuchus": "Ophiuchus",
    "Sagittarius": "Sagittarius",
    "Capricornus": "Capricornus",
    "Aquarius": "Aquarius",
    "Pisces": "Pisces",
}

def constellation_to_13sign(iau_name: str) -> str:
    return CONSTELLATION_TO_13.get(iau_name, iau_name)

_ECLIPTIC_BANDS: List[Tuple[float, float, str]] = []  # [(start_lon, end_lon, label13)]


def _build_ecliptic_bands(step_deg: float = 0.2) -> List[Tuple[float, float, str]]:
    """
    Sample β=0° points around the ecliptic at J2000 and record where the IAU constellation
    (mapped to 13-sign) changes. Returns bands covering [0,360). Frame: ecliptic J2000.
    """
    bands: List[Tuple[float, float, str]] = []
    t = Time("J2000")
    current_label: Optional[str] = None
    seg_start = 0.0

    n = int(360.0 / step_deg)
    for i in range(n + 1):  # include 360
        L = (i * step_deg) % 360.0
        sc = SkyCoord(
            lon=L * u.deg,
            lat=0.0 * u.deg,
            frame=GeocentricTrueEcliptic(equinox=t),
        ).icrs
        label = constellation_to_13sign(get_constellation(sc))
        if current_label is None:
            current_label = label
            seg_start = L
        elif label != current_label:
            bands.append((seg_start, L, current_label))
            seg_start = L
            current_label = label

    if current_label is not None:
        bands.append((seg_start, 360.0, current_label))

    # Merge wrap if first and last share label
    if bands and bands[0][2] == bands[-1][2]:
        first_s, first_e, lab = bands[0]
        last_s, last_e, _ = bands[-1]
        bands = [(last_s - 360.0, first_e, lab)] + bands[1:-1]

    bands.sort(key=lambda x: x[0])
    return bands


def _ensure_bands_built():
    global _ECLIPTIC_BANDS
    if not _ECLIPTIC_BANDS:
        _ECLIPTIC_BANDS = _build_ecliptic_bands(step_deg=0.2)


def sign_from_ecliptic_lon(lon_deg: float) -> str:
    """Return 13-sign label for an ecliptic J2000 longitude using the precomputed bands."""
    _ensure_bands_built()
    x = lon_deg % 360.0
    for s, e, lab in _ECLIPTIC_BANDS:
        # bands may include a negative start for the wrap segment
        if s < 0:
            if (0 <= x < e) or (x + 360.0 >= s and x + 360.0 < e):
                return lab
        else:
            if s <= x < e:
                return lab
    return "Pisces"  # should not happen


# ---------------------------------------------------------------------
# CORE ASTRONOMY HELPERS — all converted to ECLIPTIC J2000
# ---------------------------------------------------------------------
def ecl_lon_from_app(app) -> float:
    """
    Apparent RA/Dec from Skyfield → ICRS → ecliptic J2000 → longitude in degrees [0..360).
    Keeps the same frame as the precomputed IAU bands (J2000).
    """
    ra, dec, _ = app.radec()  # hours, degrees
    sc_icrs = SkyCoord(ra=float(ra.hours) * u.hourangle,
                       dec=float(dec.degrees) * u.deg,
                       frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


def mean_lunar_node_lon(dt_utc: datetime) -> float:
    """
    Mean lunar node longitude Ω (deg). Meeus-style approximation (good for labeling).
    """
    t = Time(dt_utc).jd
    T = (t - 2451545.0) / 36525.0
    Omega = 125.04452 - 1934.136261 * T + 0.0020708 * T * T + (T ** 3) / 450000.0
    return Omega % 360.0


def midheaven_lon(dt_utc: datetime, lon_deg: float) -> float:
    """
    MC: take local apparent sidereal time as RA on the meridian, convert to ecliptic J2000.
    """
    obstime = Time(dt_utc)
    lst_deg = obstime.sidereal_time("apparent", longitude=lon_deg * u.deg).deg
    sc_icrs = SkyCoord(ra=lst_deg * u.deg, dec=0.0 * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


def chiron_ecl_lon(dt_utc: datetime, lat_deg: float, lon_deg: float) -> Optional[float]:
    """
    Chiron ecliptic longitude via JPL Horizons (topocentric if possible; else geocentric).
    Returns None if Horizons unavailable/errors. Output is ecliptic J2000 longitude.
    """
    if not HAS_HORIZONS:
        return None
    t = Time(dt_utc)
    loc_str = f"{lon_deg},{lat_deg},0"
    try:
        obj = Horizons(id="2060", id_type="smallbody", location=loc_str, epochs=t.jd)
        eph = obj.ephemerides()
    except Exception:
        try:
            obj = Horizons(id="2060", id_type="smallbody", location="geo", epochs=t.jd)
            eph = obj.ephemerides()
        except Exception:
            return None
    ra = float(eph["RA"][0])
    dec = float(eph["DEC"][0])
    sc_icrs = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


# ---------------------------------------------------------------------
# ASCENDANT SOLVER (accurate) — returns ecliptic J2000 longitude
# ---------------------------------------------------------------------
def compute_ascendant(ts, observer) -> float:
    """
    Find ecliptic longitude (deg) of the rising point (β=0°, alt≈0°) on the EAST horizon.
    Numerical search. All in ecliptic J2000 to match band mapping.
    """
    def alt_az_of_ecl_lon_j2000(L_deg: float):
        sc_ecl = SkyCoord(lon=L_deg * u.deg, lat=0 * u.deg,
                          frame=GeocentricTrueEcliptic(equinox=Time("J2000")))
        icrs = sc_ecl.icrs
        star = Star(ra_hours=float(icrs.ra.to(u.hourangle).value),
                    dec_degrees=float(icrs.dec.deg))
        app = observer.at(ts).observe(star).apparent()
        alt, az, _ = app.altaz()
        return float(alt.degrees), float(az.degrees)

    # coarse bracket on EAST side
    best_L, best_abs = 0.0, 1e9
    for L in range(0, 360, 2):  # 2° step
        alt, az = alt_az_of_ecl_lon_j2000(L)
        if 0.0 < az < 180.0:
            a = abs(alt)
            if a < best_abs:
                best_abs, best_L = a, L

    # refine
    step = 1.0
    L0 = best_L
    for _ in range(12):  # ~millidegree
        candidates = [L0 - step, L0, L0 + step]
        best = (1e9, L0)
        for L in candidates:
            alt, az = alt_az_of_ecl_lon_j2000(L % 360.0)
            if 0.0 < az < 180.0:
                a = abs(alt)
                if a < best[0]:
                    best = (a, L % 360.0)
        L0 = best[1]
        step *= 0.5

    return L0 % 360.0


# ---------------------------------------------------------------------
# MAIN CHART
# ---------------------------------------------------------------------
def compute_chart(name: str, date: str, time_str: Optional[str], place: str) -> Tuple[Dict[str, Any], bool]:
    # 1) Geo + time
    lat, lon, tz_str = geocode_place(place)
    dt_utc, time_estimated = to_utc(date, time_str, tz_str)

    # 2) Ephemeris + timescale + observer
    eph, load = get_ephemeris()
    ts = load.timescale().from_datetime(dt_utc)
    earth = eph["earth"]
    observer = earth + Topos(latitude_degrees=lat, longitude_degrees=lon)

    # 3) Ascendant (solve) + sign via ecliptic bands
    asc_lon = compute_ascendant(ts, observer)
    asc_sign = sign_from_ecliptic_lon(asc_lon)
    ascendant = {"lon": asc_lon, "sign": asc_sign, "constellation": asc_sign}

    # 4) Planets (apparent; convert to ecliptic J2000; use Pluto center)
    targets = {
        "Sun": "sun",
        "Moon": "moon",
        "Mercury": "mercury",
        "Venus": "venus",
        "Mars": "mars",
        "Jupiter": "jupiter barycenter",
        "Saturn": "saturn barycenter",
        "Uranus": "uranus barycenter",
        "Neptune": "neptune barycenter",
        "Pluto": "pluto",  # needs de440s.bsp
    }

    planets: List[Dict[str, Any]] = []
    for label, key in targets.items():
        app = observer.at(ts).observe(eph[key]).apparent()
        L = ecl_lon_from_app(app)
        sign13 = sign_from_ecliptic_lon(L)
        planets.append(
            {"body": label, "lon": L, "sign": sign13, "constellation": sign13}
        )

    # 5) Nodes (mean) — convert Ω to sign using same mapper
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = (nn_lon + 180.0) % 360.0
    nodes = {
        "north_node": {"lon": nn_lon, "sign": sign_from_ecliptic_lon(nn_lon)},
        "south_node": {"lon": sn_lon, "sign": sign_from_ecliptic_lon(sn_lon)},
    }

    # 6) Midheaven
    mc_lon = midheaven_lon(dt_utc, lon)
    mc = {"lon": mc_lon, "sign": sign_from_ecliptic_lon(mc_lon)}

    # 7) Chiron (optional)
    chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
    if chi_lon is not None:
        planets.append(
            {
                "body": "Chiron",
                "lon": chi_lon,
                "sign": sign_from_ecliptic_lon(chi_lon),
                "constellation": sign_from_ecliptic_lon(chi_lon),
            }
        )

    # 8) Houses: whole-sign starting from ASC sign.
    # You are outputting 12 houses on a 13-sign belt. That’s fine, but be explicit in UI.
    SIGNS_13 = [
        "Aries",
        "Taurus",
        "Gemini",
        "Cancer",
        "Leo",
        "Virgo",
        "Libra",
        "Scorpius",
        "Ophiuchus",
        "Sagittarius",
        "Capricornus",
        "Aquarius",
        "Pisces",
    ]
    start_idx = SIGNS_13.index(asc_sign)
    houses = [{"house": i + 1, "sign": SIGNS_13[(start_idx + i) % 13]} for i in range(12)]

    # 9) Payload
    payload = {
        "name": name,
        "datetime_utc": dt_utc.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "location": {"lat": lat, "lon": lon, "tz": tz_str},
        "ascendant": ascendant,
        "midheaven": mc,
        "nodes": nodes,
        "planets": planets,
        "houses": houses,
        "frames": {"ecliptic_bands": "J2000", "positions": "J2000"},
        "ephemeris": "de440s.bsp",
        "house_model": "Whole-sign start from ASC, 12 houses on 13-sign belt",
        "time_estimated": time_estimated,
    }
    return payload, time_estimated
