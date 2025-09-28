# astro/calc.py

from __future__ import annotations
import os
from datetime import datetime, time as dtime
from typing import Tuple, Dict, Any, Optional

import numpy as np

# Geocoding & timezones
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
)

# Chiron (optional; we fail gracefully if not available)
try:
    from astroquery.jplhorizons import Horizons  # noqa: F401
    HAS_HORIZONS = True
except Exception:
    HAS_HORIZONS = False

# Your 13-sign helpers (mapping IAU -> 13 names and a small utility)
from .signs import constellation_to_13sign


# --------------------------
# Small utilities
# --------------------------
def cycle_from(start: str, n: int) -> list[str]:
    """Cycle 13-sign list starting from `start`, returning `n` items."""
    SIGNS_13 = [
        "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
        "Libra", "Scorpius", "Ophiuchus", "Sagittarius",
        "Capricornus", "Aquarius", "Pisces"
    ]
    i = SIGNS_13.index(start)
    out = []
    for k in range(n):
        out.append(SIGNS_13[(i + k) % 13])
    return out


# --------------------------
# Ephemeris loader (Skyfield)
# --------------------------
def get_ephemeris():
    """
    Load JPL DE421 ephemeris, cached under /tmp/skyfield-data on Railway.
    """
    eph_path = "/tmp/skyfield-data"
    os.makedirs(eph_path, exist_ok=True)
    load = Loader(eph_path)
    return load('de421.bsp')


# --------------------------
# Geo + time helpers
# --------------------------
def geocode_place(place_text: str) -> Tuple[float, float, str]:
    """
    Resolve a place name into (lat, lon, IANA timezone).
    Uses OpenStreetMap Nominatim + timezonefinder.
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
    if time_str is None:
        clock = dtime(12, 0)
        est = True
    else:
        hh, mm = map(int, time_str.split(":"))
        clock = dtime(hh, mm)

    dt_local = datetime.strptime(date_str, "%Y-%m-%d").replace(hour=clock.hour, minute=clock.minute)
    tz = tzlib(tz_str)
    dt_loc = tz.localize(dt_local)
    dt_utc = dt_loc.astimezone(utc)
    return dt_utc, est


# --------------------------
# Astronomy helpers
# --------------------------
def ecl_lon_from_skyfield(app) -> float:
    """
    Get apparent ecliptic longitude (deg 0..360) from a Skyfield apparent position.
    """
    ecl = app.ecliptic_latlon()
    return float(ecl[1].degrees % 360.0)


def mean_lunar_node_lon(dt_utc: datetime) -> float:
    """
    Mean lunar node longitude (Ω) in degrees.
    Simple Meeus-style approximation sufficient for labeling.
    """
    t = Time(dt_utc).jd
    T = (t - 2451545.0) / 36525.0
    Omega = 125.04452 - 1934.136261*T + 0.0020708*T*T + (T*T*T)/450000.0
    return Omega % 360.0


def south_node_lon(north_node_lon: float) -> float:
    return (north_node_lon + 180.0) % 360.0


def midheaven_lon(dt_utc: datetime, lon_deg: float) -> float:
    """
    Midheaven: ecliptic longitude of the upper culmination (local meridian).
    """
    eps = 23.43929111  # mean obliquity (deg)
    obstime = Time(dt_utc)
    # Local apparent sidereal time in degrees:
    lst_deg = obstime.sidereal_time('apparent', longitude=lon_deg*u.deg).deg
    # tan(Lmc) = tan(LST) / cos(eps)
    import math
    LST = math.radians(lst_deg)
    ce = math.cos(math.radians(eps))
    Lmc = math.degrees(math.atan2(math.tan(LST), ce)) % 360.0
    return Lmc


def chiron_ecl_lon(dt_utc: datetime, lat_deg: float, lon_deg: float) -> Optional[float]:
    """
    Get Chiron topocentric apparent RA/Dec from JPL Horizons, convert to ecliptic lon (deg).
    Falls back to geocentric if topocentric fails. Returns None if Horizons unavailable.
    """
    if not HAS_HORIZONS:
        return None

    loc_str = f"{lon_deg},{lat_deg},0"  # lon,lat,elv(m)
    t = Time(dt_utc)
    try:
        obj = Horizons(id='2060', id_type='smallbody', location=loc_str, epochs=t.jd)
        eph = obj.ephemerides()
    except Exception:
        try:
            obj = Horizons(id='2060', id_type='smallbody', location='geo', epochs=t.jd)
            eph = obj.ephemerides()
        except Exception:
            return None

    ra = float(eph['RA'][0])   # deg
    dec = float(eph['DEC'][0]) # deg
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    ecl = sc.transform_to(GeocentricTrueEcliptic(equinox=t))
    return float(ecl.lon.to(u.deg).value % 360.0)


# --------------------------
# Ascendant solver (accurate)
# --------------------------
def compute_ascendant(ts, observer, lat_deg: float, lon_deg: float) -> float:
    """
    Find the ecliptic longitude (deg) where the ecliptic (lat=0) rises at the local
    eastern horizon (alt≈0°, 0°<az<180°). Numeric search using Skyfield's alt/az.
    """

    def alt_az_of_ecl_lon(L_deg: float) -> Tuple[float, float]:
        t = Time(ts.utc_datetime())
        sc = SkyCoord(lon=L_deg*u.deg, lat=0*u.deg,
                      frame=GeocentricTrueEcliptic(equinox=t))
        icrs = sc.icrs
        star = Star(ra_hours=float(icrs.ra.to(u.hourangle).value),
                    dec_degrees=float(icrs.dec.deg))
        app = observer.at(ts).observe(star).apparent()
        alt, az, _ = app.altaz()
        return float(alt.degrees), float(az.degrees)

    # coarse bracket on EAST side
    best_L, best_abs = 0.0, 1e9
    for L in range(0, 360, 2):  # 2° steps
        alt, az = alt_az_of_ecl_lon(L)
        if 0.0 < az < 180.0:
            a = abs(alt)
            if a < best_abs:
                best_abs, best_L = a, L

    # refine
    step = 1.0
    L0 = best_L
    for _ in range(12):  # ~millidegree precision
        candidates = [L0 - step, L0, L0 + step]
        best = (1e9, L0)
        for L in candidates:
            alt, az = alt_az_of_ecl_lon(L % 360.0)
            if 0.0 < az < 180.0:
                a = abs(alt)
                if a < best[0]:
                    best = (a, L % 360.0)
        L0 = best[1]
        step *= 0.5

    return L0 % 360.0


# --------------------------
# Main chart computation
# --------------------------
def compute_chart(name: str, date: str, time_str: Optional[str], place: str) -> Tuple[Dict[str, Any], bool]:
    # 1) Geocode + time
    lat, lon, tz_str = geocode_place(place)
    dt_utc, est = to_utc(date, time_str, tz_str)

    # 2) Ephemeris + timescale
    eph = get_ephemeris()
    load = Loader("/tmp/skyfield-data")
    ts = load.timescale().from_datetime(dt_utc)

    # 3) Observer
    earth = eph['earth']
    observer = earth + Topos(latitude_degrees=lat, longitude_degrees=lon)

    # 4) Ascendant (accurate) + sign label from IAU constellation at that ecliptic point
    asc_lon = compute_ascendant(ts, observer, lat, lon)
    t_ast = Time(ts.utc_datetime())
    asc_sc = SkyCoord(lon=asc_lon*u.deg, lat=0*u.deg, frame=GeocentricTrueEcliptic(equinox=t_ast))
    asc_icrs = asc_sc.icrs
    asc_const_raw = get_constellation(asc_icrs)          # IAU name
    asc_sign = constellation_to_13sign(asc_const_raw)    # map to 13-sign label
    ascendant = {'lon': asc_lon, 'sign': asc_sign, 'constellation': asc_sign}

    # 5) Planets (apparent ecliptic longitudes; barycenters where needed)
    keys = {
        'Sun': 'sun',
        'Moon': 'moon',
        'Mercury': 'mercury',
        'Venus': 'venus',
        'Mars': 'mars',
        'Jupiter': 'jupiter barycenter',
        'Saturn': 'saturn barycenter',
        'Uranus': 'uranus barycenter',
        'Neptune': 'neptune barycenter',
        'Pluto': 'pluto barycenter'
    }

    planets = []
    for label, eph_key in keys.items():
        target = eph[eph_key]
        app = observer.at(ts).observe(target).apparent()
        lon_ecl = ecl_lon_from_skyfield(app)

        # Equatorial constellation (mapped) for sign label
        try:
            ra, dec, _ = app.radec()
            sc = SkyCoord(ra=ra.hours*u.hourangle, dec=dec.degrees*u.deg, frame='icrs')
            const_raw = get_constellation(sc)      # e.g., 'Libra', 'Scorpius', 'Ophiuchus', ...
            sign13 = constellation_to_13sign(const_raw)
        except Exception:
            sign13 = ""

        planets.append({
            'body': label,
            'lon': lon_ecl,
            'sign': sign13,
            'constellation': sign13
        })

    # 6) Nodes (mean)
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = south_node_lon(nn_lon)

    def sign_from_ecl_point(lon_deg: float) -> str:
        """Sign label for an ecliptic point: transform point -> ICRS -> IAU constellation -> 13-sign."""
        tpt = SkyCoord(lon=lon_deg*u.deg, lat=0*u.deg, frame=GeocentricTrueEcliptic(equinox=t_ast)).icrs
        return constellation_to_13sign(get_constellation(tpt))

    nodes = {
        'north_node': {'lon': nn_lon, 'sign': sign_from_ecl_point(nn_lon)},
        'south_node': {'lon': sn_lon, 'sign': sign_from_ecl_point(sn_lon)}
    }

    # 7) Midheaven
    mc_lon = midheaven_lon(dt_utc, lon)
    mc = {'lon': mc_lon, 'sign': sign_from_ecl_point(mc_lon)}

    # 8) Chiron (optional)
    chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
    if chi_lon is not None:
        planets.append({
            'body': 'Chiron',
            'lon': chi_lon,
            'sign': sign_from_ecl_point(chi_lon),
            'constellation': sign_from_ecl_point(chi_lon)
        })

    # 9) Houses: whole-sign from the ASC sign
    houses = [{'house': i+1, 'sign': s} for i, s in enumerate(cycle_from(asc_sign, 12))]

    # 10) Payload
    payload = {
        'name': name,
        'datetime_utc': dt_utc.strftime('%Y-%m-%dT%H:%M:%SZ'),
        'location': {'lat': lat, 'lon': lon, 'tz': tz_str},
        'ascendant': ascendant,
        'midheaven': mc,
        'nodes': nodes,
        'planets': planets,
        'houses': houses
    }

    return payload, est
