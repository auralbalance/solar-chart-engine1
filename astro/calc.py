# astro/calc.py

import os
from datetime import datetime, time as dtime
from typing import Tuple, Dict, Any

import numpy as np
from geopy.geocoders import Nominatim
from timezonefinder import TimezoneFinder
from pytz import timezone as tzlib, utc

from skyfield.api import Loader, Topos
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import SkyCoord, GeocentricTrueEcliptic, get_constellation
from astroquery.jplhorizons import Horizons

from .signs import constellation_to_13sign, cycle_from, sign_from_ecliptic_lon


# --------------------------
# EPHEMERIS LOADER (Skyfield)
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
# GEO + TIME HELPERS
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


def to_utc(date_str: str, time_str: str | None, tz_str: str) -> Tuple[datetime, bool]:
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
# ASTRONOMY HELPERS
# --------------------------
def ecl_lon_from_skyfield(app) -> float:
    """
    Get apparent ecliptic longitude (deg 0..360) from a Skyfield apparent position.
    """
    ecl = app.ecliptic_latlon()
    return float(ecl[1].degrees % 360.0)


def compute_ascendant(ts, observer, tt, lat_deg: float, lon_deg: float) -> float:
    """
    Simple MVP ascendant longitude.
    Uses local apparent sidereal time to approximate the ecliptic longitude
    rising at the eastern horizon. (Good enough for v1.)
    """
    # Skyfield time object has gmst; convert to local sidereal time (hours)
    lst_hours = (ts.gmst + (lon_deg / 15.0)) % 24.0
    # Convert to degrees along the equator, then approximate to ecliptic lon
    # A more exact solution uses horizon intersection; for MVP this is acceptable.
    asc_approx = (lst_hours * 15.0) % 360.0
    return float(asc_approx)


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


def chiron_ecl_lon(dt_utc: datetime, lat_deg: float, lon_deg: float) -> float:
    """
    Get Chiron topocentric apparent RA/Dec from JPL Horizons, convert to ecliptic lon (deg).
    Falls back to geocentric if topocentric fails.
    """
    loc_str = f"{lon_deg},{lat_deg},0"  # lon,lat,elv(m)
    t = Time(dt_utc)
    try:
        obj = Horizons(id='2060', id_type='smallbody', location=loc_str, epochs=t.jd)
        eph = obj.ephemerides()
    except Exception:
        obj = Horizons(id='2060', id_type='smallbody', location='geo', epochs=t.jd)
        eph = obj.ephemerides()

    ra = float(eph['RA'][0])   # deg
    dec = float(eph['DEC'][0]) # deg
    sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    ecl = sc.transform_to(GeocentricTrueEcliptic(equinox=t))
    return float(ecl.lon.to(u.deg).value % 360.0)


# --------------------------
# MAIN CHART COMPUTATION
# --------------------------
def compute_chart(name: str, date: str, time_str: str | None, place: str) -> Tuple[Dict[str, Any], bool]:
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

    # 4) Ascendant (MVP) + cosmetic constellation (not used for sign)
    asc_lon = compute_ascendant(ts, observer, ts.tt, lat, lon)
    asc_sign = sign_from_ecliptic_lon(asc_lon)

    # Cosmetic (equatorial) constellation label for asc — purely informative:
    try:
        dummy = observer.at(ts).observe(eph['sun']).apparent()
        ra, dec, _ = dummy.radec()
        sc = SkyCoord(ra=ra.hours*u.hourangle, dec=dec.degrees*u.deg, frame='icrs')
        asc_const = constellation_to_13sign(get_constellation(sc))
    except Exception:
        asc_const = ""

    ascendant = {'lon': asc_lon, 'sign': asc_sign, 'constellation': asc_const}

    # 5) Planets (apparent ecliptic longitudes)
    # Use barycenters for outer planets for stability; Pluto is fine as 'pluto'
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

        # Equatorial constellation (mapped/merged) for info
        try:
            ra, dec, _ = app.radec()
            sc = SkyCoord(ra=ra.hours*u.hourangle, dec=dec.degrees*u.deg, frame='icrs')
            const_raw = get_constellation(sc)
        except Exception:
            const_raw = ""

        planets.append({
            'body': label,
            'lon': lon_ecl,
            'sign': sign_from_ecliptic_lon(lon_ecl),
            'constellation': constellation_to_13sign(const_raw)
        })

    # 6) Nodes
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = south_node_lon(nn_lon)
    nodes = {
        'north_node': {'lon': nn_lon, 'sign': sign_from_ecliptic_lon(nn_lon)},
        'south_node': {'lon': sn_lon, 'sign': sign_from_ecliptic_lon(sn_lon)}
    }

    # 7) Midheaven
    mc_lon = midheaven_lon(dt_utc, lon)
    mc = {'lon': mc_lon, 'sign': sign_from_ecliptic_lon(mc_lon)}

    # 8) Chiron
    try:
        chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
        planets.append({
            'body': 'Chiron',
            'lon': chi_lon,
            'sign': sign_from_ecliptic_lon(chi_lon),
            'constellation': ''  # optional: compute via ICRS→IAU like above if desired
        })
    except Exception:
        # If Horizons fails (network, rate-limit), skip to keep API responsive
        pass

    # 9) Houses: whole-sign starting from the Asc sign
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
