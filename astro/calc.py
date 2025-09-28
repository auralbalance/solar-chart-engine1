# astro/calc.py

from __future__ import annotations
import os
from datetime import datetime, time as dtime
from typing import Tuple, Dict, Any, Optional, List

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

# Optional: Chiron via JPL Horizons
try:
    from astroquery.jplhorizons import Horizons  # type: ignore
    HAS_HORIZONS = True
except Exception:
    HAS_HORIZONS = False


# =========================
# EPHEMERIS + TARGET HELPERS
# =========================
def get_ephemeris():
    """
    Load a robust ephemeris, cached to /tmp/skyfield-data (Railway-friendly).
    Preference: de441.bsp → de440.bsp → de421.bsp
    """
    data_dir = "/tmp/skyfield-data"
    os.makedirs(data_dir, exist_ok=True)
    load = Loader(data_dir)
    for fname in ("de441.bsp", "de440.bsp", "de421.bsp"):
        try:
            eph = load(fname)
            return eph, load, fname
        except Exception:
            continue
    # Last resort: raise a clear error
    raise RuntimeError("No ephemeris available (tried de441.bsp, de440.bsp, de421.bsp).")


def pick_key(eph, preferred: str, fallback: str) -> str:
    """
    Return preferred target key if present in eph, else fallback if present,
    else raise helpful error.
    """
    try:
        _ = eph[preferred]
        return preferred
    except Exception:
        try:
            _ = eph[fallback]
            return fallback
        except Exception:
            raise KeyError(
                f"Neither '{preferred}' nor '{fallback}' available in ephemeris."
            )


# =================
# GEO / TIME
# =================
def geocode_place(place_text: str) -> Tuple[float, float, str]:
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


# ==========================
# 13-SIGN IAU BAND MAPPING
# ==========================
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
    Build bands by sampling β=0° around the ecliptic at J2000 and detecting constellation changes.
    Output bands cover [0,360). Frame: ecliptic J2000.
    """
    bands: List[Tuple[float, float, str]] = []
    t = Time("J2000")
    current_label: Optional[str] = None
    seg_start = 0.0

    n = int(360.0 / step_deg)
    for i in range(n + 1):
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

    # Merge wrap if needed
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
    _ensure_bands_built()
    x = lon_deg % 360.0
    for s, e, lab in _ECLIPTIC_BANDS:
        if s < 0:
            if (0 <= x < e) or (x + 360.0 >= s and x + 360.0 < e):
                return lab
        else:
            if s <= x < e:
                return lab
    return "Pisces"  # should not happen


# ==========================
# CORE ASTRONOMY (J2000 PIPE)
# ==========================
def ecl_lon_from_app(app) -> float:
    """
    Skyfield apparent RA/Dec → ICRS → ecliptic J2000 → λ (deg 0..360).
    Keep same frame as band map.
    """
    ra, dec, _ = app.radec()  # hours, degrees
    sc_icrs = SkyCoord(ra=float(ra.hours) * u.hourangle, dec=float(dec.degrees) * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


def mean_lunar_node_lon(dt_utc: datetime) -> float:
    t = Time(dt_utc).jd
    T = (t - 2451545.0) / 36525.0
    Omega = 125.04452 - 1934.136261 * T + 0.0020708 * T * T + (T ** 3) / 450000.0
    return Omega % 360.0


def midheaven_lon(dt_utc: datetime, lon_deg: float) -> float:
    obstime = Time(dt_utc)
    lst_deg = obstime.sidereal_time("apparent", longitude=lon_deg * u.deg).deg
    sc_icrs = SkyCoord(ra=lst_deg * u.deg, dec=0.0 * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


def chiron_ecl_lon(dt_utc: datetime, lat_deg: float, lon_deg: float) -> Optional[float]:
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
    ra = float(eph["RA"][0]); dec = float(eph["DEC"][0])
    sc_icrs = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


# ==================
# ASCENDANT SOLVER
# ==================
def compute_ascendant(ts, observer) -> float:
    """
    Solve the ecliptic J2000 longitude of the rising point (β=0°) on the EAST horizon.
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

    # coarse
    best_L, best_abs = 0.0, 1e9
    for L in range(0, 360, 2):
        alt, az = alt_az_of_ecl_lon_j2000(L)
        if 0.0 < az < 180.0:
            a = abs(alt)
            if a < best_abs:
                best_abs, best_L = a, L

    # refine
    step = 1.0
    L0 = best_L
    for _ in range(12):
        candidates = [L0 - step, L0, L0 + step]
        best = (1e9, L0)
        for L in candidates:
            alt, az = alt_az_of_ecl_lon_j2000(L % 360.0)
            if 0.0 < az < 180.0:
                a = abs(alt)
                if a < best[0]:
                    best = (a, L % 360.0)
        L0 = best[1]; step *= 0.5
    return L0 % 360.0


# ===============
# MAIN CHART
# ===============
def compute_chart(name: str, date: str, time_str: Optional[str], place: str) -> Tuple[Dict[str, Any], bool]:
    # 1) Geo + time
    lat, lon, tz_str = geocode_place(place)
    dt_utc, time_estimated = to_utc(date, time_str, tz_str)

    # 2) Ephemeris + observer
    eph, load, eph_name = get_ephemeris()
    ts = load.timescale().from_datetime(dt_utc)
    earth = eph["earth"]
    observer = earth + Topos(latitude_degrees=lat, longitude_degrees=lon)

    # 3) Resolve targets with fallbacks for this ephemeris
    targets = {
        "Sun": pick_key(eph, "sun", "10"),  # 10 is SUN id in many kernels
        "Moon": pick_key(eph, "moon", "301"),
        "Mercury": pick_key(eph, "mercury", "199"),
        "Venus": pick_key(eph, "venus", "299"),
        "Mars": pick_key(eph, "mars", "mars barycenter"),
        "Jupiter": pick_key(eph, "jupiter", "jupiter barycenter"),
        "Saturn": pick_key(eph, "saturn", "saturn barycenter"),
        "Uranus": pick_key(eph, "uranus", "uranus barycenter"),
        "Neptune": pick_key(eph, "neptune", "neptune barycenter"),
        # Prefer Pluto center if present; else barycenter
        "Pluto": pick_key(eph, "pluto", "pluto barycenter"),
    }

    # 4) Ascendant
    asc_lon = compute_ascendant(ts, observer)
    asc_sign = sign_from_ecliptic_lon(asc_lon)
    ascendant = {"lon": asc_lon, "sign": asc_sign, "constellation": asc_sign}

    # 5) Planets → ecliptic J2000
    planets: List[Dict[str, Any]] = []
    for label, key in targets.items():
        app = observer.at(ts).observe(eph[key]).apparent()
        L = ecl_lon_from_app(app)
        sign13 = sign_from_ecliptic_lon(L)
        planets.append({"body": label, "lon": L, "sign": sign13, "constellation": sign13})

    # 6) Nodes (mean)
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = (nn_lon + 180.0) % 360.0
    nodes = {
        "north_node": {"lon": nn_lon, "sign": sign_from_ecliptic_lon(nn_lon)},
        "south_node": {"lon": sn_lon, "sign": sign_from_ecliptic_lon(sn_lon)},
    }

    # 7) Midheaven
    mc_lon = midheaven_lon(dt_utc, lon)
    mc = {"lon": mc_lon, "sign": sign_from_ecliptic_lon(mc_lon)}

    # 8) Chiron (optional)
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

    # 9) Houses: whole-sign from ASC, 12 houses on 13-sign belt (explicit by design)
    SIGNS_13 = [
        "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
        "Libra", "Scorpius", "Ophiuchus", "Sagittarius",
        "Capricornus", "Aquarius", "Pisces",
    ]
    start_idx = SIGNS_13.index(asc_sign)
    houses = [{"house": i + 1, "sign": SIGNS_13[(start_idx + i) % 13]} for i in range(12)]

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
        "ephemeris": eph_name,
        "house_model": "Whole-sign start from ASC, 12 houses on 13-sign belt",
        "time_estimated": time_estimated,
    }
    return payload, time_estimated
