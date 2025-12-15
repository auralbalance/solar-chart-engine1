# astro/calc.py
from __future__ import annotations
import os, time
from datetime import datetime, time as dtime
from typing import Tuple, Dict, Any, Optional, List
from functools import lru_cache

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
# EPHEMERIS (cached)
# =========================
@lru_cache(maxsize=1)
def get_ephemeris():
    """
    Load robust ephemeris, cached in /tmp/skyfield-data.
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
    raise RuntimeError("No ephemeris available")


def pick_key(eph, preferred: str, fallback: str) -> str:
    try:
        _ = eph[preferred]
        return preferred
    except Exception:
        try:
            _ = eph[fallback]
            return fallback
        except Exception:
            raise KeyError(f"Neither '{preferred}' nor '{fallback}' in ephemeris")


# =================
# GEO / TIME  (TTL cache)
# =================
_GEO_CACHE: Dict[str, Tuple[float, float, str, float]] = {}
_GEO_TTL = float(os.getenv("GEOCODE_TTL_SECONDS", "604800"))  # 7 days

def geocode_place(place_text: str) -> Tuple[float, float, str]:
    """
    Resolve place to (lat, lon, tz). In-memory cache with TTL.
    """
    key = place_text.strip().lower()
    now = time.time()
    hit = _GEO_CACHE.get(key)
    if hit and (now - hit[3] < _GEO_TTL):
        return hit[0], hit[1], hit[2]

    geocoder = Nominatim(user_agent="solar-chart-api", timeout=15)
    loc = geocoder.geocode(place_text, addressdetails=False, language="en")
    if not loc:
        raise ValueError(f"Place '{place_text}' not found")
    lat = float(loc.latitude)
    lon = float(loc.longitude)

    tf = TimezoneFinder()
    tz = tf.timezone_at(lng=lon, lat=lat)
    if not tz:
        raise ValueError("Timezone not found for the given location")

    _GEO_CACHE[key] = (lat, lon, tz, now)
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


# ==========================
# CANONICAL 13-SIGN NAMES
# ==========================
SIGNS_13 = [
    "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
    "Libra","Scorpius","Ophiuchus","Sagittarius",
    "Capricornus","Aquarius","Pisces",
]

# 12-SIGN SIDEREAL (30° equal divisions, no Ophiuchus)
SIGNS_12_SIDEREAL = [
    "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
    "Libra","Scorpius","Sagittarius","Capricornus","Aquarius","Pisces",
]

def sign_12_sidereal(lon_deg: float) -> str:
    """
    12-sign sidereal: equal 30° divisions.
    0-30° = Aries, 30-60° = Taurus, etc.
    Uses Lahiri ayanamsa approximation (~24° for 2000).
    """
    # Lahiri ayanamsa for J2000.0 ≈ 23.85°
    ayanamsa = 23.85
    sidereal_lon = (lon_deg - ayanamsa) % 360.0
    index = int(sidereal_lon / 30.0) % 12
    return SIGNS_12_SIDEREAL[index]

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
def canon_sign(name: str) -> str:
    n = str(name).strip()
    low = n.lower()
    if low.startswith("scorp"):
        return "Scorpius"
    if low.startswith("capricor"):
        return "Capricornus"
    if low.startswith("ophiu"):
        return "Ophiuchus"
    for s in SIGNS_13:
        if low == s.lower():
            return s
    return n[:1].upper() + n[1:].lower()


# ==========================
# IAU ECLIPTIC BANDS (β=0°)
# ==========================
_ECLIPTIC_BANDS: List[Tuple[float, float, str]] = []  # (start, end, sign)

def _build_ecliptic_bands(step_deg: float = 0.05) -> List[Tuple[float, float, str]]:
    bands: List[Tuple[float, float, str]] = []
    t = Time("J2000")
    current_label: Optional[str] = None
    seg_start = 0.0
    n = int(360.0 / step_deg)

    for i in range(n + 1):
        L = (i * step_deg) % 360.0
        sc = SkyCoord(lon=L * u.deg, lat=0.0 * u.deg,
                      frame=GeocentricTrueEcliptic(equinox=t)).icrs
        iau = get_constellation(sc)
        label = CONSTELLATION_TO_13.get(iau)
        if label is None:
            label = current_label if current_label is not None else "Pisces"
        if current_label is None:
            current_label = label
            seg_start = L
        elif label != current_label:
            bands.append((seg_start, L, current_label))
            seg_start = L
            current_label = label

    if current_label is not None:
        bands.append((seg_start, 360.0, current_label))

    if bands and bands[0][2] == bands[-1][2]:
        first_s, first_e, lab = bands[0]
        last_s, last_e, _ = bands[-1]
        bands = [(last_s - 360.0, first_e, lab)] + bands[1:-1]

    bands.sort(key=lambda x: x[0])
    return bands


def _ensure_bands():
    global _ECLIPTIC_BANDS
    if not _ECLIPTIC_BANDS:
        _ECLIPTIC_BANDS = _build_ecliptic_bands(step_deg=0.05)


def get_ecliptic_bands() -> List[Dict[str, Any]]:
    """
    Return the ecliptic band table used for labeling.
    Each item: {"start": deg, "end": deg, "sign": name}
    """
    _ensure_bands()
    return [{"start": s, "end": e, "sign": lab} for (s, e, lab) in _ECLIPTIC_BANDS]


def sign_from_ecliptic_lon(lon_deg: float) -> str:
    _ensure_bands()
    x = lon_deg % 360.0
    for s, e, lab in _ECLIPTIC_BANDS:
        if s <= e:
            if s <= x < e:
                return lab
        else:
            if x >= s or x < e:
                return lab
    return "Pisces"


# ==========================
# CORE ASTRONOMY (λ on J2000 ecliptic)
# ==========================
def ecl_lon_from_app(app) -> float:
    ra, dec, _ = app.radec()
    sc_icrs = SkyCoord(ra=float(ra.hours) * u.hourangle,
                       dec=float(dec.degrees) * u.deg,
                       frame=ICRS())
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
    ra = float(eph["RA"][0])
    dec = float(eph["DEC"][0])
    sc_icrs = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame=ICRS())
    ecl = sc_icrs.transform_to(GeocentricTrueEcliptic(equinox=Time("J2000")))
    return float(ecl.lon.to(u.deg).value % 360.0)


# ==================
# ASCENDANT SOLVER
# ==================
def compute_ascendant(ts, observer) -> float:
    """
    Solve ecliptic J2000 longitude of the rising point on the EAST horizon.
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

    best_L, best_abs = 0.0, 1e9
    for L in range(0, 360, 2):
        alt, az = alt_az_of_ecl_lon_j2000(L)
        if 0.0 < az < 180.0:
            a = abs(alt)
            if a < best_abs:
                best_abs, best_L = a, L

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
        L0 = best[1]
        step *= 0.5
    return L0 % 360.0


# ===============
# 13-HOUSE MODEL
# ===============
def build_houses_13(asc_sign: str) -> List[Dict[str, Any]]:
    asc_sign = canon_sign(asc_sign)
    start_idx = SIGNS_13.index(asc_sign)
    return [{"house": i + 1, "sign": SIGNS_13[(start_idx + i) % 13]} for i in range(13)]

def house_number_for_sign_13(sign: str, asc_sign: str) -> int:
    sign = canon_sign(sign)
    asc_sign = canon_sign(asc_sign)
    asc_i = SIGNS_13.index(asc_sign)
    sign_i = SIGNS_13.index(sign)
    return ((sign_i - asc_i) % 13) + 1


# ===============
# MAIN CHART  (true sky, 13 houses)
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

    # 3) Ascendant
    asc_lon = compute_ascendant(ts, observer)
    asc_sign_12 = sign_12_sidereal(asc_lon)  # 12-sign sidereal
    asc_sign_13 = sign_from_ecliptic_lon(asc_lon)  # 13-sign constellation
    ascendant = {"lon": asc_lon, "sign": asc_sign_12, "constellation": asc_sign_13, "house": 1}

    # 4) Planets — classify by ecliptic longitude bands
    targets = {
        "Sun": pick_key(eph, "sun", "10"),
        "Moon": pick_key(eph, "moon", "301"),
        "Mercury": pick_key(eph, "mercury", "199"),
        "Venus": pick_key(eph, "venus", "299"),
        "Mars": pick_key(eph, "mars", "mars barycenter"),
        "Jupiter": pick_key(eph, "jupiter", "jupiter barycenter"),
        "Saturn": pick_key(eph, "saturn", "saturn barycenter"),
        "Uranus": pick_key(eph, "uranus", "uranus barycenter"),
        "Neptune": pick_key(eph, "neptune", "neptune barycenter"),
        "Pluto": pick_key(eph, "pluto", "pluto barycenter"),
    }

    planets: List[Dict[str, Any]] = []
    for label, key in targets.items():
        app = observer.at(ts).observe(eph[key]).apparent()
        L = ecl_lon_from_app(app)
        sign_12 = sign_12_sidereal(L)  # 12-sign sidereal (30° divisions)
        sign_13 = sign_from_ecliptic_lon(L)  # 13-sign constellation (IAU boundaries)
        planets.append({
            "body": label,
            "lon": L,
            "sign": sign_12,  # 12-sign sidereal
            "constellation": sign_13,  # 13-sign constellation
            "house": house_number_for_sign_13(sign_13, asc_sign_13),
        })

    # 5) Nodes — mean; classify by both systems
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = (nn_lon + 180.0) % 360.0
    nn_sign_12 = sign_12_sidereal(nn_lon)
    nn_sign_13 = sign_from_ecliptic_lon(nn_lon)
    sn_sign_12 = sign_12_sidereal(sn_lon)
    sn_sign_13 = sign_from_ecliptic_lon(sn_lon)
    nodes = {
        "north_node": {"lon": nn_lon, "sign": nn_sign_12, "constellation": nn_sign_13, "house": house_number_for_sign_13(nn_sign_13, asc_sign_13)},
        "south_node": {"lon": sn_lon, "sign": sn_sign_12, "constellation": sn_sign_13, "house": house_number_for_sign_13(sn_sign_13, asc_sign_13)},
    }

    # 6) Midheaven — classify by both systems
    mc_lon = midheaven_lon(dt_utc, lon)
    mc_sign_12 = sign_12_sidereal(mc_lon)
    mc_sign_13 = sign_from_ecliptic_lon(mc_lon)
    mc = {"lon": mc_lon, "sign": mc_sign_12, "constellation": mc_sign_13, "house": house_number_for_sign_13(mc_sign_13, asc_sign_13)}

    # 7) Chiron (Horizons) — classify by both systems
    chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
    if chi_lon is not None:
        chi_sign_12 = sign_12_sidereal(chi_lon)
        chi_sign_13 = sign_from_ecliptic_lon(chi_lon)
        planets.append({
            "body": "Chiron",
            "lon": chi_lon,
            "sign": chi_sign_12,  # 12-sign sidereal
            "constellation": chi_sign_13,  # 13-sign constellation
            "house": house_number_for_sign_13(chi_sign_13, asc_sign_13),
        })

    # 8) Houses
    houses = build_houses_13(asc_sign_13)

    payload = {
        "name": name,
        "datetime_utc": dt_utc.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "location": {"lat": lat, "lon": lon, "tz": tz_str},
        "ascendant": ascendant,
        "midheaven": mc,
        "nodes": nodes,
        "planets": planets,
        "houses": houses,
        "frames": {"positions": "ICRS→Ecliptic J2000 for λ", "signs": "IAU ecliptic bands (β=0°)"},
        "ephemeris": eph_name,
        "house_model": "13 houses, whole-sign, one house per constellation sign",
        "time_estimated": time_estimated,
    }
    return payload, time_estimated
