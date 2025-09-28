# astro/calc.py
from __future__ import annotations
import os, json
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
    raise RuntimeError("No ephemeris available (tried de441.bsp, de440.bsp, de421.bsp).")


def pick_key(eph, preferred: str, fallback: str) -> str:
    """Return preferred target key if present; else fallback; else raise."""
    try:
        _ = eph[preferred]
        return preferred
    except Exception:
        try:
            _ = eph[fallback]
            return fallback
        except Exception:
            raise KeyError(f"Neither '{preferred}' nor '{fallback}' available in ephemeris.")


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
# ZODIAC PROFILES (J2000)
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
SIGNS_13 = [
    "Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
    "Libra", "Scorpius", "Ophiuchus", "Sagittarius",
    "Capricornus", "Aquarius", "Pisces",
]

def constellation_to_13sign(iau_name: str) -> str:
    return CONSTELLATION_TO_13.get(iau_name, iau_name)

# Built-in CUSTOM_13 tuned to your MTZ screenshot
# Key differences vs IAU:
# - Scorpius ends at 252° so Mercury 250.27° is Scorpius
# - Sagittarius ends at 306° so North Node 301.47° is Sagittarius
DEFAULT_CUSTOM_BANDS: List[Tuple[float, float, str]] = [
    (0.0,   30.0,  "Aries"),
    (30.0,  72.0,  "Taurus"),
    (72.0,  90.0,  "Gemini"),
    (90.0,  138.0, "Cancer"),
    (138.0, 168.0, "Leo"),
    (168.0, 204.0, "Virgo"),
    (204.0, 241.0, "Libra"),
    (241.0, 252.0, "Scorpius"),
    (252.0, 266.0, "Ophiuchus"),
    (266.0, 306.0, "Sagittarius"),
    (306.0, 326.0, "Capricornus"),
    (326.0, 353.0, "Aquarius"),
    (353.0, 360.0, "Pisces"),
]

# In-memory band store
_ECLIPTIC_BANDS: List[Tuple[float, float, str]] = []  # (start, end, name)
_BUILT_FOR_PROFILE: Optional[str] = None


def _build_iau_bands(step_deg: float = 0.2) -> List[Tuple[float, float, str]]:
    """
    Build bands by sampling β=0° around the ecliptic at J2000 and detecting constellation changes.
    """
    bands: List[Tuple[float, float, str]] = []
    t = Time("J2000")
    current_label: Optional[str] = None
    seg_start = 0.0

    n = int(360.0 / step_deg)
    for i in range(n + 1):
        L = (i * step_deg) % 360.0
        sc = SkyCoord(lon=L * u.deg, lat=0.0 * u.deg,
                      frame=GeocentricTrueEcliptic(equinox=t)).icrs
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

    if bands and bands[0][2] == bands[-1][2]:
        first_s, first_e, lab = bands[0]
        last_s, last_e, _ = bands[-1]
        bands = [(last_s - 360.0, first_e, lab)] + bands[1:-1]

    bands.sort(key=lambda x: x[0])
    return [(s, e, name) for (s, e, name) in bands]


def _load_custom_bands() -> List[Tuple[float, float, str]]:
    """
    Load custom bands from env CUSTOM_BANDS (JSON string).
    We intentionally ignore any file to avoid stale overrides.
    Format:
      [{"name":"Aries","start":0.0,"end":30.0}, ...]
    """
    raw = os.getenv("CUSTOM_BANDS")
    if not raw:
        return DEFAULT_CUSTOM_BANDS
    data = json.loads(raw)
    bands: List[Tuple[float, float, str]] = []
    for item in data:
        bands.append((float(item["start"]), float(item["end"]), str(item["name"])))
    bands.sort(key=lambda x: x[0])
    return bands


def _ensure_bands(profile: str):
    global _ECLIPTIC_BANDS, _BUILT_FOR_PROFILE
    if (_BUILT_FOR_PROFILE != profile) or (not _ECLIPTIC_BANDS):
        if profile == "IAU_13":
            _ECLIPTIC_BANDS = _build_iau_bands(step_deg=0.2)
        elif profile == "CUSTOM_13":
            _ECLIPTIC_BANDS = _load_custom_bands()
        else:
            raise ValueError(f"Unknown ZODIAC_PROFILE '{profile}'")
        _BUILT_FOR_PROFILE = profile


def sign_from_ecliptic_lon(lon_deg: float, profile: str) -> str:
    _ensure_bands(profile)
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
# CORE ASTRONOMY (J2000 PIPE)
# ==========================
def ecl_lon_from_app(app) -> float:
    """
    Skyfield apparent RA/Dec → ICRS → ecliptic J2000 → λ (deg 0..360).
    """
    ra, dec, _ = app.radec()  # hours, degrees
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
# HOUSES (13 default)
# ===============
def build_houses(asc_sign: str, mode: str) -> List[Dict[str, Any]]:
    """
    Whole-sign logic, aligned to 13-sign belt.
    mode="12": 12 houses on 13-sign belt (ASC sign = House 1, then next 11 signs)
    mode="13": 13 houses on 13-sign belt (ASC sign = House 1, full cycle)
    """
    if mode not in ("12", "13"):
        mode = "13"
    start_idx = SIGNS_13.index(asc_sign)
    count = 13 if mode == "13" else 12
    return [{"house": i + 1, "sign": SIGNS_13[(start_idx + i) % 13]} for i in range(count)]


def house_number_for_sign(sign: str, asc_sign: str, mode: str) -> int:
    asc_i = SIGNS_13.index(asc_sign)
    sign_i = SIGNS_13.index(sign)
    delta = (sign_i - asc_i) % 13
    count = 13 if mode == "13" else 12
    return (delta % count) + 1


# ===============
# MAIN CHART
# ===============
def compute_chart(
    name: str,
    date: str,
    time_str: Optional[str],
    place: str,
    house_mode: Optional[str] = None,
    zodiac_profile: Optional[str] = None,
) -> Tuple[Dict[str, Any], bool]:
    """
    house_mode: "12" or "13". If None, env HOUSE_MODE; default "13".
    zodiac_profile: "IAU_13" or "CUSTOM_13". If None, env ZODIAC_PROFILE; default "CUSTOM_13".
    """
    mode = (house_mode or os.getenv("HOUSE_MODE") or "13").strip().lower()
    mode = "13" if mode in ("13", "h13", "thirteen") else "12"

    profile = (zodiac_profile or os.getenv("ZODIAC_PROFILE") or "CUSTOM_13").strip().upper()
    if profile not in ("IAU_13", "CUSTOM_13"):
        profile = "CUSTOM_13"

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

    # 4) Ascendant
    asc_lon = compute_ascendant(ts, observer)
    asc_sign = sign_from_ecliptic_lon(asc_lon, profile)
    ascendant = {"lon": asc_lon, "sign": asc_sign, "constellation": asc_sign, "house": 1}

    # 5) Planets → ecliptic J2000 + houses
    planets: List[Dict[str, Any]] = []
    for label, key in targets.items():
        app = observer.at(ts).observe(eph[key]).apparent()
        L = ecl_lon_from_app(app)
        sign13 = sign_from_ecliptic_lon(L, profile)
        planets.append({
            "body": label,
            "lon": L,
            "sign": sign13,
            "constellation": sign13,
            "house": house_number_for_sign(sign13, asc_sign, mode),
        })

    # 6) Nodes (mean)
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = (nn_lon + 180.0) % 360.0
    nn_sign = sign_from_ecliptic_lon(nn_lon, profile)
    sn_sign = sign_from_ecliptic_lon(sn_lon, profile)
    nodes = {
        "north_node": {"lon": nn_lon, "sign": nn_sign, "house": house_number_for_sign(nn_sign, asc_sign, mode)},
        "south_node": {"lon": sn_lon, "sign": sn_sign, "house": house_number_for_sign(sn_sign, asc_sign, mode)},
    }

    # 7) Midheaven
    mc_lon = midheaven_lon(dt_utc, lon)
    mc_sign = sign_from_ecliptic_lon(mc_lon, profile)
    mc = {"lon": mc_lon, "sign": mc_sign, "house": house_number_for_sign(mc_sign, asc_sign, mode)}

    # 8) Chiron (optional)
    chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
    if chi_lon is not None:
        chi_sign = sign_from_ecliptic_lon(chi_lon, profile)
        planets.append(
            {
                "body": "Chiron",
                "lon": chi_lon,
                "sign": chi_sign,
                "constellation": chi_sign,
                "house": house_number_for_sign(chi_sign, asc_sign, mode),
            }
        )

    # 9) Houses
    houses = build_houses(asc_sign, mode)

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
        "house_model": f"Whole-sign; {len(houses)} houses on 13-sign belt",
        "house_mode": mode,
        "zodiac_profile": profile,
        "time_estimated": time_estimated,
    }
    return payload, time_estimated
