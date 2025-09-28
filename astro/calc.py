from skyfield.api import Loader, Topos
from timezonefinder import TimezoneFinder
from geopy.geocoders import Nominatim
from astroquery.jplhorizons import Horizons
from astropy.time import Time
from datetime import datetime, time as dtime
from pytz import timezone, utc
import numpy as np
from astropy.coordinates import SkyCoord, get_constellation, AltAz, EarthLocation, GeocentricTrueEcliptic
from astropy import units as u
from astropy.time import Time
from .signs import constellation_to_13sign, cycle_from, sign_from_ecliptic_lon

LOADER = Loader("./skyfield_data")  # caches ephemeris here

def _eph():
    return LOADER("de421.bsp")

def geocode(place:str):
    g=Nominatim(user_agent="solar-chart-engine")
    loc=g.geocode(place, timeout=15)
    if not loc: raise ValueError(f"Place '{place}' not found")
    tf=TimezoneFinder(); tz=tf.timezone_at(lng=loc.longitude, lat=loc.latitude)
    if not tz: raise ValueError("Timezone not found")
    return float(loc.latitude), float(loc.longitude), tz

def ecl_lon_from_skyfield(app):
    ecl = app.ecliptic_latlon()
    return float(ecl[1].degrees % 360.0)

def to_utc(date:str, time_str:str|None, tz_str:str):
    if time_str is None: t=dtime(12,0); estimated=True
    else: hh,mm=map(int,time_str.split(':')); t=dtime(hh,mm); estimated=False
    dt=datetime.strptime(date,'%Y-%m-%d').replace(hour=t.hour, minute=t.minute)
    tz=timezone(tz_str); return tz.localize(dt).astimezone(utc), estimated

def ecl_lon_from_skyfield(apparent):
    lat, lon, dist = apparent.ecliptic_latlon()
    return float(lon.degrees % 360.0)

def mean_lunar_node_lon(dt_utc):
    # Meeus approximation, sufficient for chart labeling (mean node)
    # T in Julian centuries since J2000.0
    t = Time(dt_utc).jd
    T = (t - 2451545.0) / 36525.0
    # Longitude of ascending node of the Moon's mean orbit (Ω), degrees
    # Meeus (1998) ch. 22 truncated:
    Omega = 125.04452 - 1934.136261*T + 0.0020708*T*T + (T*T*T)/450000.0
    return Omega % 360.0

def south_node_lon(north_node_lon):
    return (north_node_lon + 180.0) % 360.0

def midheaven_lon(dt_utc, lon_deg):
    # MC longitude: ecliptic longitude where ecliptic intersects local MERIDIAN (upper culmination).
    # Use local apparent sidereal time and obliquity.
    obstime = Time(dt_utc)
    eps = 23.43929111  # mean obliquity deg (OK for labeling)
    # Local sidereal angle in degrees:
    lst = obstime.sidereal_time('apparent', longitude=lon_deg*u.deg).deg  # 0..360
    # MC ecliptic longitude formula:
    # tan(Lmc) = tan(LST) / cos(eps)
    import math
    LST = math.radians(lst)
    ce = math.cos(math.radians(eps))
    Lmc = math.degrees(math.atan2(math.tan(LST), ce)) % 360.0
    return Lmc

def ascendant_lon_accurate(dt_utc, lat, lon):
    # Alt=0°, Az=90° (due East) transformed to ecliptic longitude
    obstime = Time(dt_utc)
    loc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg)
    east_horizon = SkyCoord(alt=0*u.deg, az=90*u.deg, frame=AltAz(obstime=obstime, location=loc))
    ecl = east_horizon.transform_to(GeocentricTrueEcliptic(equinox=obstime))
    return float(ecl.lon.to(u.deg).value % 360.0)

def chiron_ecl_lon(dt_utc, lat_deg, lon_deg):
    """
    Fetch Chiron topocentric apparent RA/Dec from JPL Horizons and return ecliptic longitude (deg).
    Uses observer location (lon,lat, elev=0).
    """
    # Horizons expects "lon,lat,elv" in *degrees* and meters
    loc_str = f"{lon_deg},{lat_deg},0"
    t = Time(dt_utc)
    try:
        obj = Horizons(id='2060', id_type='smallbody', location=loc_str, epochs=t.jd)
        eph = obj.ephemerides()
        ra = float(eph['RA'][0])   # degrees
        dec = float(eph['DEC'][0]) # degrees
        sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        ecl = sc.transform_to(GeocentricTrueEcliptic(equinox=t))
        return float(ecl.lon.to(u.deg).value % 360.0)
    except Exception as e:
        # Fallback to geocentric if topocentric query fails
        obj = Horizons(id='2060', id_type='smallbody', location='geo', epochs=t.jd)
        eph = obj.ephemerides()
        ra = float(eph['RA'][0])
        dec = float(eph['DEC'][0])
        sc = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
        ecl = sc.transform_to(GeocentricTrueEcliptic(equinox=t))
        return float(ecl.lon.to(u.deg).value % 360.0)

def compute_chart(name:str, date:str, time_str:str|None, place:str):
    lat, lon, tz = geocode(place)
    dt_utc, est = to_utc(date, time_str, tz)
    eph=_eph(); ts=LOADER.timescale().from_datetime(dt_utc)
    obs = eph['earth'] + Topos(latitude_degrees=lat, longitude_degrees=lon)
        # ... planets loop above ...

    # --- Chiron ---
    chi_lon = chiron_ecl_lon(dt_utc, lat, lon)
    planets.append({
        'body': 'Chiron',
        'lon': chi_lon,
        'sign': sign_from_ecliptic_lon(chi_lon),
        'constellation': constellation_to_13sign(
            get_constellation(
                SkyCoord(
                    lon=chi_lon*u.deg, lat=0*u.deg,
                    frame=GeocentricTrueEcliptic(equinox=Time(dt_utc))
                ).transform_to('icrs')
            )
        )
    })

    # Nodes
    nn_lon = mean_lunar_node_lon(dt_utc)
    sn_lon = south_node_lon(nn_lon)
    nodes = {
        'north_node': {'lon': nn_lon, 'sign': sign_from_ecliptic_lon(nn_lon)},
        'south_node': {'lon': sn_lon, 'sign': sign_from_ecliptic_lon(sn_lon)}
    }

    # Midheaven (MC)
    mc_lon = midheaven_lon(dt_utc, lon)
    mc = {'lon': mc_lon, 'sign': sign_from_ecliptic_lon(mc_lon)}

    bodies = {
        'Sun': eph['sun'],
        'Moon': eph['moon'],
        'Mercury': eph['mercury'],
        'Venus': eph['venus'],
        'Mars': eph['mars'],
        'Jupiter': eph['jupiter barycenter'],
        'Saturn': eph['saturn barycenter'],
        'Uranus': eph['uranus barycenter'],
        'Neptune': eph['neptune barycenter'],
        'Pluto': eph['pluto barycenter']
    }

    planets=[]
    for name_b, target in bodies.items():
        app = obs.at(ts).observe(target).apparent()
        lon_ecl = ecl_lon_from_skyfield(app)
        ra, dec, _ = app.radec()
        sc = SkyCoord(ra=ra.hours*u.hourangle, dec=dec.degrees*u.deg, frame='icrs')
        const = get_constellation(sc)
        # Map Serpens -> Ophiuchus then ALSO compute sign from ecliptic longitude
        const_mapped = constellation_to_13sign(const)
        sign = sign_from_ecliptic_lon(lon_ecl)
        planets.append({'body': name_b, 'lon': lon_ecl, 'sign': sign, 'constellation': const_mapped})

    asc_lon = ascendant_lon_accurate(dt_utc, lat, lon)
    ecl_sc = SkyCoord(lon=asc_lon*u.deg, lat=0*u.deg, frame=GeocentricTrueEcliptic(equinox=Time(dt_utc)))
    icrs_sc = ecl_sc.transform_to('icrs')
    asc_const = get_constellation(icrs_sc)
    asc_sign = sign_from_ecliptic_lon(asc_lon)
    ascendant = {'lon': asc_lon, 'sign': asc_sign, 'constellation': constellation_to_13sign(asc_const)}


    houses = [{'house': i+1, 'sign': cycle_from(asc_sign,12)[i]} for i in range(12)]

    return {
        'name': name,
        'datetime_utc': dt_utc.strftime('%Y-%m-%dT%H:%M:%SZ'),
        'location': {'lat': lat, 'lon': lon, 'tz': tz},
        'ascendant': ascendant,
        'midheaven': mc,
        'nodes': nodes,
        'planets': planets,
        'houses': houses
    }, est
