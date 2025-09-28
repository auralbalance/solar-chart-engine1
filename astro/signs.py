# astro/signs.py

SIGNS_13 = [
    "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
    "Libra","Scorpius","Ophiuchus","Sagittarius",
    "Capricornus","Aquarius","Pisces"
]

IAU_TO_13 = {
    "Aries":"Aries","Taurus":"Taurus","Gemini":"Gemini","Cancer":"Cancer","Leo":"Leo","Virgo":"Virgo",
    "Libra":"Libra","Scorpius":"Scorpius","Ophiuchus":"Ophiuchus","Sagittarius":"Sagittarius",
    "Capricornus":"Capricornus","Aquarius":"Aquarius","Pisces":"Pisces",
    # Force merges for zodiac practice:
    "Serpens":"Ophiuchus"  # Serpens shares Ophiuchus region along the ecliptic
}

def constellation_to_13sign(name:str)->str:
    return IAU_TO_13.get(name, name)

def cycle_from(start:str, n:int=12)->list[str]:
    i = SIGNS_13.index(start)
    return [SIGNS_13[(i+k)%13] for k in range(n)]

# Ecliptic longitude spans (degrees, IAU-like along the ecliptic; simple MVP)
# These are approximate sector edges in ecliptic longitude increasing from 0° Aries.
# Source: widely used 13-sign boundary approximations on the ecliptic.
SIGN_EDGES = [
    ("Aries",        0.0),
    ("Taurus",      33.0),
    ("Gemini",      60.0),
    ("Cancer",      90.0),
    ("Leo",        118.0),
    ("Virgo",      150.0),
    ("Libra",      180.0),
    ("Scorpius",   210.0),
    ("Ophiuchus",  229.0),
    ("Sagittarius",247.0),
    ("Capricornus",276.0),
    ("Aquarius",   306.0),
    ("Pisces",     336.0),
    ("Aries",     360.0)  # wrap
]

def sign_from_ecliptic_lon(lon_deg: float) -> str:
    x = lon_deg % 360.0
    for i in range(len(SIGN_EDGES)-1):
        name, a = SIGN_EDGES[i]
        _, b = SIGN_EDGES[i+1]
        if a <= x < b:
            return name
    return "Aries"
