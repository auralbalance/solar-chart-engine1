SIGNS_13 = [
    "Aries","Taurus","Gemini","Cancer","Leo","Virgo",
    "Libra","Scorpius","Ophiuchus","Sagittarius",
    "Capricornus","Aquarius","Pisces"
]

IAU_TO_13 = {
    "Aries":"Aries","Taurus":"Taurus","Gemini":"Gemini","Cancer":"Cancer","Leo":"Leo","Virgo":"Virgo",
    "Libra":"Libra","Scorpius":"Scorpius","Ophiuchus":"Ophiuchus","Sagittarius":"Sagittarius",
    "Capricornus":"Capricornus","Aquarius":"Aquarius","Pisces":"Pisces"
}

def constellation_to_13sign(name:str)->str:
    return IAU_TO_13.get(name, name)

def cycle_from(start:str, n:int=12)->list[str]:
    i = SIGNS_13.index(start)
    return [SIGNS_13[(i+k)%13] for k in range(n)]
