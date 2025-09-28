import math, svgwrite
from .signs import SIGNS_13

# ecliptic spans approx; for labels only (not classification)
SPANS = {'Aries':25,'Taurus':37,'Gemini':29,'Cancer':21,'Leo':36,'Virgo':44,
         'Libra':23,'Scorpius':22,'Ophiuchus':18,'Sagittarius':31,'Capricornus':27,'Aquarius':24,'Pisces':38}
COLORS = ["#ffb3ba","#ffdfba","#ffffba","#baffc9","#bae1ff","#d7bde2","#fadbd8",
          "#f5b7b1","#a9dfbf","#aed6f1","#f9e79f","#a3e4d7","#d2b4de"]

def _arc(cx,cy,r,a1,a2):
    from math import radians,cos,sin,pi
    s,e = radians(a1-90), radians(a2-90)
    x1,y1 = cx+r*cos(s), cy+r*sin(s); x2,y2 = cx+r*cos(e), cy+r*sin(e)
    large = 1 if (e-s)% (2*pi) > pi else 0
    return f"M{x1},{y1} A{r},{r} 0 {large},1 {x2},{y2}"

def make_svg_wheel(data:dict, outfile:str, R:int=200):
    size=(2*R+40,2*R+40); cx=cy=R+20
    dwg=svgwrite.Drawing(outfile,size=size); g=dwg.g()
    g.add(dwg.circle(center=(cx,cy),r=R,fill="white",stroke="#222",stroke_width=2))
    total=sum(SPANS[s] for s in SIGNS_13); scale=360/total; ang=0.0
    for i,s in enumerate(SIGNS_13):
        span=SPANS[s]*scale; path=_arc(cx,cy,R-10,ang,ang+span)
        g.add(dwg.path(d=path, fill="none", stroke=COLORS[i%len(COLORS)], stroke_width=18))
        mid=math.radians((ang+ang+span)/2 - 90); lx=cx+(R-36)*math.cos(mid); ly=cy+(R-36)*math.sin(mid)
        g.add(dwg.text(s, insert=(lx,ly+4), text_anchor="middle", font_size="12px", fill="#111")); ang+=span
    for p in data.get("planets",[]):
        lon=p["lon"]%360; th=math.radians(lon-90); x=cx+(R-55)*math.cos(th); y=cy+(R-55)*math.sin(th)
        g.add(dwg.circle(center=(x,y), r=6, fill="#222")); g.add(dwg.text(p["body"][:2], insert=(x,y+4), fill="#fff", font_size="9px", text_anchor="middle"))
    asc=data.get("ascendant",{}).get("lon",0)%360; ax=cx+(R-6)*math.cos(math.radians(asc-90)); ay=cy+(R-6)*math.sin(math.radians(asc-90))
    g.add(dwg.line(start=(cx,cy), end=(ax,ay), stroke="#d33", stroke_width=2)); g.add(dwg.text("ASC", insert=(ax,ay-10), text_anchor="middle", font_size="10px", fill="#d33"))
    dwg.add(g); dwg.save(); return outfile
