# astro/wheel.py

import math
import numpy as np
import svgwrite
from .signs import SIGNS_13

# Approximate ecliptic spans for 13-sign sectors (used only for drawing labels)
SIGN_SPANS = {
    'Aries': 33.0, 'Taurus': 27.0, 'Gemini': 30.0, 'Cancer': 28.0, 'Leo': 32.0,
    'Virgo': 30.0, 'Libra': 30.0, 'Scorpius': 19.0, 'Ophiuchus': 18.0,
    'Sagittarius': 29.0, 'Capricornus': 30.0, 'Aquarius': 30.0, 'Pisces': 24.0
}

SIGN_COLORS = [
    "#e06666","#ffe599","#93c47d","#6fa8dc","#c27ba0","#ffd966","#f6b26b",
    "#76a5af","#783f04","#674ea7","#e69138","#3d85c6","#b4a7d6"
]

def describe_arc(x, y, radius, start_angle, end_angle):
    # Draw one arc segment of the outer ring
    start_angle, end_angle = math.radians(start_angle-90), math.radians(end_angle-90)
    x1 = x + radius * math.cos(start_angle)
    y1 = y + radius * math.sin(start_angle)
    x2 = x + radius * math.cos(end_angle)
    y2 = y + radius * math.sin(end_angle)
    large_arc = "1" if end_angle - start_angle > math.pi else "0"
    d = f"M{x1},{y1} A{radius},{radius} 0 {large_arc},1 {x2},{y2}"
    return d

def make_svg_wheel(chart_data, outfile, radius=180):
    """
    Draw a simple 13-sign wheel, planet dots, and markers for ASC/MC/NN/SN.
    """
    size = (2*radius+20, 2*radius+20)
    dwg = svgwrite.Drawing(outfile, size=size)
    cx = radius+10
    cy = radius+10
    g = dwg.g(transform=f"translate({cx},{cy})")

    # 1) Sign ring
    angle = 0.0
    for i, sign in enumerate(SIGNS_13):
        span = SIGN_SPANS.get(sign, 30.0)
        next_angle = angle + (span/360.0)*360.0
        g.add(dwg.path(
            d=describe_arc(0, 0, radius, angle, next_angle),
            fill="none",
            stroke=SIGN_COLORS[i % len(SIGN_COLORS)],
            stroke_width=28
        ))
        mid = (angle + next_angle) / 2.0
        mx = 0.80 * radius * math.cos(math.radians(mid-90))
        my = 0.80 * radius * math.sin(math.radians(mid-90))
        g.add(dwg.text(
            sign,
            insert=(mx, my+5),
            text_anchor="middle",
            font_size="14"
        ))
        angle = next_angle

    # 2) Planets (dots + 2-letter label)
    for p in chart_data.get("planets", []):
        lon = float(p["lon"])
        x = 0.70*radius * math.cos(math.radians(lon-90))
        y = 0.70*radius * math.sin(math.radians(lon-90))
        g.add(dwg.circle(center=(x, y), r=7, fill="#111"))
        g.add(dwg.text(
            p["body"][:2],
            insert=(x, y+4),
            text_anchor="middle",
            font_size="10",
            fill="#fff"
        ))

    # Utility for markers
    def mark(tag, lon_deg, color="#e33"):
        ax = 0.98*radius * math.cos(math.radians(lon_deg-90))
        ay = 0.98*radius * math.sin(math.radians(lon_deg-90))
        g.add(dwg.text(tag, insert=(ax, ay-10), text_anchor="middle", font_size="11", fill=color))
        g.add(dwg.line(start=(0,0), end=(ax, ay), stroke=color, stroke_width=2))

    # 3) ASC / MC / NN / SN markers (guard if missing)
    asc = chart_data.get("ascendant", {})
    if "lon" in asc:
        mark("ASC", float(asc["lon"]), "#e33")

    mc = chart_data.get("midheaven", {})
    if "lon" in mc:
        mark("MC", float(mc["lon"]), "#2ad")

    nodes = chart_data.get("nodes", {})
    nn = nodes.get("north_node", {})
    sn = nodes.get("south_node", {})
    if "lon" in nn:
        mark("NN", float(nn["lon"]), "#2ad")
    if "lon" in sn:
        mark("SN", float(sn["lon"]), "#2ad")

    dwg.add(g)
    dwg.save()
