"""Draw a Ventsim network to SVG (no plotting library needed)."""
from __future__ import annotations

import os
from html import escape

from . import vsm

VIEWS = {"plan": ("X", "Y"), "section_x": ("X", "Z"), "section_y": ("Y", "Z")}
_STOPS = [(0.0, (49, 99, 206)), (0.33, (40, 170, 110)), (0.66, (240, 190, 40)), (1.0, (214, 48, 39))]


def _ramp(t: float) -> str:
    t = min(max(t, 0.0), 1.0)
    for (t0, c0), (t1, c1) in zip(_STOPS, _STOPS[1:]):
        if t <= t1:
            f = (t - t0) / (t1 - t0)
            return "#%02x%02x%02x" % tuple(int(a + (b - a) * f) for a, b in zip(c0, c1))
    return "#d63027"


def svg(model: vsm.Model, view: str = "plan", color_by: str = "Velocity m/s",
        highlight: set[str] | None = None, width_px: int = 1600) -> tuple[str, int]:
    """SVG text and number of airways drawn (airways of the simulated stage)."""
    if view not in VIEWS:
        raise ValueError(f"view must be one of {list(VIEWS)}")
    h, v = VIEWS[view]
    ws = [w for w in vsm.airways(model) if w["In Simulated Stage"]]
    if not ws:
        raise ValueError("no airways to draw")
    if color_by not in ws[0]:
        raise ValueError(f"color_by must be one of {list(ws[0])}")
    hl = highlight or set()

    def pt(w, end):
        return (w[h + end] or 0.0, w[v + end] or 0.0)

    xs = [pt(w, e)[0] for w in ws for e in "12"]
    ys = [pt(w, e)[1] for w in ws for e in "12"]
    x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
    pad, W = 40, max(int(width_px), 300)
    sx = (W - 2 * pad) / (x1 - x0) if x1 > x0 else 1.0
    H = min(max(int((y1 - y0) * sx) + 2 * pad, 300), 4 * W)
    sc = min(sx, (H - 2 * pad) / (y1 - y0)) if y1 > y0 else sx

    def X(x):
        return pad + (x - x0) * sc

    def Y(y):
        return H - pad - (y - y0) * sc

    vals = sorted(abs(w[color_by]) for w in ws if isinstance(w[color_by], (int, float)))
    lo = vals[int(0.02 * (len(vals) - 1))] if vals else 0.0
    hi = vals[int(0.98 * (len(vals) - 1))] if vals else 1.0
    qmax = max((abs(w["Quantity"] or 0) for w in ws), default=0) or 1.0

    def colour(val):
        if not isinstance(val, (int, float)):
            return "#bbbbbb"
        return _ramp(0.0 if hi <= lo else (abs(val) - lo) / (hi - lo))

    out = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
           f'<rect width="{W}" height="{H}" fill="white"/>', '<g stroke-linecap="round" fill="none">']
    for w in sorted(ws, key=lambda w: abs(w["Quantity"] or 0)):
        (ax, ay), (bx, by) = pt(w, "1"), pt(w, "2")
        sw = 0.6 + 4.0 * (abs(w["Quantity"] or 0) / qmax) ** 0.5
        tip = escape(f'{w["Unique No"]} {w["Branch Name"]} {color_by}={w[color_by]} Q={w["Quantity"]}')
        out.append(f'<line x1="{X(ax):.1f}" y1="{Y(ay):.1f}" x2="{X(bx):.1f}" y2="{Y(by):.1f}" '
                   f'stroke="{colour(w[color_by])}" stroke-width="{sw:.2f}"><title>{tip}</title></line>')
    for w in ws:
        if w["Unique No"] in hl:
            (ax, ay), (bx, by) = pt(w, "1"), pt(w, "2")
            out.append(f'<line x1="{X(ax):.1f}" y1="{Y(ay):.1f}" x2="{X(bx):.1f}" y2="{Y(by):.1f}" '
                       f'stroke="#d0189a" stroke-width="6" opacity="0.8"/>')
    out.append("</g>")
    for w in ws:
        if w["Fan Type"] not in ("", "0"):
            (ax, ay), (bx, by) = pt(w, "1"), pt(w, "2")
            tip = escape(f'fan {w["Unique No"]} {w["Fan Name"]}')
            out.append(f'<circle cx="{X((ax + bx) / 2):.1f}" cy="{Y((ay + by) / 2):.1f}" r="6" '
                       f'fill="white" stroke="black" stroke-width="2"><title>{tip}</title></circle>')
    label = escape(f"{os.path.basename(model.path)}: {view}, colour = {color_by} "
                   f"(blue {lo:.3g} to red {hi:.3g}), line width = airflow, circles = fans")
    out.append(f'<text x="{pad}" y="24" font-family="sans-serif" font-size="16">{label}</text>')
    out.append("</svg>")
    return "\n".join(out), len(ws)
