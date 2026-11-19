"""Put the tables of analysis/zoom_report.py between the markers in README.md."""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zoom_report import zoom_md, load, ROOT
p = os.path.join(ROOT, "README.md"); s = open(p).read(); md = zoom_md()
open(os.path.join(ROOT, "analysis", "zoom.md"), "w").write(md)
s2 = re.sub(r"<!-- zoom tables -->.*?<!-- /zoom tables -->", lambda m: "<!-- zoom tables -->\n" + md + "<!-- /zoom tables -->", s, flags=re.S)
Zo = load("zoom_out.json")
if Zo:
    s1, s2_ = Zo["summary"]["step1"], Zo["summary"]["step2"]; f = lambda v: f"{v[0]:.3f} ± {v[2]:.3f}"
    t = ["| window | real: log2 of the growth in cells | its spread | generated | its spread | correlation with the count inside: real | generated |", "|---|---:|---:|---:|---:|---:|---:|",
         f"| twice as wide | {f(s1['real'])} | {s1['real'][1]:.2f} | {f(s1['gen'])} | {s1['gen'][1]:.2f} | {s1['corr_with_inside_real']:.2f} | {s1['corr_with_inside_gen']:.2f} |",
         f"| four times, over twice | {f(s2_['real'])} | {s2_['real'][1]:.2f} | {f(s2_['gen'])} | {s2_['gen'][1]:.2f} | {s2_['corr_with_inside_real']:.2f} | {s2_['corr_with_inside_gen']:.2f} |", "",
         "Generated over real, cell for cell, by ring of the first new level (as shares of the way from its middle to its edge): " + ", ".join(f"{r[0]}–{r[1]}: {r[2]:.2f}" for r in Zo["rings"]) + ". "
         f"{Zo['roots']} cells in {len({r['file'] for r in Zo['rows']})} samples. Run {Zo['levels']} levels outward from a {Zo['G_far']:,} window, the fully relaxed ring of each new level holds "
         f"{Zo['first_half'][0]:,.0f} cells on average over the first half of the levels and {Zo['second_half'][0]:,.0f} over the second, with perimeter per cell {Zo['first_half'][1]:.4f} and "
         f"{Zo['second_half'][1]:.4f}; a level takes {Zo['seconds_per_level']:.2f} s. The generator puts a few per cent too much shape around a window, most of it in the innermost tenth of the new ring, where the links are made. The test is of the amount of shape at each distance, not of whether its arrangement is the right one for the given middle."]
    s2 = re.sub(r"<!-- zoom-out table -->.*?<!-- /zoom-out table -->", lambda m: "<!-- zoom-out table -->\n" + "\n".join(t) + "\n<!-- /zoom-out table -->", s2, flags=re.S)
open(p, "w").write(s2)
