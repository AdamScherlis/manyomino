"""Put the tables of analysis/zoom_report.py between the markers in README.md."""
import os, re, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zoom_report import zoom_md, ROOT
p = os.path.join(ROOT, "README.md"); s = open(p).read(); md = zoom_md()
open(os.path.join(ROOT, "analysis", "zoom.md"), "w").write(md)
s2 = re.sub(r"<!-- zoom tables -->.*?<!-- /zoom tables -->", lambda m: "<!-- zoom tables -->\n" + md + "<!-- /zoom tables -->", s, flags=re.S)
assert s2 != s or md in s; open(p, "w").write(s2)
