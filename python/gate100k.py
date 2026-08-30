"""Pool ALL n=100k chain segments by lineage (bar* vs rect*) and run the
convergence gate.  Usage: python3 gate100k.py [burn]"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fourseed import run

ROOT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "prod")
allcsv = glob.glob(os.path.join(ROOT, "n100000_*.csv"))
stringy = sorted(f for f in allcsv if "_bar" in f)
compact = sorted(f for f in allcsv if "_rect" in f)
burn = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
print(f"{len(stringy)} stringy + {len(compact)} compact segments")
sys.exit(0 if run(stringy, compact, burn) else 1)
