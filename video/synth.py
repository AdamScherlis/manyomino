"""Soundtrack for the music video, synthesised from the chain's own event log.

  python3 video/synth.py events.json out.wav

138 BPM.  A mellow half-time opening builds into two trance drops.  What the chain plays:
  * every accepted relocation is a tiny high tick;
  * on each transplant burst the largest accepted branch is a mallet note - the bigger the branch,
    the lower and longer the note (A minor pentatonic), panned to where the branch landed;
  * in the build-ups and drops the same branch also fires the lead synth: one stab per sixteenth
    in which a branch moved (no move, no note), brighter and louder the bigger the branch, and
    stepping up the chord for branches of 24+ and 150+ cells;
  * branches of 350+ cells get a low boom in the quiet sections and a zap in the drops.
Pad, bass, drums, risers and the chord sequence are a fixed arrangement around that.
"""
import json, sys, wave
import numpy as np

SR, FPS, BPM = 44100, 30, 138
ev = json.load(open(sys.argv[1]))
NF = len(ev); NS = int(NF / FPS * SR) + 1
rng = np.random.default_rng(7)
SIXT = 15.0 / BPM; BEAT, BAR = 4 * SIXT, 16 * SIXT
SECT = {'intro': (0, 8), 'build': (8, 16), 'drop1': (16, 32), 'break': (32, 40), 'drop2': (40, 56), 'outro': (56, 62), 'end': (62, 64)}
NBARS = 64
DROPS = ('drop1', 'drop2')
mtof = lambda m: 440.0 * 2 ** ((m - 69) / 12)
tt = np.arange(NS) / SR


def sect(bar):
    for k, (a, b) in SECT.items():
        if a <= bar < b: return k
    return 'end'


class Bus:
    def __init__(self): self.a = np.zeros((NS + SR * 6, 2), np.float32)
    def add(self, t, mono, pan=0.0, gain=1.0):
        i = int(round(t * SR))
        if i < 0 or i >= NS: return
        th = (pan + 1) * np.pi / 4
        self.a[i:i + len(mono), 0] += mono * (gain * np.cos(th)); self.a[i:i + len(mono), 1] += mono * (gain * np.sin(th))
    def add2(self, t, st, gain=1.0):
        i = int(round(t * SR))
        if i < 0 or i >= NS: return
        self.a[i:i + len(st)] += st * gain
    def out(self): return self.a[:NS].astype(np.float64)


def filt(x, lo=None, hi=None, order=2):
    """Zero-phase Butterworth-shaped filter in the frequency domain (mono or stereo)."""
    X = np.fft.rfft(x, axis=0); f = np.fft.rfftfreq(len(x), 1 / SR); g = np.ones_like(f)
    if lo: g *= 1 / np.sqrt(1 + (lo / np.maximum(f, 1e-9)) ** (2 * order))
    if hi: g *= 1 / np.sqrt(1 + (f / hi) ** (2 * order))
    return np.fft.irfft(X * (g[:, None] if x.ndim == 2 else g), len(x), axis=0)


# ---------------------------------------------------------------- arrangement
CHORDS = [  # two bars each: bass note, pad voicing, lead tones for small / medium / large branches (MIDI)
    (33, [45, 52, 59, 60, 64, 67], (69, 72, 76)),   # Am9
    (29, [41, 48, 57, 60, 64, 67], (69, 72, 77)),   # Fmaj9
    (36, [48, 55, 59, 62, 64, 67], (67, 72, 76)),   # Cmaj9
    (31, [43, 50, 57, 59, 62, 64], (67, 71, 74)),   # G6/9
]
chord_at = lambda bar: CHORDS[(int(bar) // 2) % 4] if bar < 62 else CHORDS[0]
PENTA = [m for m in range(45, 89) if (m - 45) % 12 in (0, 3, 5, 7, 10)]       # A2 .. E6

pad, bass, kickb, drums, lead, stabs, ticks, fx = Bus(), Bus(), Bus(), Bus(), Bus(), Bus(), Bus(), Bus()
FILL = (23, 31, 47, 55)          # last bar of each 8-bar phrase in a drop: kick drops out for a snare fill

kicks = []
for bar in range(NBARS):
    s = sect(bar)
    for q in range(4):
        t = bar * BAR + q * BEAT
        if s in DROPS:
            if not (bar in FILL and q == 3): kicks.append((t, 1.0))
        elif s == 'intro' and bar >= 4 and q == 0: kicks.append((t, 0.3))
        elif s == 'build' and bar < 12 and q in (0, 2): kicks.append((t, 0.38 if q == 0 else 0.25))
        elif (s == 'build' and 12 <= bar < 15) or (s == 'break' and 36 <= bar < 39):
            kicks.append((t, 0.45 + 0.15 * (bar % 4) + 0.03 * q))
        elif s == 'break' and bar < 36 and q == 0: kicks.append((t, 0.3))
        elif s == 'outro' and bar < 58: kicks.append((t, 0.85 - 0.2 * (bar - 56) - 0.04 * q))
        elif s == 'outro' and bar < 60 and q == 0: kicks.append((t, 0.3))
# sidechain: gain dips on every kick
duck = np.ones(NS)
for t, v in kicks:
    i = int(t * SR); n = int(0.36 * SR); x = np.arange(n) / n
    duck[i:i + n] = np.minimum(duck[i:i + n], 1 - 0.8 * min(1, v * 1.1) * (1 - x) ** 2.0)[:len(duck[i:i + n])]

# ---- pad: soft in the quiet sections, a bright pumping wall in the drops
PAD_LVL = {'intro': 0.9, 'build': 0.95, 'drop1': 1.5, 'break': 1.0, 'drop2': 1.7, 'outro': 1.0, 'end': 0.7}
PAD_BR = {'intro': 1.6, 'build': 2.6, 'drop1': 7.0, 'break': 3.0, 'drop2': 9.0, 'outro': 3.0, 'end': 1.8}
for slot in range(0, NBARS, 2):
    s = sect(slot); _, notes, _ = chord_at(slot)
    dur = 2 * BAR + 1.3; n = int(dur * SR); t = np.arange(n) / SR
    e = np.minimum(1, t / 0.9) ** 1.5 * np.minimum(1, (dur - t) / 1.3) ** 1.5
    br = PAD_BR[s] * (1 + 0.8 * (slot - 8) / 8 if s == 'build' else 1)
    L = np.zeros(n); R = np.zeros(n)
    for j, m in enumerate(notes):
        f0 = mtof(m + (12 if s in DROPS and j >= 3 else 0))
        for det, side in ((-6, 0), (6, 1), (0.7, j % 2)):
            f = f0 * 2 ** (det / 1200); ph = rng.uniform(0, 6.28); y = np.zeros(n)
            for k in range(1, 15):
                if f * k > 9000 or np.exp(-(k - 1) / br) < 0.01: break
                y += np.sin(2 * np.pi * f * k * t + ph * k) * np.exp(-(k - 1) / br) / k
            y *= 1 + 0.12 * np.sin(2 * np.pi * (0.13 + 0.05 * j) * t + j)
            if side == 0: L += y
            else: R += y
    st = np.stack([L * 0.8 + R * 0.2, R * 0.8 + L * 0.2], 1) * e[:, None] * 0.03 * PAD_LVL[s]
    pad.add2(slot * BAR - 0.4 if slot else 0, st)

# ---- bass: a held sub in the quiet sections; in the drops, offbeat notes with a quieter roll between
def sub_note(f, dur):
    n = int(dur * SR); t = np.arange(n) / SR
    y = np.sin(2 * np.pi * f * t) + 0.3 * np.sin(2 * np.pi * 2 * f * t)
    return np.tanh(1.2 * y * np.minimum(1, t / 0.4) * np.minimum(1, (dur - t) / 0.3)) * 0.75
def drop_bass(f, dur, vel):
    n = int((dur + 0.03) * SR); t = np.arange(n) / SR
    fc = 170 + 1100 * vel * np.exp(-t / 0.07)
    k = np.arange(1, 15)[:, None]
    y = (np.sin(2 * np.pi * f * k * t[None]) / k / np.sqrt(1 + (k * f / fc[None]) ** 4)).sum(0) + 0.6 * np.sin(2 * np.pi * f * t)
    e = np.minimum(1, t / 0.004) * np.minimum(1, (dur + 0.03 - t) / 0.03)
    return np.tanh(1.6 * y * e) * vel
for bar in range(NBARS):
    s = sect(bar); b0 = chord_at(bar)[0]
    if s in DROPS:
        for i in range(16):
            if i % 4 == 0 or (bar in FILL and i >= 12): continue
            off = i % 4 == 2
            f = mtof(b0 + (12 if (i == 14 and bar % 2 == 1) else 0))
            bass.add(bar * BAR + i * SIXT, drop_bass(f, SIXT * (1.7 if off else 0.8), 1.0 if off else 0.5), 0, 0.2)
    elif bar % 2 == 0 and bar >= 2:
        lvl = {'intro': 0.1, 'build': 0.14, 'break': 0.14, 'outro': 0.12, 'end': 0.08}[s]
        bass.add(bar * BAR, sub_note(mtof(b0 + 12 if b0 < 31 else b0), (2 if s != 'end' else 1.5) * BAR), 0, lvl)

# ---- drums
def kick(v):
    n = int(0.42 * SR); t = np.arange(n) / SR
    ph = 2 * np.pi * (54 * t + 190 * 0.016 * (1 - np.exp(-t / 0.016)) + 45 * 0.06 * (1 - np.exp(-t / 0.06)))
    body = np.sin(ph) * (0.12 + 0.88 * np.exp(-t / 0.1)) * np.clip((0.3 - t) / 0.1, 0, 1)
    click = filt(rng.standard_normal(n) * np.exp(-t / 0.003), 1500, 9000) * 0.5
    return np.tanh(2.2 * body * v) / 1.0 + click * v
def clap():
    n = int(0.35 * SR); t = np.arange(n) / SR
    e = sum(np.exp(-np.maximum(t - d, 0) / 0.012) * (t >= d) for d in (0, 0.011, 0.023)) * 0.5 + np.exp(-t / 0.09)
    return filt(rng.standard_normal(n) * e, 900, 5000) * 0.6
def snare():
    n = int(0.22 * SR); t = np.arange(n) / SR
    tone = np.sin(2 * np.pi * (185 * t + 70 * 0.02 * (1 - np.exp(-t / 0.02)))) * np.exp(-t / 0.05)
    return 0.5 * tone + filt(rng.standard_normal(n) * np.exp(-t / 0.07), 1400, 9500) * 0.9
def hat(d, lo=7000):
    n = int(d * 6 * SR); t = np.arange(n) / SR
    return filt(rng.standard_normal(n) * np.exp(-t / d), lo, None, 3)
HC = [hat(0.02) for _ in range(6)]; HO = [hat(0.075, 6000) for _ in range(4)]; RIDE = [hat(0.05, 9000) for _ in range(4)]
CL = [clap() for _ in range(4)]; SN = [snare() for _ in range(6)]
for t, v in kicks: kickb.add(t, kick(v), 0, 0.5)
for bar in range(NBARS):
    s = sect(bar); t0 = bar * BAR
    if s in DROPS:
        big = s == 'drop2'
        for q in (1, 3):
            if bar in FILL and q == 3: continue
            c = CL[(bar + q) % 4]; drums.add2(t0 + q * BEAT, np.stack([c, CL[(bar + q + 1) % 4]], 1), 0.2)
            drums.add(t0 + q * BEAT, SN[q], 0, 0.1)
        for i in range(16):
            if i % 4 == 2: drums.add(t0 + i * SIXT, HO[i % 4], 0.1, 0.2 if big else 0.17)
            else: drums.add(t0 + i * SIXT, HC[(i + bar) % 6], 0.3 * (-1) ** i, (0.05 + 0.03 * (i % 2)) * (1.2 if big else 1))
            if big: drums.add(t0 + i * SIXT, RIDE[i % 4], -0.2, 0.035 + 0.03 * (i % 4 == 0))
        if bar in FILL:
            for i in range(12, 16): drums.add(t0 + i * SIXT, SN[i % 6], 0.15 * (-1) ** i, 0.12 + 0.03 * (i - 12))
    elif (s == 'build' and bar < 12) or (s == 'break' and bar < 36) or (s == 'intro' and bar >= 6):
        for i in (4, 12): drums.add(t0 + i * SIXT, HC[i % 6], 0.2, 0.05)       # soft half-time backbeat tick
        if s != 'intro': drums.add(t0 + 8 * SIXT, CL[bar % 4], 0, 0.07)
    elif s == 'outro' and bar < 58:
        for i in range(16):
            if i % 4 == 2: drums.add(t0 + i * SIXT, HO[i % 4], 0.1, 0.12 - 0.04 * (bar - 56))
        for q in (1, 3): drums.add(t0 + q * BEAT, CL[q], 0, 0.14 - 0.04 * (bar - 56))
# build-ups: four bars of snare roll that doubles its rate, a noise riser, a rising tone
for b0 in (12, 36):
    hits = [(b0 + bb, i) for bb in (0, 1) for i in range(0, 16, 4)] + [(b0 + 2, i) for i in range(0, 16, 2)] + [(b0 + 3, i) for i in range(0, 8)]
    hits = [(bb * BAR + i * SIXT) for bb, i in hits] + [(b0 + 3) * BAR + 8 * SIXT + j * SIXT / 2 for j in range(16)]
    T0, T1 = b0 * BAR, (b0 + 4) * BAR
    for j, t in enumerate(hits):
        x = (t - T0) / (T1 - T0)
        drums.add(t, SN[j % 6], 0.2 * (-1) ** j, 0.045 + 0.13 * x ** 1.5)
    n = int(4 * BAR * SR); x = np.arange(n) / n
    nz = rng.standard_normal((n, 2))
    lo = filt(nz, 300, 2500); hi_ = filt(nz, 2500, None)
    fx.add2(T0, (lo * (1 - x)[:, None] + hi_ * x[:, None] * 1.3) * (x ** 2.2)[:, None], 0.085)
    tr = np.arange(n) / SR; f = 220 * 2 ** (2 * x ** 1.5)                 # two octaves up, accelerating
    ph = 2 * np.pi * np.cumsum(f) / SR
    tone = (np.sin(ph) + 0.4 * np.sin(2 * ph) + 0.2 * np.sin(3 * ph)) * x ** 2 * (1 + 0.3 * np.sin(2 * np.pi * tr / SIXT))
    fx.add(T0, tone, 0, 0.035)
# impacts on the drops; a falling wash where each drop ends
for b in (16, 40):
    n = int(3.0 * SR); t = np.arange(n) / SR
    fx.add2(b * BAR, filt(rng.standard_normal((n, 2)), 3500, None) * np.exp(-t / 0.8)[:, None], 0.09)
    fx.add(b * BAR, np.sin(2 * np.pi * (38 * t + 40 * 0.25 * (1 - np.exp(-t / 0.25)))) * np.exp(-t / 0.6), 0, 0.3)
for b in (32, 56):
    n = int(2 * BAR * SR); x = np.arange(n) / n
    nz = rng.standard_normal((n, 2))
    fx.add2(b * BAR, (filt(nz, 2500, None) * (1 - x)[:, None] ** 2 + filt(nz, 300, 2500) * (x * (1 - x) * 3)[:, None]) * ((1 - x) ** 1.5)[:, None], 0.07)

# ---------------------------------------------------------------- what the chain plays
def mallet(f, m, long=1.0):
    d = (0.22 + 0.11 * np.log2(m + 1)) * long
    n = int(min(4.0, d * 6) * SR); t = np.arange(n) / SR
    idx = min(2.6, 2.4 * (220 / f) ** 0.35) * np.exp(-t / 0.11)
    y = np.sin(2 * np.pi * f * t + idx * np.sin(2 * np.pi * 2 * f * t)) + 0.25 * np.sin(2 * np.pi * 4.01 * f * t) * np.exp(-t / 0.05)
    return y * (1 - np.exp(-t / 0.002)) * np.exp(-t / d)
def pitch_of(m):
    k = int(round(len(PENTA) - 1 - np.log2(m) * (len(PENTA) - 1) / np.log2(1250)))
    return PENTA[max(0, min(len(PENTA) - 1, k))]
DET = np.array([-19, -12, -5, 0, 5, 12, 19]) / 1200.0
def stab(midi, bright, gate):
    """Seven detuned saws through a low-pass that snaps shut: one stereo trance stab."""
    f0 = mtof(midi); n = int((gate + 0.07) * SR); t = np.arange(n) / SR
    fc = 350 + bright * (0.45 + 0.55 * np.exp(-t / 0.06))
    k = np.arange(1, int(min(44, 15000 // f0)) + 1)[:, None]
    amp = 1 / k / np.sqrt(1 + (k * f0 / fc[None]) ** 4)
    L = np.zeros(n); R = np.zeros(n)
    for i, d in enumerate(DET):
        y = (np.sin(2 * np.pi * f0 * 2 ** d * k * t[None] + k * rng.uniform(0, 6.28)) * amp).sum(0)
        p = i / 6.0
        L += y * (1 - p) ** 0.5; R += y * p ** 0.5
    e = np.minimum(1, t / 0.003) * (0.7 + 0.3 * np.exp(-t / 0.05)) * np.clip((gate + 0.07 - t) / 0.07, 0, 1)
    return np.stack([L, R], 1) * e[:, None] / 7
def openness(b):
    """How far the lead's filter is open, by position in the song (bars)."""
    s = sect(b)
    if s in DROPS: return 1.0
    if s == 'build' and b >= 12: return 0.03 * (0.8 / 0.03) ** ((b - 12) / 4)
    if s == 'break' and b >= 36: return 0.03 * (0.8 / 0.03) ** ((b - 36) / 4)
    if s == 'outro' and b < 58: return 0.6 * (0.04 / 0.6) ** ((b - 56) / 2)
    return 0.0
nnotes = nstab = 0
for e in ev:
    f = e['f']; s = e['name']; k = e.get('k', -1)
    t = k * SIXT if k >= 0 else f / FPS
    drop = s in DROPS
    # relocations: up to 8 ticks per frame
    if e['rel']:
        nt = min(e['rel'], 8); amp = 0.085 * min(1.0, 2.2 / np.sqrt(nt)) * (1 + 0.5 * min(1, e['rel'] / 60))
        if drop: amp *= 0.5
        for _ in range(nt):
            fq = mtof(rng.choice(PENTA[-6:]) + 12 * rng.integers(1, 3)); n = int(0.05 * SR); x = np.arange(n) / SR
            ticks.add(f / FPS + rng.uniform(0, 1 / FPS), np.sin(2 * np.pi * fq * x) * np.exp(-x / 0.008) * (1 - np.exp(-x / 0.0004)), rng.uniform(-0.8, 0.8), amp * rng.uniform(0.5, 1))
    if k < 0 or not e['top']: continue
    b = k / 16.0; pos16 = k % 16
    # transplants: mallets
    for j, tr in enumerate(e['top']):
        if j == 1 and pos16 % 2: continue
        if j == 2 and pos16 % 4: continue
        if j >= 3 or (drop and j >= 1 and pos16 % 4): continue
        m = tr['m']; p = pitch_of(m)
        long = 1.9 if s in ('break', 'outro') else 1.0
        g = (0.042 + 0.008 * np.log2(m + 1)) * (1.0 if j == 0 else 0.55) * (1.2 if s == 'build' else 1.0) * (0.6 if drop else 1.0)
        lead.add(t, mallet(mtof(p), m, long), tr['pan'] * 0.8, g); nnotes += 1
        if m >= 350 and j == 0:
            n = int(1.6 * SR); x = np.arange(n) / SR
            if drop:      # zap: a fast falling chirp and a burst of hiss
                n = int(0.3 * SR); x = np.arange(n) / SR
                ph = 2 * np.pi * (300 * x + 3500 * 0.03 * (1 - np.exp(-x / 0.03)))
                fx.add(t, np.sin(ph) * np.exp(-x / 0.05), tr['pan'] * 0.7, 0.05 + 0.05 * min(1, m / 900))
                fx.add2(t, filt(rng.standard_normal((n, 2)), 3000, None) * np.exp(-x / 0.06)[:, None], 0.035)
            else:
                fx.add(t, np.sin(2 * np.pi * (41 * x + 35 * 0.12 * (1 - np.exp(-x / 0.12)))) * np.exp(-x / 0.5), 0, 0.2 * min(1, m / 900) + 0.08)
                fx.add2(t, filt(rng.standard_normal((n, 2)), 2000, None) * np.exp(-x / 0.35)[:, None], 0.03)
    # transplants: the lead synth
    op = openness(b)
    if op > 0:
        m = e['top'][0]['m']; lv = min(1.0, np.log2(m + 1) / 10)
        tone = chord_at(b)[2][0 if m < 24 else (1 if m < 150 else 2)]
        vel = (0.5 + 0.5 * lv) * (0.45 + 0.55 * op) * (1.15 if pos16 % 4 == 0 else 1.0)
        gate = SIXT * (0.72 if m < 150 else 1.5)
        y = stab(tone, (2500 + 9000 * lv ** 1.2) * op, gate)
        if s == 'drop2':
            y = y + 0.4 * stab(tone + 12, (2500 + 9000 * lv ** 1.2) * 0.8, gate)
        stabs.add2(t, y, vel); nstab += 1
print('mallet notes', nnotes, 'stabs', nstab, file=sys.stderr)

# ---------------------------------------------------------------- mix
def delay(x, time, fb, wet):
    """ping-pong delay"""
    d = int(time * SR); y = np.zeros_like(x); cur = x.copy()
    for k in range(1, 7):
        cur = np.roll(cur, d, axis=0); cur[:d] = 0; cur = cur[:, ::-1] * fb
        y += cur
    return x + wet * filt(y, 250, 5000)
def reverb(x, rt=2.4, wet=0.25):
    n = int(rt * 1.4 * SR); t = np.arange(n) / SR
    ir = rng.standard_normal((n, 2)) * np.exp(-t * 6.9 / rt)[:, None]; ir[:int(0.012 * SR)] = 0
    ir = filt(ir, 250, 6500); ir /= np.sqrt((ir ** 2).sum(0)).mean()
    N = 1 << int(np.ceil(np.log2(len(x) + n)))
    y = np.stack([np.fft.irfft(np.fft.rfft(x[:, c], N) * np.fft.rfft(ir[:, c], N), N)[:len(x)] for c in (0, 1)], 1)
    return x + wet * y

GAIN = dict(pad=1.2, bass=1.2, kick=0.95, drums=1.6, lead=1.0, stabs=1.0, ticks=1.0, fx=0.75)
if len(sys.argv) > 3: GAIN.update(json.loads(sys.argv[3]))
D = duck[:, None]
stems = {
    'pad': pad.out() * (0.12 + 0.88 * D),
    'bass': filt(bass.out(), None, 3000) * (0.08 + 0.92 * D),
    'kick': kickb.out(),
    'drums': drums.out(),
    'lead': delay(lead.out(), 3 * SIXT, 0.45, 0.33),
    'stabs': filt(delay(stabs.out(), 3 * SIXT, 0.4, 0.22), 150, None) * (0.5 + 0.5 * D),
    'ticks': ticks.out(),
    'fx': fx.out(),
}
for k in stems: stems[k] = stems[k] * GAIN[k]
wet = reverb(stems['lead'] + 0.8 * stems['ticks'] + 0.3 * stems['fx'] + 0.45 * stems['stabs'] + 0.25 * stems['drums'], 2.4, 1.0) \
    - (stems['lead'] + 0.8 * stems['ticks'] + 0.3 * stems['fx'] + 0.45 * stems['stabs'] + 0.25 * stems['drums'])
stems['reverb'] = 0.3 * wet * (0.6 + 0.4 * D)
mix = filt(sum(stems.values()), 28, None)

def db(x): return 20 * np.log10(np.sqrt((x ** 2).mean()) + 1e-9)
def seg(x, b0, b1): return x[int(b0 * BAR * SR):int(b1 * BAR * SR)]
print('stem levels (dB rms)   ' + ' '.join(f'{k:>6s}' for k in stems), file=sys.stderr)
for name, (b0, b1) in (('intro 4-8', (4, 8)), ('build 8-12', (8, 12)), ('build 12-16', (12, 16)), ('drop1', (16, 32)), ('break 32-36', (32, 36)), ('drop2', (40, 56)), ('outro 58-62', (58, 62))):
    print(f'{name:22s} ' + ' '.join(f'{db(seg(v, b0, b1)):6.1f}' for v in stems.values()), file=sys.stderr)

# master: set the level by the drops, soft-clip the peaks, normalise
ref = np.percentile(np.abs(seg(mix, 16, 32)), 99.5)
mix = np.tanh(mix / ref * 1.05)
fade = np.minimum(1, tt / 0.4) * np.minimum(1, np.maximum(0, (NS / SR - 0.15 - tt) / 2.5))
mix *= fade[:, None]
mix *= 0.89 / np.abs(mix).max()
for bar in range(0, NBARS, 2):
    s2 = seg(mix, bar, bar + 2); X = np.abs(np.fft.rfft(s2.mean(1))) ** 2; fr = np.fft.rfftfreq(len(s2), 1 / SR)
    bands = [10 * np.log10(X[(fr >= a) & (fr < b)].sum() / X.sum() + 1e-12) for a, b in ((20, 70), (70, 250), (250, 2000), (2000, 8000), (8000, 20000))]
    print(f'bars {bar:2d}-{bar + 2:2d} {sect(bar):6s} rms {db(s2):6.1f} dB  crest {20 * np.log10(np.abs(s2).max() / (np.sqrt((s2 ** 2).mean()) + 1e-9)):4.1f}  bands ' + ' '.join(f'{v:5.1f}' for v in bands), file=sys.stderr)
w = wave.open(sys.argv[2], 'wb'); w.setnchannels(2); w.setsampwidth(2); w.setframerate(SR)
w.writeframes((mix * 32767).astype('<i2').tobytes()); w.close()
