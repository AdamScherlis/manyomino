"""Soundtrack for the music video, synthesised from the chain's own event log.

  python3 video/synth.py events.json out.wav

90 BPM, locked to the video's frame grid (a sixteenth = 5 frames at 30 fps).  What the chain plays:
  * every accepted relocation is a tiny high tick;
  * on each transplant burst the largest accepted branch is a mallet note - the bigger the branch,
    the lower and longer the note (A minor pentatonic), panned to where the branch landed;
  * branches of 350+ cells also get a low boom.
Pad, bass and drums are a fixed arrangement around that.
"""
import json, sys, wave
import numpy as np

SR, FPS, F16, BARF = 44100, 30, 5, 80
ev = json.load(open(sys.argv[1]))
NF = len(ev); NS = int(NF / FPS * SR) + 1
rng = np.random.default_rng(7)
SIXT, BEAT, BAR = F16 / FPS, 4 * F16 / FPS, BARF / FPS
SECT = {'intro': (0, 4), 'build': (4, 8), 'drop1': (8, 16), 'break': (16, 20), 'drop2': (20, 28), 'outro': (28, 32), 'end': (32, 34)}
NBARS = NF // BARF
mtof = lambda m: 440.0 * 2 ** ((m - 69) / 12)
tt = np.arange(NS) / SR


def sect(bar):
    for k, (a, b) in SECT.items():
        if a <= bar < b: return k
    return 'end'


class Bus:
    def __init__(self): self.a = np.zeros((NS + SR * 6, 2))
    def add(self, t, mono, pan=0.0, gain=1.0):
        i = int(round(t * SR))
        if i < 0 or i >= NS: return
        th = (pan + 1) * np.pi / 4
        self.a[i:i + len(mono), 0] += mono * (gain * np.cos(th)); self.a[i:i + len(mono), 1] += mono * (gain * np.sin(th))
    def add2(self, t, st, gain=1.0):
        i = int(round(t * SR)); self.a[i:i + len(st)] += st * gain
    def out(self): return self.a[:NS]


def env_ad(n, a, d):
    t = np.arange(n) / SR
    return (1 - np.exp(-t / a)) * np.exp(-t / d)


def filt(x, lo=None, hi=None, order=2):
    """Zero-phase Butterworth-shaped filter in the frequency domain (mono or stereo)."""
    X = np.fft.rfft(x, axis=0); f = np.fft.rfftfreq(len(x), 1 / SR); g = np.ones_like(f)
    if lo: g *= 1 / np.sqrt(1 + (lo / np.maximum(f, 1e-9)) ** (2 * order))
    if hi: g *= 1 / np.sqrt(1 + (f / hi) ** (2 * order))
    return np.fft.irfft(X * (g[:, None] if x.ndim == 2 else g), len(x), axis=0)


# ---------------------------------------------------------------- arrangement
CHORDS = [  # two bars each: bass note, pad voicing (MIDI)
    (33, [45, 52, 59, 60, 64, 67]),   # Am9
    (29, [41, 48, 57, 60, 64, 67]),   # Fmaj9
    (36, [48, 55, 59, 62, 64, 67]),   # Cmaj9
    (28, [40, 47, 55, 59, 62, 66 - 2]),   # Em7
]
chord_at = lambda bar: CHORDS[(int(bar) // 2) % 4] if bar < 32 else CHORDS[0]
PENTA = [m for m in range(45, 89) if (m - 45) % 12 in (0, 3, 5, 7, 10)]       # A2 .. E6

pad, bass, drums, lead, ticks, fx = Bus(), Bus(), Bus(), Bus(), Bus(), Bus()

# sidechain: gain dips on every kick
kicks = []
for bar in range(NBARS):
    s = sect(bar)
    for q in range(4):
        t = bar * BAR + q * BEAT
        if s in ('drop1', 'drop2'): kicks.append((t, 1.0))
        elif s == 'build' and bar >= 6: kicks.append((t, 0.45 + 0.25 * (bar - 6) + 0.06 * q))
        elif s == 'outro' and bar < 30 and q % 2 == 0: kicks.append((t, 0.6 - 0.15 * (bar - 28)))
duck = np.ones(NS)
for t, v in kicks:
    i = int(t * SR); n = int(0.42 * SR); x = np.arange(n) / n
    duck[i:i + n] = np.minimum(duck[i:i + n], 1 - 0.62 * v * (1 - x) ** 2.2)

# ---- pad
PAD_LVL = {'intro': 0.9, 'build': 0.95, 'drop1': 1.0, 'break': 0.9, 'drop2': 1.05, 'outro': 0.9, 'end': 0.7}
PAD_BR = {'intro': 1.6, 'build': 2.6, 'drop1': 4.5, 'break': 3.0, 'drop2': 5.5, 'outro': 2.6, 'end': 1.8}
for slot in range(0, NBARS, 2):
    s = sect(slot); _, notes = chord_at(slot)
    dur = 2 * BAR + 1.6; n = int(dur * SR); t = np.arange(n) / SR
    e = np.minimum(1, t / 1.1) ** 1.5 * np.minimum(1, (dur - t) / 1.6) ** 1.5
    br = PAD_BR[s] * (1 + 0.5 * t / dur if s == 'build' else 1)
    L = np.zeros(n); R = np.zeros(n)
    for j, m in enumerate(notes):
        f0 = mtof(m)
        for det, side in ((-5, 0), (5, 1), (0.7, j % 2)):
            f = f0 * 2 ** (det / 1200); ph = rng.uniform(0, 6.28); y = np.zeros(n)
            for k in range(1, 9):
                if f * k > 9000: break
                y += np.sin(2 * np.pi * f * k * t + ph * k) * np.exp(-(k - 1) / br) / k
            y *= 1 + 0.12 * np.sin(2 * np.pi * (0.13 + 0.05 * j) * t + j)
            if side == 0: L += y
            else: R += y
    st = np.stack([L * 0.8 + R * 0.2, R * 0.8 + L * 0.2], 1) * e[:, None] * 0.03 * PAD_LVL[s]
    pad.add2(slot * BAR - 0.5 if slot else 0, st)

# ---- bass
BASS_PAT = [0, 3, 6, 8, 11, 14]
def bass_note(f, dur, punch=True):
    n = int(dur * SR); t = np.arange(n) / SR
    y = np.sin(2 * np.pi * f * t) + 0.45 * np.sin(2 * np.pi * 2 * f * t) * np.exp(-t / 0.25) + 0.2 * np.sin(2 * np.pi * 3 * f * t) * np.exp(-t / 0.12)
    e = np.minimum(1, t / (0.006 if punch else 0.4)) * np.minimum(1, (dur - t) / 0.05)
    if punch: e *= 0.55 + 0.45 * np.exp(-t / 0.2)
    return np.tanh(1.4 * y * e) * 0.75
for bar in range(NBARS):
    s = sect(bar); root = mtof(chord_at(bar)[0] + 12 if chord_at(bar)[0] < 31 else chord_at(bar)[0])
    if s in ('drop1', 'drop2'):
        for a, b in zip(BASS_PAT, BASS_PAT[1:] + [16]):
            f = root * (2 if (a == 14 and bar % 2 == 1) else 1)
            bass.add(bar * BAR + a * SIXT, bass_note(f, (b - a) * SIXT * 0.86), 0, 0.2)
    elif (s == 'intro' and bar >= 2) or s in ('build', 'break', 'outro', 'end'):
        lvl = {'intro': 0.1, 'build': 0.14, 'break': 0.14, 'outro': 0.12, 'end': 0.08}[s]
        if bar % 2 == 0 or s == 'end': bass.add(bar * BAR, bass_note(root, (2 if s != 'end' else 1) * BAR, False), 0, lvl)

# ---- drums
def kick(v):
    n = int(0.5 * SR); t = np.arange(n) / SR
    ph = 2 * np.pi * (52 * t + (160 - 52) * 0.03 * (1 - np.exp(-t / 0.03)))
    return (np.sin(ph) * np.exp(-t / 0.16) + 0.25 * rng.standard_normal(n) * np.exp(-t / 0.004)) * v
noise_long = rng.standard_normal(SR)
def clap():
    n = int(0.35 * SR); t = np.arange(n) / SR
    e = sum(np.exp(-np.maximum(t - d, 0) / 0.012) * (t >= d) for d in (0, 0.011, 0.023)) * 0.5 + np.exp(-t / 0.09)
    return filt(rng.standard_normal(n) * e, 900, 5000) * 0.6
def hat(d, lo=7000):
    n = int(d * 6 * SR); t = np.arange(n) / SR
    return filt(rng.standard_normal(n) * np.exp(-t / d), lo, None, 3)
for t, v in kicks: drums.add(t, kick(v), 0, 0.5)
for bar in range(NBARS):
    s = sect(bar); t0 = bar * BAR
    if s in ('drop1', 'drop2'):
        for q in (1, 3): drums.add2(t0 + q * BEAT, np.stack([clap(), clap()], 1), 0.21)
        for i in range(16):
            off = i % 4 == 2
            v = (0.17 if off else 0.055 + 0.03 * (i % 2)) * (1.0 if s == 'drop2' else 0.8)
            drums.add(t0 + i * SIXT, hat(0.07 if (off and s == 'drop2') else 0.022), 0.25 * (-1) ** i, v)
    if bar in (7, 19):        # roll into each drop
        for i in range(8, 16):
            drums.add(t0 + i * SIXT, clap(), 0, 0.04 + 0.012 * (i - 8))
            if i >= 12: drums.add(t0 + (i + 0.5) * SIXT, clap(), 0, 0.03 + 0.012 * (i - 12))
# risers into the drops, crash on the drops
for b0, b1 in ((6, 8), (18.5, 20)):
    n = int((b1 - b0) * BAR * SR); x = np.arange(n) / n
    nz = rng.standard_normal((n, 2))
    lo = filt(nz, 300, 2500); hi_ = filt(nz, 2500, None)
    fx.add2(b0 * BAR, (lo * (1 - x)[:, None] + hi_ * x[:, None] * 1.3) * (x ** 2.5)[:, None], 0.07)
for b in (8, 20):
    n = int(2.5 * SR); t = np.arange(n) / SR
    fx.add2(b * BAR, filt(rng.standard_normal((n, 2)), 3500, None) * np.exp(-t / 0.7)[:, None], 0.07)
    fx.add(b * BAR, np.sin(2 * np.pi * (38 * t + 30 * 0.25 * (1 - np.exp(-t / 0.25)))) * np.exp(-t / 0.7), 0, 0.3)

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
nnotes = 0
for e in ev:
    f = e['f']; t = f / FPS; s = e['name']; bar = f // BARF
    # relocations: up to 8 ticks per frame
    if e['rel']:
        k = min(e['rel'], 8); amp = 0.085 * min(1.0, 2.2 / np.sqrt(k)) * (1 + 0.5 * min(1, e['rel'] / 60))
        if s in ('drop1', 'drop2'): amp *= 0.55
        for _ in range(k):
            fq = mtof(rng.choice(PENTA[-6:]) + 12 * rng.integers(1, 3)); n = int(0.05 * SR); x = np.arange(n) / SR
            ticks.add(t + rng.uniform(0, 1 / FPS), np.sin(2 * np.pi * fq * x) * np.exp(-x / 0.008) * (1 - np.exp(-x / 0.0004)), rng.uniform(-0.8, 0.8), amp * rng.uniform(0.5, 1))
    # transplants
    pos16 = (f % BARF) // F16
    for j, tr in enumerate(e['top']):
        if j == 1 and pos16 % 2: continue
        if j == 2 and pos16 % 4: continue
        if j >= 3: continue
        m = tr['m']; p = pitch_of(m)
        long = 1.9 if s in ('break', 'outro') else 1.0
        g = (0.042 + 0.008 * np.log2(m + 1)) * (1.0 if j == 0 else 0.55) * (1.2 if s == 'build' else 1.0)
        lead.add(t, mallet(mtof(p), m, long), tr['pan'] * 0.8, g); nnotes += 1
        if s == 'drop2' and j == 0 and m >= 30: lead.add(t, mallet(mtof(p + 12), m, 0.6), -tr['pan'] * 0.6, g * 0.3)
        if m >= 350 and j == 0:
            n = int(1.6 * SR); x = np.arange(n) / SR
            fx.add(t, np.sin(2 * np.pi * (41 * x + 35 * 0.12 * (1 - np.exp(-x / 0.12)))) * np.exp(-x / 0.5), 0, 0.2 * min(1, m / 900) + 0.08)
            fx.add2(t, filt(rng.standard_normal((n, 2)), 2000, None) * np.exp(-x / 0.35)[:, None], 0.03)
print('notes', nnotes, file=sys.stderr)

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

D = duck[:, None]
L = delay(lead.out(), 3 * SIXT, 0.45, 0.33)
wetbus = reverb(L + 0.8 * ticks.out() + 0.25 * fx.out(), 2.6, 0.3)
mix = pad.out() * (0.35 + 0.65 * D) + filt(bass.out(), None, 2500) * (0.2 + 0.8 * D) + drums.out() + wetbus + 0.75 * fx.out()
mix = filt(mix, 28, None)
# level: gentle bus compression by section-agnostic soft clip, then normalise
mix = np.tanh(mix * 1.6) / 1.6
fade = np.minimum(1, tt / 0.4) * np.minimum(1, np.maximum(0, (NS / SR - 0.15 - tt) / 2.2))
mix *= fade[:, None]
mix *= 0.89 / np.abs(mix).max()
for name, b in (('pad', pad), ('bass', bass), ('drums', drums), ('lead', lead), ('ticks', ticks), ('fx', fx)):
    print(f'{name:6s} rms {np.sqrt((b.out() ** 2).mean()):.4f} peak {np.abs(b.out()).max():.3f}', file=sys.stderr)
for bar in range(0, NBARS, 2):
    seg = mix[int(bar * BAR * SR):int((bar + 2) * BAR * SR)]
    print(f'bars {bar:2d}-{bar + 2:2d} {sect(bar):6s} rms {20 * np.log10(np.sqrt((seg ** 2).mean()) + 1e-9):6.1f} dB', file=sys.stderr)
w = wave.open(sys.argv[2], 'wb'); w.setnchannels(2); w.setsampwidth(2); w.setframerate(SR)
w.writeframes((mix * 32767).astype('<i2').tobytes()); w.close()
