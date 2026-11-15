// The chain and its schedule for the music video.  Shared by node (dry run) and the renderer.
// The music is 138 BPM; the picture is 30 fps, so a sixteenth note is 3.26 frames and each grid
// point is drawn on the nearest frame (at most 17 ms from the sound).  Relocation proposals are
// spread over every frame; transplant proposals come in bursts on the grid, so branches jump on
// the beat.  Any fixed sequence of the two symmetric kernels keeps the uniform distribution
// stationary, so this is still the same chain.
(function (root) {
  const FPS = 30, BPM = 138, S16F = FPS * 15 / BPM, BARS = 64, NF = Math.round(BARS * 16 * S16F);
  const SIDE = 50, N = SIDE * SIDE, SEED = 20261114;
  const ACC = [4, 1, 2, 1, 3, 1, 2, 1, 3.5, 1, 2, 1, 3, 1, 2, 1.5];   // burst size by sixteenth
  const lerpExp = (a, b, t) => a * Math.pow(b / a, Math.max(0, Math.min(1, t)));
  const GRID = new Int32Array(NF).fill(-1);                           // frame -> sixteenth index, or -1
  for (let k = 0; k < BARS * 16; k++) { const f = Math.round(k * S16F); if (f < NF) GRID[f] = k; }
  // position of a frame in bars; a grid frame sits exactly on its sixteenth
  const barOf = f => GRID[f] >= 0 ? GRID[f] / 16 : f / (16 * S16F);

  // SECTIONS: [name, first bar, last bar)
  const SECTIONS = [['intro', 0, 8], ['build', 8, 16], ['drop1', 16, 32], ['break', 32, 40], ['drop2', 40, 56], ['outro', 56, 62], ['end', 62, 64]];
  function section(b) { for (const s of SECTIONS) if (b >= s[1] && b < s[2]) return s; return SECTIONS[SECTIONS.length - 1]; }

  // -> { R: relocation proposals this frame, T: transplant proposals this frame, k: sixteenth index or -1 }
  function schedule(f) {
    const b = barOf(f), [name, b0, b1] = section(b), u = (b - b0) / (b1 - b0), k = GRID[f];
    const p16 = k % 16, on16 = k >= 0, on8 = on16 && k % 2 === 0, on4 = on16 && k % 4 === 0, on2 = on16 && k % 8 === 0;
    let R = 0, T = 0;
    if (name === 'intro') { R = Math.round(lerpExp(3, 500, u)); }
    else if (name === 'build') {      // half notes, quarters, eighths, then sixteenths into the drop
      R = Math.round(lerpExp(500, 2500, u));
      const g = u < 0.25 ? on2 : u < 0.5 ? on4 : u < 0.75 ? on8 : on16;
      if (g) T = Math.round(lerpExp(200, 1400, u) * (u >= 0.75 ? ACC[p16] / 2 : 1));
    }
    else if (name === 'drop1') { R = 3000; if (on16) T = Math.round(750 * ACC[p16]); }
    else if (name === 'break') {      // slow and large, then a second run-up
      R = u < 0.5 ? 40 : Math.round(lerpExp(40, 2500, (u - 0.5) * 2));
      const g = u < 0.25 ? on2 : u < 0.5 ? on4 : u < 0.75 ? on8 : on16;
      if (g) T = u < 0.5 ? (on2 ? 2000 : 900) : Math.round(lerpExp(500, 1600, (u - 0.5) * 2) * (u >= 0.75 ? ACC[p16] / 2 : 1));
    }
    else if (name === 'drop2') { R = 3500; if (on16) T = Math.round(1050 * ACC[p16]); }
    else if (name === 'outro') {
      R = Math.round(lerpExp(2000, 4, u));
      const g = u < 1 / 3 ? on8 : u < 2 / 3 ? on4 : on2;
      if (g) T = Math.round(lerpExp(1200, 150, u));
    }
    return { R, T, name, k };
  }

  function makeSim(Sampler) {
    const cells = [];
    for (let y = 0; y < SIDE; y++) for (let x = 0; x < SIDE; x++) cells.push([x + 400, y + 400]);
    const s = Sampler(cells, 10, SEED);
    const n = s.n;
    const tFrame = new Int32Array(n).fill(-1000), tSize = new Float32Array(n);   // last transplant: frame, log size
    const rFrame = new Int32Array(n).fill(-1000);                                // last relocation
    const pX = new Int32Array(n), pY = new Int32Array(n);
    let moves = 0;
    // Advance one frame.  Returns what happened, for the picture and for the score.
    function step(f) {
      const { R, T, name, k } = schedule(f);
      pX.set(s.X); pY.set(s.Y);
      let relAcc = 0;
      // relocations: note which cells moved by diffing afterwards
      for (let i = 0; i < R; i++) if (s.relocate()) relAcc++;
      if (relAcc) for (let i = 0; i < n; i++) if (s.X[i] !== pX[i] || s.Y[i] !== pY[i]) rFrame[i] = f;
      moves += R + T;
      let burst = null;
      if (T) {
        burst = [];
        let cx = 0, cy = 0; for (let i = 0; i < n; i++) { cx += s.X[i]; cy += s.Y[i]; } cx /= n; cy /= n;
        const rg = Math.sqrt(s.rg2());
        for (let i = 0; i < T; i++) {
          if (!s.transplant()) continue;
          const B = s.branch, m = B.length;
          let ox = 0, oy = 0, nx = 0, ny = 0;
          for (let k = 0; k < m; k++) { const j = B[k]; ox += pX[j]; oy += pY[j]; nx += s.X[j]; ny += s.Y[j]; tFrame[j] = f; tSize[j] = Math.log2(m + 1); }
          burst.push({ m, from: [ox / m, oy / m], to: [nx / m, ny / m], pan: Math.max(-1, Math.min(1, (nx / m - cx) / (2.2 * rg + 1))), cells: m >= 8 ? Array.from(B) : null });
        }
        burst.sort((a, b) => b.m - a.m);
      }
      return { f, k, name, R, T, relAcc, burst, moves };
    }
    return { s, n, step, tFrame, tSize, rFrame, pX, pY };
  }
  const api = { FPS, BPM, S16F, BARS, NF, SIDE, N, GRID, barOf, SECTIONS, section, schedule, makeSim };
  if (typeof module !== 'undefined') module.exports = api; else root.SIM = api;
})(this);
