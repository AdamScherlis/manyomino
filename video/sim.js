// The chain and its schedule for the music video.  Shared by node (dry run) and the renderer.
// Time is locked to a 90 BPM grid at 30 fps: a sixteenth note is 5 frames, a bar is 80.
// Relocation proposals are spread over every frame; transplant proposals come in bursts on
// the grid, so branches jump on the beat.  Any fixed sequence of the two symmetric kernels
// keeps the uniform distribution stationary, so this is still the same chain.
(function (root) {
  const FPS = 30, F16 = 5, BAR = 80, BARS = 34, NF = BARS * BAR;
  const SIDE = 50, N = SIDE * SIDE, SEED = 20261114;
  const ACC = [4, 1, 1.5, 1, 2.5, 1, 1.5, 1, 3, 1, 1.5, 1, 2.5, 1, 1.5, 2];   // burst size by sixteenth
  const lerpExp = (a, b, t) => a * Math.pow(b / a, Math.max(0, Math.min(1, t)));

  // SECTIONS: [name, first bar, last bar)
  const SECTIONS = [['intro', 0, 4], ['build', 4, 8], ['drop1', 8, 16], ['break', 16, 20], ['drop2', 20, 28], ['outro', 28, 32], ['end', 32, 34]];
  function section(f) { const b = f / BAR; for (const s of SECTIONS) if (b >= s[1] && b < s[2]) return s; return SECTIONS[SECTIONS.length - 1]; }

  // -> { R: relocation proposals this frame, T: transplant proposals this frame }
  function schedule(f) {
    const [name, b0, b1] = section(f), b = f / BAR, u = (b - b0) / (b1 - b0);
    const p16 = (f % BAR) / F16, on16 = f % F16 === 0, on8 = f % (2 * F16) === 0, on4 = f % (4 * F16) === 0;
    let R = 0, T = 0;
    if (name === 'intro') { R = Math.round(lerpExp(3, 500, u)); }
    else if (name === 'build') { R = Math.round(lerpExp(500, 2500, u)); const g = u < 0.5 ? on4 : on8; if (g) T = Math.round(lerpExp(150, 2500, u)); }
    else if (name === 'drop1') { R = 3000; if (on16) T = Math.round(900 * ACC[p16]); }
    else if (name === 'break') { R = 40; const g = u < 0.5 ? on4 : on8; if (g) T = on4 ? 1500 : 400; }
    else if (name === 'drop2') { R = 3500; if (on16) T = Math.round(1300 * ACC[p16]); }
    else if (name === 'outro') { R = Math.round(lerpExp(1500, 4, u)); const g = u < 0.5 ? on4 : (u < 0.75 ? f % (8 * F16) === 0 : f % BAR === 0); if (g) T = Math.round(lerpExp(900, 150, u)); }
    return { R, T, name };
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
      const { R, T, name } = schedule(f);
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
      return { f, name, R, T, relAcc, burst, moves };
    }
    return { s, n, step, tFrame, tSize, rFrame, pX, pY };
  }
  const api = { FPS, F16, BAR, BARS, NF, SIDE, N, SECTIONS, section, schedule, makeSim };
  if (typeof module !== 'undefined') module.exports = api; else root.SIM = api;
})(this);
