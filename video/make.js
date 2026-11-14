// node video/make.js OUTDIR            -> OUTDIR/silent.mp4
// node video/make.js OUTDIR 100,900    -> OUTDIR/f100.jpg, OUTDIR/f900.jpg (stills, for checking)
const { chromium } = require('playwright'), path = require('path'), fs = require('fs'), { spawn } = require('child_process');
(async () => {
  const out = process.argv[2], stills = process.argv[3] ? process.argv[3].split(',').map(Number) : null;
  fs.mkdirSync(out, { recursive: true });
  const browser = await chromium.launch(), page = await browser.newPage({ viewport: { width: 1920, height: 1080 } });
  page.on('pageerror', e => { console.error('page error', e); process.exit(1); });
  await page.goto('file://' + path.join(__dirname, 'render.html'), { timeout: 300000 });
  await page.evaluate(() => document.fonts.load('300 62px "Inter Display"').then(() => document.fonts.load('400 22px "JetBrains Mono"')));
  const NF = await page.evaluate(() => SIM.NF);
  let ff = null;
  if (!stills) ff = spawn('ffmpeg', ['-y', '-loglevel', 'error', '-f', 'image2pipe', '-framerate', '30', '-c:v', 'mjpeg', '-i', '-', '-c:v', 'libx264', '-preset', 'medium', '-crf', '17', '-pix_fmt', 'yuv420p', '-threads', '2', path.join(out, 'silent.mp4')], { stdio: ['pipe', 'inherit', 'inherit'] });
  const last = stills ? Math.max(...stills) : NF - 1, t0 = Date.now();
  for (let f = 0; f <= last; f++) {
    const want = !stills || stills.includes(f);
    const url = await page.evaluate(([f, d]) => frame(f, d), [f, want]);
    if (!want) continue;
    const buf = Buffer.from(url.slice(url.indexOf(',') + 1), 'base64');
    if (stills) fs.writeFileSync(path.join(out, `f${f}.jpg`), buf);
    else if (!ff.stdin.write(buf)) await new Promise(r => ff.stdin.once('drain', r));
    if (f % 300 === 0) console.error(`frame ${f}/${NF}  ${((Date.now() - t0) / 1000).toFixed(0)}s`);
  }
  if (ff) { ff.stdin.end(); await new Promise(r => ff.on('close', r)); }
  await browser.close();
})();
