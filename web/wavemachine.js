/*
 * Virtual Wave Machine: browser version.
 *
 * Model, parameters and controls follow the original Easy Java Simulations applet
 * by Oleg Tchernyshyov (Johns Hopkins University), see ../original/_library/TorsionalOscillator.xml.
 *
 *   dθ_i/dt = ω_i
 *   dω_i/dt = -k0 θ_i - b ω_i + k (θ_{i-1} - θ_i) + k (θ_{i+1} - θ_i)  [+ τ0 cos φ for i = 0]
 *   dφ/dt   = ω_drive
 *   dω_drive/dt = sweep rate
 *
 * Plain JavaScript, no libraries.
 */
'use strict';

(() => {
  // ---- Model constants (from the original applet) ----
  const N = 12;        // number of rods
  const K = 1.0;       // coupling through the wire, s^2 = C/I
  const K0 = 0.25;     // gravity restoring term, omega0^2
  const D = 0.75;      // rod spacing along the wire (drawing units)
  const H = 0.05;      // maximum RK4 step (the original used RKF45 with dt = 0.1)
  const FMAX = 2.5;    // top of the frequency range
  const SAMPLE_DT = 0.5; // simulated time between points on the energy plot
  const MAX_POINTS = 40000;

  const DEFAULTS = { freq: 0.2, sweep: 0.0005, ampl: 0.125, damp: 0.025, speed: 12 };

  // ---- State: y = [θ_0..θ_{N-1}, ω_0..ω_{N-1}, φ, ω_drive] ----
  const M = 2 * N + 2;
  const PHI = 2 * N, FREQ = 2 * N + 1;
  const y = new Float64Array(M);
  const k1 = new Float64Array(M), k2 = new Float64Array(M), k3 = new Float64Array(M),
        k4 = new Float64Array(M), tmp = new Float64Array(M);

  const p = { ...DEFAULTS };
  let t = 0;
  let clamped = false;
  let driveOn = true;
  let running = true;
  let sampleClock = 0;

  // Energy plot data
  let ptsF = [], ptsE = [], ptsDir = [];
  let picked = null; // {f, e} chosen by clicking the plot

  function deriv(s, ds) {
    for (let i = 0; i < N; i++) ds[i] = s[N + i];
    for (let i = 0; i < N; i++) {
      let f = -K0 * s[i] - p.damp * s[N + i];
      if (i > 0) f += K * (s[i - 1] - s[i]);
      if (i < N - 1) f += K * (s[i + 1] - s[i]);
      if (i === 0 && driveOn) f += p.ampl * Math.cos(s[PHI]);
      ds[N + i] = f;
    }
    if (clamped) { ds[N - 1] = 0; ds[2 * N - 1] = 0; }
    ds[PHI] = s[FREQ];
    ds[FREQ] = p.sweep;
  }

  function rk4(h) {
    deriv(y, k1);
    for (let j = 0; j < M; j++) tmp[j] = y[j] + 0.5 * h * k1[j];
    deriv(tmp, k2);
    for (let j = 0; j < M; j++) tmp[j] = y[j] + 0.5 * h * k2[j];
    deriv(tmp, k3);
    for (let j = 0; j < M; j++) tmp[j] = y[j] + h * k3[j];
    deriv(tmp, k4);
    for (let j = 0; j < M; j++) y[j] += (h / 6) * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
    // Keep the drive frequency inside the plotted range.
    if (y[FREQ] < 0) y[FREQ] = 0;
    if (y[FREQ] > FMAX) y[FREQ] = FMAX;
    if (clamped) { y[N - 1] = 0; y[2 * N - 1] = 0; }
    t += h;
  }

  function energy() {
    let e = 0;
    for (let i = 0; i < N; i++) {
      const w = y[N + i], th = y[i];
      e += 0.5 * w * w + 0.5 * K0 * th * th;
      if (i > 0) { const dth = th - y[i - 1]; e += 0.5 * K * dth * dth; }
    }
    return e;
  }

  function advance(simTime) {
    const steps = Math.max(1, Math.ceil(simTime / H));
    const h = simTime / steps;
    for (let s = 0; s < steps; s++) {
      rk4(h);
      sampleClock += h;
      if (sampleClock >= SAMPLE_DT) {
        sampleClock -= SAMPLE_DT;
        addPoint(y[FREQ], energy(), Math.sign(p.sweep));
      }
    }
  }

  function addPoint(f, e, dir) {
    ptsF.push(f); ptsE.push(e); ptsDir.push(dir);
    if (ptsF.length > MAX_POINTS) { // thin out: keep every other point
      ptsF = ptsF.filter((_, i) => i % 2 === 0);
      ptsE = ptsE.filter((_, i) => i % 2 === 0);
      ptsDir = ptsDir.filter((_, i) => i % 2 === 0);
    }
  }

  function clearPlot() { ptsF = []; ptsE = []; ptsDir = []; picked = null; setPickedText(); }

  function resetAll() {
    Object.assign(p, DEFAULTS);
    y.fill(0);
    y[FREQ] = p.freq;
    t = 0; sampleClock = 0;
    clamped = false; driveOn = true;
    clearPlot();
    syncControls(true);
  }

  // Normal-mode frequencies of the 12-rod chain (exact for the discrete model).
  function predictedModes() {
    const out = [];
    if (!clamped) {
      // Both ends free: θ_i ∝ cos(q (i + 1/2)), q = π m / N.
      for (let m = 0; m < N; m++) {
        const s = Math.sin(Math.PI * m / (2 * N));
        out.push(Math.sqrt(K0 + 4 * K * s * s));
      }
    } else {
      // Left free, last rod clamped: N-1 moving rods, q = (2m+1) π / (2N - 1).
      for (let m = 0; m < N - 1; m++) {
        const s = Math.sin((2 * m + 1) * Math.PI / (2 * (2 * N - 1)));
        out.push(Math.sqrt(K0 + 4 * K * s * s));
      }
    }
    return out;
  }

  // ---- DOM ----
  const $ = (id) => document.getElementById(id);
  const machine = $('machine'), plot = $('plot');
  const mctx = machine.getContext('2d'), pctx = plot.getContext('2d');
  let view3d = true;
  let yaw = -0.6, elev = 0.32; // camera angles (radians)

  const sliders = [
    ['freq', 'freq-range', 'freq-num', 3],
    ['sweep', 'sweep-range', 'sweep-num', 5],
    ['ampl', 'ampl-range', 'ampl-num', 3],
    ['damp', 'damp-range', 'damp-num', 3],
    ['speed', 'speed-range', 'speed-num', 0],
  ];

  function fmt(v, digits) { return digits === 0 ? String(Math.round(v)) : v.toFixed(digits); }

  function setParam(key, v) {
    if (!Number.isFinite(v)) return;
    if (key === 'freq') { v = Math.min(FMAX, Math.max(0, v)); y[FREQ] = v; }
    if (key === 'speed') v = Math.max(1, v);
    if (key === 'damp' || key === 'ampl') v = Math.max(0, v);
    p[key] = v;
  }

  for (const [key, rid, nid, digits] of sliders) {
    const r = $(rid), n = $(nid);
    r.addEventListener('input', () => { setParam(key, parseFloat(r.value)); n.value = fmt(p[key], digits); });
    n.addEventListener('change', () => { setParam(key, parseFloat(n.value)); syncControls(true); });
  }

  // Write values into the controls. Fields being edited are left alone unless forced.
  function syncControls(force) {
    p.freq = y[FREQ];
    for (const [key, rid, nid, digits] of sliders) {
      const r = $(rid), n = $(nid);
      if (force || document.activeElement !== r) r.value = p[key];
      if (force || document.activeElement !== n) n.value = fmt(p[key], digits);
    }
    $('clamp').checked = clamped;
    $('drive-on').checked = driveOn;
    $('r-end').textContent = clamped ? 'clamped' : 'free';
  }

  function toggleClamp(v) {
    clamped = (v === undefined) ? !clamped : v;
    if (clamped) { y[N - 1] = 0; y[2 * N - 1] = 0; }
    syncControls(false);
  }

  $('clamp').addEventListener('change', (e) => toggleClamp(e.target.checked));
  $('drive-on').addEventListener('change', (e) => { driveOn = e.target.checked; });
  $('logy').addEventListener('change', drawPlot);
  $('theory').addEventListener('change', drawPlot);
  $('clear-plot').addEventListener('click', () => { clearPlot(); drawPlot(); });
  $('reset').addEventListener('click', () => { resetAll(); drawAll(); });
  $('kick').addEventListener('click', () => { y[0] = 0.6; y[N] = 0; drawAll(); });

  const playBtn = $('play');
  function setRunning(v) {
    running = v;
    playBtn.textContent = running ? 'Pause' : 'Play';
    last = performance.now();
  }
  playBtn.addEventListener('click', () => setRunning(!running));
  document.addEventListener('keydown', (e) => {
    if (e.code === 'Space' && !/INPUT|BUTTON|TEXTAREA|SELECT/.test(document.activeElement.tagName)) {
      e.preventDefault(); setRunning(!running);
    }
  });

  function setView(is3d) {
    view3d = is3d;
    $('view3d').setAttribute('aria-pressed', String(is3d));
    $('view2d').setAttribute('aria-pressed', String(!is3d));
    machine.style.cursor = is3d ? '' : 'pointer';
    $('machine-hint').textContent = is3d
      ? 'Drag to turn the view. The red rod is driven. Click the last rod to clamp or free that end.'
      : 'Rotation angle θ of each rod, in radians. The red rod is driven. Click the last rod to clamp or free that end.';
    drawMachine();
  }
  $('view3d').addEventListener('click', () => setView(true));
  $('view2d').addEventListener('click', () => setView(false));

  // ---- Canvas sizing (sharp on high-DPI screens) ----
  function fitCanvas(c, ctx) {
    const dpr = window.devicePixelRatio || 1;
    const w = c.clientWidth, h = c.clientHeight;
    if (c.width !== Math.round(w * dpr) || c.height !== Math.round(h * dpr)) {
      c.width = Math.round(w * dpr); c.height = Math.round(h * dpr);
    }
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { w, h };
  }

  let colors = {};
  function readColors() {
    const cs = getComputedStyle(document.documentElement);
    for (const k of ['fg', 'muted', 'line', 'grid', 'drive', 'rod', 'clamp', 'wire', 'up', 'down', 'theory', 'band', 'panel']) {
      colors[k] = cs.getPropertyValue('--' + k).trim();
    }
  }

  // ---- 3D wave machine ----
  // World: wire along Y, rods in the X–Z plane. θ is measured from horizontal (X).
  const LEN = (N - 1) * D;
  const CAM = 16;
  let proj = { S: 1, cx: 0, cy: 0 };

  function project(X, Y, Z) {
    const Yc = Y - LEN / 2;
    const ca = Math.cos(yaw), sa = Math.sin(yaw), ce = Math.cos(elev), se = Math.sin(elev);
    const u = Yc * ca + X * sa;
    const w = -Yc * sa + X * ca;
    const vert = Z * ce - w * se;
    const depth = w * ce + Z * se;
    const f = CAM / (CAM - depth);
    return { x: proj.cx + proj.S * f * u, y: proj.cy - proj.S * f * vert, f, depth };
  }

  function rodEnds(i) {
    const th = y[i], c = Math.cos(th), s = Math.sin(th), Y = i * D;
    return [project(-c, Y, -s), project(c, Y, s), project(0, Y, 0)];
  }

  function draw3D(w, h) {
    proj.S = Math.min(w / (LEN + 3.2), h / 4.2);
    proj.cx = w / 2; proj.cy = h / 2;
    const ctx = mctx;
    ctx.lineCap = 'round';

    // Wire, drawn a little past both end rods
    const a = project(0, -0.4, 0), b = project(0, LEN + 0.4, 0);
    ctx.strokeStyle = colors.wire; ctx.lineWidth = 1.5;
    ctx.beginPath(); ctx.moveTo(a.x, a.y); ctx.lineTo(b.x, b.y); ctx.stroke();

    // Rods, far to near
    const order = [...Array(N).keys()].map((i) => ({ i, e: rodEnds(i) }));
    order.sort((p1, p2) => p1.e[2].depth - p2.e[2].depth);
    for (const { i, e } of order) {
      const [p1, p2, pc] = e;
      const col = i === 0 ? colors.drive : (i === N - 1 && clamped ? colors.clamp : colors.rod);
      ctx.strokeStyle = col; ctx.fillStyle = col;
      ctx.lineWidth = 2.2 * pc.f;
      ctx.beginPath(); ctx.moveTo(p1.x, p1.y); ctx.lineTo(p2.x, p2.y); ctx.stroke();
      for (const q of [p1, p2]) {
        ctx.beginPath(); ctx.arc(q.x, q.y, 3.6 * q.f, 0, 2 * Math.PI); ctx.fill();
      }
      if (i === N - 1 && clamped) drawClamp(ctx, pc);
    }
  }

  function drawClamp(ctx, pc) {
    // A small bracket under the clamped rod's hub
    const r = 9 * pc.f;
    ctx.strokeStyle = colors.clamp; ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(pc.x - r, pc.y - r); ctx.lineTo(pc.x - r, pc.y + r); ctx.lineTo(pc.x + r, pc.y + r); ctx.lineTo(pc.x + r, pc.y - r);
    ctx.stroke();
    ctx.beginPath(); ctx.moveTo(pc.x, pc.y + r); ctx.lineTo(pc.x, pc.y + 2 * r); ctx.stroke();
  }

  // ---- 2D angle view ----
  let angleScale = 0.3;
  let view2dGeom = null;

  function draw2D(w, h) {
    const ctx = mctx;
    let peak = 0;
    for (let i = 0; i < N; i++) peak = Math.max(peak, Math.abs(y[i]));
    // Scale grows at once and shrinks slowly, so the picture doesn't jump around.
    angleScale = Math.max(0.25, peak * 1.15, angleScale * 0.995);
    const left = 44, right = 16, top = 12, bottom = 26;
    const pw = w - left - right, ph = h - top - bottom;
    const X = (i) => left + (i + 0.5) * pw / N;
    const Y = (th) => top + ph / 2 - (th / angleScale) * (ph / 2);
    view2dGeom = { X, Y };

    ctx.font = '12px system-ui, sans-serif';
    ctx.fillStyle = colors.muted; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    const step = niceStep(2 * angleScale);
    ctx.strokeStyle = colors.grid; ctx.lineWidth = 1;
    for (let v = -Math.floor(angleScale / step) * step; v <= angleScale + 1e-9; v += step) {
      const yy = Y(v);
      ctx.beginPath(); ctx.moveTo(left, yy); ctx.lineTo(w - right, yy); ctx.stroke();
      ctx.fillText(fmtTick(v, step), left - 6, yy);
    }
    ctx.strokeStyle = colors.muted;
    ctx.beginPath(); ctx.moveTo(left, Y(0)); ctx.lineTo(w - right, Y(0)); ctx.stroke();
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let i = 0; i < N; i++) ctx.fillText(String(i + 1), X(i), h - bottom + 8);

    ctx.strokeStyle = colors.wire; ctx.lineWidth = 1.5;
    ctx.beginPath();
    for (let i = 0; i < N; i++) (i ? ctx.lineTo : ctx.moveTo).call(ctx, X(i), Y(y[i]));
    ctx.stroke();
    for (let i = 0; i < N; i++) {
      ctx.fillStyle = i === 0 ? colors.drive : (i === N - 1 && clamped ? colors.clamp : colors.rod);
      ctx.beginPath(); ctx.arc(X(i), Y(y[i]), 6, 0, 2 * Math.PI); ctx.fill();
    }
  }

  function drawMachine() {
    const { w, h } = fitCanvas(machine, mctx);
    mctx.clearRect(0, 0, w, h);
    if (view3d) draw3D(w, h); else draw2D(w, h);
  }

  // Pointer: drag to rotate (3D), click the last rod to clamp/free.
  let drag = null;
  machine.addEventListener('pointerdown', (e) => {
    drag = { x: e.clientX, y: e.clientY, yaw, elev, moved: false };
    machine.setPointerCapture(e.pointerId);
  });
  machine.addEventListener('pointermove', (e) => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 4) drag.moved = true;
    if (drag.moved && view3d) {
      machine.classList.add('dragging');
      yaw = drag.yaw + dx * 0.01;
      elev = Math.max(-1.4, Math.min(1.4, drag.elev + dy * 0.01));
      if (!running) drawMachine();
    }
  });
  machine.addEventListener('pointerup', (e) => {
    machine.classList.remove('dragging');
    if (drag && !drag.moved) {
      const r = machine.getBoundingClientRect();
      if (hitsLastRod(e.clientX - r.left, e.clientY - r.top)) { toggleClamp(); drawMachine(); }
    }
    drag = null;
  });

  function hitsLastRod(px, py) {
    if (view3d) {
      const [p1, p2] = rodEnds(N - 1);
      return distToSegment(px, py, p1, p2) < 14;
    }
    if (!view2dGeom) return false;
    const { X, Y } = view2dGeom;
    return Math.hypot(px - X(N - 1), py - Y(y[N - 1])) < 16;
  }

  function distToSegment(px, py, a, b) {
    const vx = b.x - a.x, vy = b.y - a.y;
    const L2 = vx * vx + vy * vy || 1;
    const s = Math.max(0, Math.min(1, ((px - a.x) * vx + (py - a.y) * vy) / L2));
    return Math.hypot(px - (a.x + s * vx), py - (a.y + s * vy));
  }

  // ---- Energy vs frequency plot ----
  let plotGeom = null;

  function niceStep(range) {
    const raw = range / 5;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const n = raw / mag;
    return (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag;
  }
  function fmtTick(v, step) {
    const d = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    return (Math.abs(v) < step / 1e6 ? 0 : v).toFixed(d);
  }

  function drawPlot() {
    const { w, h } = fitCanvas(plot, pctx);
    const ctx = pctx;
    ctx.clearRect(0, 0, w, h);
    const log = $('logy').checked;
    const eNow = energy();

    let emax = eNow, emin = Infinity;
    for (let i = 0; i < ptsE.length; i++) {
      if (ptsE[i] > emax) emax = ptsE[i];
      if (ptsE[i] > 0 && ptsE[i] < emin) emin = ptsE[i];
    }
    emax = Math.max(emax, 1e-3);

    const left = 64, right = 14, top = 12, bottom = 40;
    const pw = w - left - right, ph = h - top - bottom;
    let y0, y1, Ymap;
    if (log) {
      y1 = Math.ceil(Math.log10(emax * 1.2));
      y0 = Math.floor(Math.log10(Math.max(Math.min(emin, emax), emax * 1e-6)));
      if (y1 - y0 < 1) y0 = y1 - 1;
      Ymap = (e) => top + ph * (1 - (Math.log10(Math.max(e, 1e-300)) - y0) / (y1 - y0));
    } else {
      y0 = 0; y1 = emax * 1.1;
      Ymap = (e) => top + ph * (1 - e / y1);
    }
    const Xmap = (f) => left + pw * f / FMAX;
    plotGeom = { left, top, pw, ph, log, y0, y1, Xmap };

    ctx.save();
    ctx.beginPath(); ctx.rect(left, top, pw, ph); ctx.clip();

    const showTheory = $('theory').checked;
    if (showTheory) {
      ctx.fillStyle = colors.band;
      ctx.fillRect(left, top, Xmap(Math.sqrt(K0)) - left, ph);
    }

    // Grid
    ctx.strokeStyle = colors.grid; ctx.lineWidth = 1;
    const fStep = 0.5;
    for (let f = 0; f <= FMAX + 1e-9; f += fStep) {
      ctx.beginPath(); ctx.moveTo(Xmap(f), top); ctx.lineTo(Xmap(f), top + ph); ctx.stroke();
    }
    let yTicks = [];
    if (log) { for (let d = y0; d <= y1; d++) yTicks.push(Math.pow(10, d)); }
    else { const s = niceStep(y1); for (let v = 0; v <= y1 + 1e-12; v += s) yTicks.push(v); }
    for (const v of yTicks) {
      ctx.beginPath(); ctx.moveTo(left, Ymap(v)); ctx.lineTo(left + pw, Ymap(v)); ctx.stroke();
    }

    if (showTheory) {
      ctx.strokeStyle = colors.theory; ctx.setLineDash([4, 4]); ctx.lineWidth = 1.2;
      for (const f of predictedModes()) {
        ctx.beginPath(); ctx.moveTo(Xmap(f), top); ctx.lineTo(Xmap(f), top + ph); ctx.stroke();
      }
      ctx.setLineDash([]);
    }

    // Data, as short segments colored by sweep direction
    ctx.lineWidth = 1.5;
    for (let i = 1; i < ptsF.length; i++) {
      if (ptsDir[i] !== ptsDir[i - 1] || Math.abs(ptsF[i] - ptsF[i - 1]) > 0.05) continue;
      ctx.strokeStyle = ptsDir[i] < 0 ? colors.down : colors.up;
      ctx.beginPath(); ctx.moveTo(Xmap(ptsF[i - 1]), Ymap(ptsE[i - 1])); ctx.lineTo(Xmap(ptsF[i]), Ymap(ptsE[i])); ctx.stroke();
    }
    for (let i = 0; i < ptsF.length; i++) { // isolated points (e.g. sweep rate 0)
      const lone = (i === 0 || ptsDir[i] !== ptsDir[i - 1]) && (i === ptsF.length - 1 || ptsDir[i + 1] !== ptsDir[i]);
      if (lone || ptsDir[i] === 0) {
        ctx.fillStyle = ptsDir[i] < 0 ? colors.down : colors.up;
        ctx.fillRect(Xmap(ptsF[i]) - 1.5, Ymap(ptsE[i]) - 1.5, 3, 3);
      }
    }

    // Current state
    ctx.fillStyle = colors.drive;
    ctx.beginPath(); ctx.arc(Xmap(y[FREQ]), Ymap(eNow), 4.5, 0, 2 * Math.PI); ctx.fill();

    if (picked) {
      ctx.strokeStyle = colors.fg; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
      const px = Xmap(picked.f), py = Ymap(picked.e);
      ctx.beginPath(); ctx.moveTo(px, top); ctx.lineTo(px, top + ph); ctx.moveTo(left, py); ctx.lineTo(left + pw, py); ctx.stroke();
      ctx.setLineDash([]);
    }
    ctx.restore();

    // Axes and labels
    ctx.strokeStyle = colors.muted; ctx.lineWidth = 1;
    ctx.strokeRect(left + 0.5, top + 0.5, pw, ph);
    ctx.fillStyle = colors.muted; ctx.font = '12px system-ui, sans-serif';
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (let f = 0; f <= FMAX + 1e-9; f += fStep) ctx.fillText(f.toFixed(1), Xmap(f), top + ph + 4);
    ctx.fillText('drive frequency ω (rad/s)', left + pw / 2, top + ph + 20);
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const v of yTicks) {
      const label = log ? (v >= 0.01 && v <= 100 ? String(+v.toPrecision(1)) : '1e' + Math.round(Math.log10(v))) : fmtTick(v, yTicks[1] - yTicks[0] || 1);
      ctx.fillText(label, left - 6, Ymap(v));
    }
    ctx.save();
    ctx.translate(12, top + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText('energy', 0, 0);
    ctx.restore();
  }

  plot.addEventListener('pointerdown', (e) => {
    if (!plotGeom) return;
    const r = plot.getBoundingClientRect();
    const px = e.clientX - r.left, py = e.clientY - r.top;
    const { left, top, pw, ph, log, y0, y1 } = plotGeom;
    if (px < left || px > left + pw || py < top || py > top + ph) return;
    const f = FMAX * (px - left) / pw;
    const frac = 1 - (py - top) / ph;
    const en = log ? Math.pow(10, y0 + frac * (y1 - y0)) : frac * y1;
    picked = { f, e: en };
    setPickedText();
    drawPlot();
  });

  function setPickedText() {
    $('plot-readout').textContent = picked
      ? `Marked: ω = ${picked.f.toFixed(3)} rad/s, energy = ${picked.e.toPrecision(3)}`
      : 'Click the plot to read off a frequency and energy.';
  }

  // ---- Main loop ----
  function updateReadouts() {
    $('r-t').textContent = t.toFixed(0);
    $('r-w').textContent = y[FREQ].toFixed(3);
    $('r-e').textContent = energy().toPrecision(3);
  }

  function drawAll() {
    readColors();
    drawMachine();
    drawPlot();
    updateReadouts();
    syncControls(false);
  }

  let last = performance.now();
  function frame(now) {
    const dt = Math.min(0.05, (now - last) / 1000);
    last = now;
    if (running) advance(p.speed * dt);
    drawAll();
    requestAnimationFrame(frame);
  }

  window.addEventListener('resize', drawAll);

  // Expose a small hook for testing from the console.
  window.waveMachine = { state: y, params: p, energy, predictedModes, advance, N };

  resetAll();
  requestAnimationFrame(frame);
})();
