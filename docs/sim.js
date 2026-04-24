const els = {
  canvas: document.getElementById('simCanvas'),
  flavourCanvas: document.getElementById('flavourCanvas'),
  playPauseBtn: document.getElementById('playPauseBtn'),
  resetBtn: document.getElementById('resetBtn'),
  defaultPresetBtn: document.getElementById('defaultPresetBtn'),
  nearThresholdPresetBtn: document.getElementById('nearThresholdPresetBtn'),
  g: document.getElementById('g'),
  gNumber: document.getElementById('gNumber'),
  kappa1: document.getElementById('kappa1'),
  kappa1Number: document.getElementById('kappa1Number'),
  kappa2: document.getElementById('kappa2'),
  kappa2Number: document.getElementById('kappa2Number'),
  L: document.getElementById('L'),
  LNumber: document.getElementById('LNumber'),
  N: document.getElementById('N'),
  NNumber: document.getElementById('NNumber'),
  dt: document.getElementById('dt'),
  dtNumber: document.getElementById('dtNumber'),
  A: document.getElementById('A'),
  ANumber: document.getElementById('ANumber'),
  sigma: document.getElementById('sigma'),
  sigmaNumber: document.getElementById('sigmaNumber'),
  substeps: document.getElementById('substeps'),
  substepsNumber: document.getElementById('substepsNumber'),
  viewMode: document.getElementById('viewMode'),
  timeReadout: document.getElementById('timeReadout'),
  energyReadout: document.getElementById('energyReadout'),
  dxReadout: document.getElementById('dxReadout'),
  phiCReadout: document.getElementById('phiCReadout'),
  phiSReadout: document.getElementById('phiSReadout'),
  consistencyReadout: document.getElementById('consistencyReadout'),
  messageReadout: document.getElementById('messageReadout')
};

const ctx = els.canvas.getContext('2d');
const fctx = els.flavourCanvas.getContext('2d');
const ALPHA = Math.sqrt(2 * Math.PI);
const DELTA_CS = Math.sqrt(2 * Math.PI);
const SQRT2 = Math.sqrt(2);
const SQRT_PI = Math.sqrt(Math.PI);
const BETA = 2 * SQRT_PI;
let sim = null;
let running = false;
let rafId = null;

const pairs = [
  ['g', 'gNumber'],
  ['kappa1', 'kappa1Number'],
  ['kappa2', 'kappa2Number'],
  ['L', 'LNumber'],
  ['N', 'NNumber'],
  ['dt', 'dtNumber'],
  ['A', 'ANumber'],
  ['sigma', 'sigmaNumber'],
  ['substeps', 'substepsNumber']
];

function setPreset() {
  els.g.value = els.gNumber.value = '1.0';
  els.kappa1.value = els.kappa1Number.value = '5.0';
  els.kappa2.value = els.kappa2Number.value = '0.5';
  els.L.value = els.LNumber.value = '200';
  els.N.value = els.NNumber.value = '2000';
  els.dt.value = els.dtNumber.value = '0.02';
  els.A.value = els.ANumber.value = '4.0';
  els.sigma.value = els.sigmaNumber.value = '3.0';
  els.substeps.value = els.substepsNumber.value = '2';
  els.viewMode.value = 'both';
}

function setNearThresholdPreset() {
  els.g.value = els.gNumber.value = '1.58';
  els.kappa1.value = els.kappa1Number.value = '1.0';
  els.kappa2.value = els.kappa2Number.value = '1.0';
  els.L.value = els.LNumber.value = '240';
  els.N.value = els.NNumber.value = '2800';
  els.dt.value = els.dtNumber.value = '0.008';
  els.A.value = els.ANumber.value = '2.5';
  els.sigma.value = els.sigmaNumber.value = '4.0';
  els.substeps.value = els.substepsNumber.value = '2';
  els.viewMode.value = 'both';
}

function params() {
  return {
    g: Number(els.gNumber.value),
    kappa1: Number(els.kappa1Number.value),
    kappa2: Number(els.kappa2Number.value),
    L: Number(els.LNumber.value),
    N: Math.round(Number(els.NNumber.value)),
    dt: Number(els.dtNumber.value),
    A: Number(els.ANumber.value),
    sigma: Number(els.sigmaNumber.value),
    substeps: Math.round(Number(els.substepsNumber.value)),
    viewMode: els.viewMode.value
  };
}

function syncPair(from, to) {
  els[to].value = els[from].value;
}

function bindPair(rangeId, numberId) {
  els[rangeId].addEventListener('input', () => {
    syncPair(rangeId, numberId);
    reinitialize();
  });
  els[numberId].addEventListener('input', () => {
    syncPair(numberId, rangeId);
    reinitialize();
  });
}

function reinitialize() {
  running = false;
  els.playPauseBtn.textContent = 'Play';
  initialize();
}

function initialize() {
  const p = params();
  const dx = p.L / p.N;
  const x = new Float64Array(p.N);
  const phiC = new Float64Array(p.N);
  const phiS = new Float64Array(p.N);
  const piC = new Float64Array(p.N);
  const piS = new Float64Array(p.N);
  const phi1 = new Float64Array(p.N);
  const phi2 = new Float64Array(p.N);
  const pi1 = new Float64Array(p.N);
  const pi2 = new Float64Array(p.N);

  for (let i = 0; i < p.N; i++) {
    x[i] = -p.L / 2 + i * dx;
    const kick = p.A * Math.exp(-((x[i] / p.sigma) ** 2));
    piC[i] = kick;
    piS[i] = kick;
    pi1[i] = SQRT2 * kick;
    pi2[i] = 0;
  }

  sim = {
    ...p,
    dx,
    x,
    phiC,
    phiS,
    piC,
    piS,
    phi1,
    phi2,
    pi1,
    pi2,
    t: 0
  };

  updateReadout('initialized');
  draw();
  drawFlavourCheck();
}

function laplacianPeriodic(field, dx) {
  const n = field.length;
  const out = new Float64Array(n);
  const inv = 1 / (dx * dx);
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const im = i === 0 ? n - 1 : i - 1;
    out[i] = (field[ip] - 2 * field[i] + field[im]) * inv;
  }
  return out;
}

function accelerationsCS(phiC, phiS) {
  const n = phiC.length;
  const lapC = laplacianPeriodic(phiC, sim.dx);
  const lapS = laplacianPeriodic(phiS, sim.dx);
  const aC = new Float64Array(n);
  const aS = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const arg1 = ALPHA * (phiC[i] + phiS[i]);
    const arg2 = ALPHA * (phiC[i] - phiS[i]);
    aC[i] = lapC[i] - sim.g * sim.g * phiC[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1) - 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2);
    aS[i] = lapS[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(arg1) + 0.5 * ALPHA * sim.kappa2 * Math.sin(arg2);
  }
  return { aC, aS };
}

function accelerations12(phi1, phi2) {
  const n = phi1.length;
  const lap1 = laplacianPeriodic(phi1, sim.dx);
  const lap2 = laplacianPeriodic(phi2, sim.dx);
  const a1 = new Float64Array(n);
  const a2 = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const chargeTerm = 0.5 * sim.g * sim.g * (phi1[i] + phi2[i]);
    a1[i] = lap1[i] - chargeTerm - SQRT_PI * sim.kappa1 * Math.sin(BETA * phi1[i]);
    a2[i] = lap2[i] - chargeTerm - SQRT_PI * sim.kappa2 * Math.sin(BETA * phi2[i]);
  }
  return { a1, a2 };
}

function stepLeapfrog() {
  const n = sim.phiC.length;
  const dt = sim.dt;
  const { aC, aS } = accelerationsCS(sim.phiC, sim.phiS);
  const { a1, a2 } = accelerations12(sim.phi1, sim.phi2);

  const piCHalf = new Float64Array(n);
  const piSHalf = new Float64Array(n);
  const pi1Half = new Float64Array(n);
  const pi2Half = new Float64Array(n);
  const phiCNew = new Float64Array(n);
  const phiSNew = new Float64Array(n);
  const phi1New = new Float64Array(n);
  const phi2New = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    piCHalf[i] = sim.piC[i] + 0.5 * dt * aC[i];
    piSHalf[i] = sim.piS[i] + 0.5 * dt * aS[i];
    phiCNew[i] = sim.phiC[i] + dt * piCHalf[i];
    phiSNew[i] = sim.phiS[i] + dt * piSHalf[i];

    pi1Half[i] = sim.pi1[i] + 0.5 * dt * a1[i];
    pi2Half[i] = sim.pi2[i] + 0.5 * dt * a2[i];
    phi1New[i] = sim.phi1[i] + dt * pi1Half[i];
    phi2New[i] = sim.phi2[i] + dt * pi2Half[i];
  }

  const { aC: aCNew, aS: aSNew } = accelerationsCS(phiCNew, phiSNew);
  const { a1: a1New, a2: a2New } = accelerations12(phi1New, phi2New);
  const piCNew = new Float64Array(n);
  const piSNew = new Float64Array(n);
  const pi1New = new Float64Array(n);
  const pi2New = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    piCNew[i] = piCHalf[i] + 0.5 * dt * aCNew[i];
    piSNew[i] = piSHalf[i] + 0.5 * dt * aSNew[i];
    pi1New[i] = pi1Half[i] + 0.5 * dt * a1New[i];
    pi2New[i] = pi2Half[i] + 0.5 * dt * a2New[i];
  }

  sim.phiC = phiCNew;
  sim.phiS = phiSNew;
  sim.piC = piCNew;
  sim.piS = piSNew;
  sim.phi1 = phi1New;
  sim.phi2 = phi2New;
  sim.pi1 = pi1New;
  sim.pi2 = pi2New;
  sim.t += dt;
}

function maxAbs(arr) {
  let m = 0;
  for (let i = 0; i < arr.length; i++) m = Math.max(m, Math.abs(arr[i]));
  return m;
}

function minValue(arr) {
  let m = Infinity;
  for (let i = 0; i < arr.length; i++) m = Math.min(m, arr[i]);
  return m;
}

function maxValue(arr) {
  let m = -Infinity;
  for (let i = 0; i < arr.length; i++) m = Math.max(m, arr[i]);
  return m;
}

function convertedFields() {
  const n = sim.phi1.length;
  const phiCFrom12 = new Float64Array(n);
  const phiSFrom12 = new Float64Array(n);
  const diffC = new Float64Array(n);
  const diffS = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    phiCFrom12[i] = (sim.phi1[i] + sim.phi2[i]) / SQRT2;
    phiSFrom12[i] = (sim.phi1[i] - sim.phi2[i]) / SQRT2;
    diffC[i] = phiCFrom12[i] - sim.phiC[i];
    diffS[i] = phiSFrom12[i] - sim.phiS[i];
  }
  return { phiCFrom12, phiSFrom12, diffC, diffS };
}

function consistencyError() {
  const { diffC, diffS } = convertedFields();
  return Math.max(maxAbs(diffC), maxAbs(diffS));
}

function energyDensity() {
  const n = sim.phiC.length;
  const ed = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const im = i === 0 ? n - 1 : i - 1;
    const gradC = (sim.phiC[ip] - sim.phiC[im]) / (2 * sim.dx);
    const gradS = (sim.phiS[ip] - sim.phiS[im]) / (2 * sim.dx);
    const arg1 = ALPHA * (sim.phiC[i] + sim.phiS[i]);
    const arg2 = ALPHA * (sim.phiC[i] - sim.phiS[i]);
    const potential = 0.5 * sim.g * sim.g * sim.phiC[i] * sim.phiC[i] - 0.5 * sim.kappa1 * Math.cos(arg1) - 0.5 * sim.kappa2 * Math.cos(arg2);
    ed[i] = 0.5 * sim.piC[i] * sim.piC[i] + 0.5 * sim.piS[i] * sim.piS[i] + 0.5 * gradC * gradC + 0.5 * gradS * gradS + potential;
  }
  return ed;
}

function totalEnergy(ed) {
  let sum = 0;
  for (let i = 0; i < ed.length; i++) sum += ed[i];
  return sum * sim.dx;
}

function updateReadout(message) {
  const ed = energyDensity();
  els.timeReadout.textContent = sim.t.toFixed(2);
  els.energyReadout.textContent = totalEnergy(ed).toFixed(4);
  els.dxReadout.textContent = sim.dx.toFixed(4);
  els.phiCReadout.textContent = (maxAbs(sim.phiC) / DELTA_CS).toFixed(3);
  els.phiSReadout.textContent = (maxAbs(sim.phiS) / DELTA_CS).toFixed(3);
  els.consistencyReadout.textContent = consistencyError().toExponential(3);
  els.messageReadout.textContent = message;
}

function formatAxisValue(v) {
  const av = Math.abs(v);
  if (av >= 100) return v.toFixed(0);
  if (av >= 10) return v.toFixed(1);
  if (av >= 1e-2) return v.toFixed(2);
  return v.toExponential(1);
}

function boundsForArrays(arrays, sectorStep = null) {
  let yMinRaw = Infinity;
  let yMaxRaw = -Infinity;
  arrays.forEach(arr => {
    yMinRaw = Math.min(yMinRaw, minValue(arr));
    yMaxRaw = Math.max(yMaxRaw, maxValue(arr));
  });

  if (Math.abs(yMaxRaw - yMinRaw) < 1e-12) {
    yMinRaw -= 1;
    yMaxRaw += 1;
  }

  if (sectorStep === null) {
    const pad = 0.18 * (yMaxRaw - yMinRaw);
    return { ymin: yMinRaw - pad, ymax: yMaxRaw + pad, sectorStep: null, nLow: 1, nHigh: 0 };
  }

  let nLow = Math.ceil(yMinRaw / sectorStep);
  let nHigh = Math.floor(yMaxRaw / sectorStep);
  if (nLow > nHigh) {
    const nearest = Math.round(0.5 * (yMinRaw + yMaxRaw) / sectorStep);
    nLow = nearest;
    nHigh = nearest;
  }
  return {
    ymin: (nLow - 1) * sectorStep,
    ymax: (nHigh + 1) * sectorStep,
    sectorStep,
    nLow,
    nHigh
  };
}

function drawAxesOn(c, x0, y0, w, h, bounds, title) {
  c.strokeStyle = '#cbd5e1';
  c.lineWidth = 1;
  c.strokeRect(x0, y0, w, h);
  c.fillStyle = '#0f172a';
  c.font = '20px sans-serif';
  c.fillText(title, x0, y0 - 12);
  c.fillStyle = '#475569';
  c.font = '15px sans-serif';
  c.fillText(formatAxisValue(bounds.ymax), x0 - 44, y0 + 14);
  c.fillText(formatAxisValue(bounds.ymin), x0 - 44, y0 + h);
  c.fillText((-sim.L / 2).toFixed(0), x0 - 4, y0 + h + 24);
  c.fillText((sim.L / 2).toFixed(0), x0 + w - 34, y0 + h + 24);

  if (bounds.sectorStep !== null) {
    for (let n = bounds.nLow; n <= bounds.nHigh; n++) {
      const yVal = n * bounds.sectorStep;
      const yy = y0 + h - ((yVal - bounds.ymin) / (bounds.ymax - bounds.ymin)) * h;
      c.setLineDash([7, 7]);
      c.strokeStyle = '#94a3b8';
      c.beginPath();
      c.moveTo(x0, yy);
      c.lineTo(x0 + w, yy);
      c.stroke();
      c.setLineDash([]);
    }
  }
}

function drawCurveOn(c, arr, x0, y0, w, h, ymin, ymax, color, lineWidth = 2.4, dash = []) {
  c.strokeStyle = color;
  c.lineWidth = lineWidth;
  c.setLineDash(dash);
  c.beginPath();
  for (let i = 0; i < arr.length; i++) {
    const xx = x0 + (i / (arr.length - 1)) * w;
    const yy = y0 + h - ((arr[i] - ymin) / (ymax - ymin)) * h;
    if (i === 0) c.moveTo(xx, yy);
    else c.lineTo(xx, yy);
  }
  c.stroke();
  c.setLineDash([]);
}

function drawLegendOn(c, entries, x, y) {
  c.font = '14px sans-serif';
  entries.forEach((entry, idx) => {
    const yy = y + idx * 20;
    c.strokeStyle = entry.color;
    c.lineWidth = 3;
    c.setLineDash(entry.dash || []);
    c.beginPath();
    c.moveTo(x, yy);
    c.lineTo(x + 20, yy);
    c.stroke();
    c.setLineDash([]);
    c.fillStyle = '#0f172a';
    c.fillText(entry.label, x + 28, yy + 5);
  });
}

function drawSinglePanel(c, arr, x0, y0, w, h, title, color, sectorStep) {
  const bounds = boundsForArrays([arr], sectorStep);
  drawAxesOn(c, x0, y0, w, h, bounds, title);
  drawCurveOn(c, arr, x0, y0, w, h, bounds.ymin, bounds.ymax, color);
}

function drawComparePanel(c, arrA, arrB, x0, y0, w, h, title, colorA, colorB, sectorStep, labelA, labelB) {
  const bounds = boundsForArrays([arrA, arrB], sectorStep);
  drawAxesOn(c, x0, y0, w, h, bounds, title);
  drawCurveOn(c, arrA, x0, y0, w, h, bounds.ymin, bounds.ymax, colorA);
  drawCurveOn(c, arrB, x0, y0, w, h, bounds.ymin, bounds.ymax, colorB, 2.2, [10, 5]);
  drawLegendOn(c, [
    { color: colorA, label: labelA },
    { color: colorB, label: labelB, dash: [10, 5] }
  ], x0 + w - 172, y0 + 22);
}

function drawResidualPanel(c, arr, x0, y0, w, h, title, color) {
  const bounds = boundsForArrays([arr], null);
  drawAxesOn(c, x0, y0, w, h, bounds, title);
  drawCurveOn(c, arr, x0, y0, w, h, bounds.ymin, bounds.ymax, color);
}

function draw() {
  ctx.clearRect(0, 0, els.canvas.width, els.canvas.height);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, els.canvas.width, els.canvas.height);

  const padX = 64;
  const padY = 52;
  const gapX = 38;
  const gapY = 52;
  const plotW = (els.canvas.width - 2 * padX - gapX) / 2;
  const plotH = (els.canvas.height - 2 * padY - gapY) / 2;

  if (sim.viewMode === 'both' || sim.viewMode === 'charge_only') {
    drawSinglePanel(ctx, sim.phiC, padX, padY, plotW, plotH, 'phi_c(x,t)', '#2563eb', DELTA_CS);
    drawSinglePanel(ctx, sim.phiS, padX + plotW + gapX, padY, plotW, plotH, 'phi_s(x,t)', '#dc2626', DELTA_CS);
  }

  if (sim.viewMode === 'both' || sim.viewMode === 'flavour_only') {
    const y0 = padY + plotH + gapY;
    drawSinglePanel(ctx, sim.phi1, padX, y0, plotW, plotH, 'phi_1(x,t)', '#7c3aed', SQRT_PI);
    drawSinglePanel(ctx, sim.phi2, padX + plotW + gapX, y0, plotW, plotH, 'phi_2(x,t)', '#059669', SQRT_PI);
  }
}

function drawFlavourCheck() {
  fctx.clearRect(0, 0, els.flavourCanvas.width, els.flavourCanvas.height);
  fctx.fillStyle = '#ffffff';
  fctx.fillRect(0, 0, els.flavourCanvas.width, els.flavourCanvas.height);

  const { phiCFrom12, phiSFrom12, diffC, diffS } = convertedFields();
  const padX = 64;
  const padY = 52;
  const gapX = 38;
  const gapY = 54;
  const plotW = (els.flavourCanvas.width - 2 * padX - gapX) / 2;
  const plotH = (els.flavourCanvas.height - 2 * padY - gapY) / 2;

  drawComparePanel(
    fctx,
    sim.phiC,
    phiCFrom12,
    padX,
    padY,
    plotW,
    plotH,
    'phi_c: direct vs reconstructed',
    '#2563eb',
    '#f97316',
    DELTA_CS,
    'direct',
    'reconstructed'
  );

  drawComparePanel(
    fctx,
    sim.phiS,
    phiSFrom12,
    padX + plotW + gapX,
    padY,
    plotW,
    plotH,
    'phi_s: direct vs reconstructed',
    '#dc2626',
    '#0ea5e9',
    DELTA_CS,
    'direct',
    'reconstructed'
  );

  const y0 = padY + plotH + gapY;
  drawResidualPanel(fctx, diffC, padX, y0, plotW, plotH, 'Delta phi_c', '#111827');
  drawResidualPanel(fctx, diffS, padX + plotW + gapX, y0, plotW, plotH, 'Delta phi_s', '#ef4444');
}

function tick() {
  if (running) {
    for (let k = 0; k < sim.substeps; k++) stepLeapfrog();
    updateReadout('running');
    draw();
    drawFlavourCheck();
  }
  rafId = requestAnimationFrame(tick);
}

els.playPauseBtn.addEventListener('click', () => {
  running = !running;
  els.playPauseBtn.textContent = running ? 'Pause' : 'Play';
  updateReadout(running ? 'running' : 'paused');
});

els.resetBtn.addEventListener('click', () => {
  reinitialize();
});

els.defaultPresetBtn.addEventListener('click', () => {
  setPreset();
  reinitialize();
});

els.nearThresholdPresetBtn.addEventListener('click', () => {
  setNearThresholdPreset();
  reinitialize();
});

pairs.forEach(([rangeId, numberId]) => bindPair(rangeId, numberId));
els.viewMode.addEventListener('input', reinitialize);

setPreset();
initialize();
tick();