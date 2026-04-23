const els = {
  canvas: document.getElementById('simCanvas'),
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
  messageReadout: document.getElementById('messageReadout')
};

const ctx = els.canvas.getContext('2d');
const ALPHA = Math.sqrt(2 * Math.PI);
const DELTA = Math.sqrt(2 * Math.PI);
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

function setPreset(name) {
  if (name === 'default') {
    els.g.value = els.gNumber.value = '1.0';
    els.kappa1.value = els.kappa1Number.value = '5.0';
    els.kappa2.value = els.kappa2Number.value = '0.5';
    els.L.value = els.LNumber.value = '200';
    els.N.value = els.NNumber.value = '2000';
    els.dt.value = els.dtNumber.value = '0.02';
    els.A.value = els.ANumber.value = '4.0';
    els.sigma.value = els.sigmaNumber.value = '4.0';
    els.substeps.value = els.substepsNumber.value = '2';
    els.viewMode.value = 'both';
  } else {
    els.g.value = els.gNumber.value = '1.0';
    els.kappa1.value = els.kappa1Number.value = '5.0';
    els.kappa2.value = els.kappa2Number.value = '0.5';
    els.L.value = els.LNumber.value = '200';
    els.N.value = els.NNumber.value = '2000';
    els.dt.value = els.dtNumber.value = '0.02';
    els.A.value = els.ANumber.value = '4.2';
    els.sigma.value = els.sigmaNumber.value = '4.0';
    els.substeps.value = els.substepsNumber.value = '2';
    els.viewMode.value = 'both';
  }
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

  for (let i = 0; i < p.N; i++) {
    x[i] = -p.L / 2 + i * dx;
    const kick = p.A * Math.exp(-((x[i] / p.sigma) ** 2));
    piC[i] = kick;
    piS[i] = kick;
  }

  sim = {
    ...p,
    dx,
    x,
    phiC,
    phiS,
    piC,
    piS,
    t: 0
  };

  updateReadout('initialized');
  draw();
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

function accelerations(phiC, phiS) {
  const n = phiC.length;
  const lapC = laplacianPeriodic(phiC, sim.dx);
  const lapS = laplacianPeriodic(phiS, sim.dx);
  const aC = new Float64Array(n);
  const aS = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const a1 = ALPHA * (phiC[i] + phiS[i]);
    const a2 = ALPHA * (phiC[i] - phiS[i]);
    aC[i] = lapC[i] - sim.g * sim.g * phiC[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(a1) - 0.5 * ALPHA * sim.kappa2 * Math.sin(a2);
    aS[i] = lapS[i] - 0.5 * ALPHA * sim.kappa1 * Math.sin(a1) + 0.5 * ALPHA * sim.kappa2 * Math.sin(a2);
  }
  return { aC, aS };
}

function stepLeapfrog() {
  const { phiC, phiS, piC, piS, dt } = sim;
  const n = phiC.length;
  const { aC, aS } = accelerations(phiC, phiS);
  const piCHalf = new Float64Array(n);
  const piSHalf = new Float64Array(n);
  const phiCNew = new Float64Array(n);
  const phiSNew = new Float64Array(n);

  for (let i = 0; i < n; i++) {
    piCHalf[i] = piC[i] + 0.5 * dt * aC[i];
    piSHalf[i] = piS[i] + 0.5 * dt * aS[i];
    phiCNew[i] = phiC[i] + dt * piCHalf[i];
    phiSNew[i] = phiS[i] + dt * piSHalf[i];
  }

  const { aC: aCNew, aS: aSNew } = accelerations(phiCNew, phiSNew);
  const piCNew = new Float64Array(n);
  const piSNew = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    piCNew[i] = piCHalf[i] + 0.5 * dt * aCNew[i];
    piSNew[i] = piSHalf[i] + 0.5 * dt * aSNew[i];
  }

  sim.phiC = phiCNew;
  sim.phiS = phiSNew;
  sim.piC = piCNew;
  sim.piS = piSNew;
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

function energyDensity() {
  const n = sim.phiC.length;
  const ed = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const ip = i + 1 === n ? 0 : i + 1;
    const im = i === 0 ? n - 1 : i - 1;
    const gradC = (sim.phiC[ip] - sim.phiC[im]) / (2 * sim.dx);
    const gradS = (sim.phiS[ip] - sim.phiS[im]) / (2 * sim.dx);
    const a1 = ALPHA * (sim.phiC[i] + sim.phiS[i]);
    const a2 = ALPHA * (sim.phiC[i] - sim.phiS[i]);
    const potential = 0.5 * sim.g * sim.g * sim.phiC[i] * sim.phiC[i] - 0.5 * sim.kappa1 * Math.cos(a1) - 0.5 * sim.kappa2 * Math.cos(a2);
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
  const currentAbsPhiC = maxAbs(sim.phiC);
  const currentAbsPhiS = maxAbs(sim.phiS);
  els.timeReadout.textContent = sim.t.toFixed(2);
  els.energyReadout.textContent = totalEnergy(ed).toFixed(4);
  els.dxReadout.textContent = sim.dx.toFixed(4);
  els.phiCReadout.textContent = (currentAbsPhiC / DELTA).toFixed(3);
  els.phiSReadout.textContent = (currentAbsPhiS / DELTA).toFixed(3);
  els.messageReadout.textContent = message;
}

function yBounds(mode) {
  const arr = mode === 'phi_c' ? sim.phiC : sim.phiS;
  const yminRaw = minValue(arr);
  const ymaxRaw = maxValue(arr);
  const span = Math.max(ymaxRaw - yminRaw, 0.25 * DELTA, 0.2);
  const pad = 0.12 * span + 0.05 * DELTA;
  return {
    ymin: yminRaw - pad,
    ymax: ymaxRaw + pad
  };
}

function drawAxes(x0, y0, w, h, ymin, ymax, title) {
  ctx.strokeStyle = '#cbd5e1';
  ctx.lineWidth = 1;
  ctx.strokeRect(x0, y0, w, h);
  ctx.fillStyle = '#0f172a';
  ctx.font = '14px sans-serif';
  ctx.fillText(title, x0, y0 - 10);
  ctx.fillStyle = '#475569';
  ctx.font = '12px sans-serif';
  ctx.fillText(ymax.toFixed(2), 8, y0 + 12);
  ctx.fillText(ymin.toFixed(2), 8, y0 + h);
  ctx.fillText((-sim.L / 2).toFixed(0), x0 - 4, y0 + h + 18);
  ctx.fillText((sim.L / 2).toFixed(0), x0 + w - 26, y0 + h + 18);

  const nMin = Math.ceil(ymin / DELTA);
  const nMax = Math.floor(ymax / DELTA);
  for (let n = nMin; n <= nMax; n++) {
    const yVal = n * DELTA;
    const yy = y0 + h - ((yVal - ymin) / (ymax - ymin)) * h;
    ctx.setLineDash([6, 6]);
    ctx.strokeStyle = '#94a3b8';
    ctx.beginPath();
    ctx.moveTo(x0, yy);
    ctx.lineTo(x0 + w, yy);
    ctx.stroke();
    ctx.setLineDash([]);
  }
}

function drawCurve(arr, x0, y0, w, h, ymin, ymax, color) {
  ctx.strokeStyle = color;
  ctx.lineWidth = 2;
  ctx.beginPath();
  for (let i = 0; i < arr.length; i++) {
    const xx = x0 + (i / (arr.length - 1)) * w;
    const yy = y0 + h - ((arr[i] - ymin) / (ymax - ymin)) * h;
    if (i === 0) ctx.moveTo(xx, yy);
    else ctx.lineTo(xx, yy);
  }
  ctx.stroke();
}

function drawLegend(entries, x, y) {
  ctx.font = '13px sans-serif';
  entries.forEach((entry, idx) => {
    ctx.strokeStyle = entry.color;
    ctx.lineWidth = 3;
    ctx.beginPath();
    ctx.moveTo(x, y + idx * 20);
    ctx.lineTo(x + 18, y + idx * 20);
    ctx.stroke();
    ctx.fillStyle = '#0f172a';
    ctx.fillText(entry.label, x + 26, y + 4 + idx * 20);
  });
}

function draw() {
  ctx.clearRect(0, 0, els.canvas.width, els.canvas.height);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(0, 0, els.canvas.width, els.canvas.height);

  const pad = 48;
  const gap = 28;
  const plotW = els.canvas.width - 2 * pad;
  const mode = sim.viewMode;
  const panels = mode === 'both' ? 2 : 1;
  const panelH = (els.canvas.height - 2 * pad - (panels - 1) * gap) / panels;

  if (mode === 'phi_c' || mode === 'both') {
    const boundsC = yBounds('phi_c');
    drawAxes(pad, pad, plotW, panelH, boundsC.ymin, boundsC.ymax, 'phi_c(x,t)');
    drawCurve(sim.phiC, pad, pad, plotW, panelH, boundsC.ymin, boundsC.ymax, '#2563eb');
  }

  if (mode === 'phi_s' || mode === 'both') {
    const y0 = mode === 'both' ? pad + panelH + gap : pad;
    const boundsS = yBounds('phi_s');
    drawAxes(pad, y0, plotW, panelH, boundsS.ymin, boundsS.ymax, 'phi_s(x,t)');
    drawCurve(sim.phiS, pad, y0, plotW, panelH, boundsS.ymin, boundsS.ymax, '#dc2626');
  }

  if (mode === 'both') {
    drawLegend([
      { color: '#2563eb', label: 'phi_c' },
      { color: '#dc2626', label: 'phi_s' }
    ], els.canvas.width - 160, 32);
  }
}

function tick() {
  if (running) {
    for (let k = 0; k < sim.substeps; k++) stepLeapfrog();
    updateReadout('running');
    draw();
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
  setPreset('default');
  reinitialize();
});

els.nearThresholdPresetBtn.addEventListener('click', () => {
  setPreset('near');
  reinitialize();
});

pairs.forEach(([rangeId, numberId]) => bindPair(rangeId, numberId));
els.viewMode.addEventListener('input', reinitialize);

setPreset('default');
initialize();
tick();