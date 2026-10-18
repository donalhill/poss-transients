const tip = document.getElementById("tip");
const V = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt = d3.format(",");
function hov(sel, html) {
  sel.on("mousemove", (ev, d) => {
    tip.innerHTML = html(d); tip.hidden = false;
    const x = Math.min(ev.clientX + 14, window.innerWidth - tip.offsetWidth - 10);
    const y = Math.min(ev.clientY + 14, window.innerHeight - tip.offsetHeight - 10);
    tip.style.left = x + "px"; tip.style.top = y + "px";
  }).on("mouseleave", () => { tip.hidden = true; });
}
function svg(id) { const s = d3.select("#" + id); s.selectAll("*").remove(); return s; }
function axisStyle(g) { g.attr("class", "ax"); g.select(".domain").remove(); }

/* ---------- 1. sky map (equal-area disc centred on the north celestial pole, as seen from the ground) ---------- */
const DEG = Math.PI / 180, SHADOW = "#ff4a2e";
function galToEq(l, b) {
  const L = l * DEG, B = b * DEG, aG = 192.85948 * DEG, dG = 27.12825 * DEG, lN = 122.93192 * DEG;
  const dec = Math.asin(Math.sin(dG) * Math.sin(B) + Math.cos(dG) * Math.cos(B) * Math.cos(lN - L));
  const y = Math.cos(B) * Math.sin(lN - L), x = Math.cos(dG) * Math.sin(B) - Math.sin(dG) * Math.cos(B) * Math.cos(lN - L);
  return [((aG + Math.atan2(y, x)) / DEG + 360) % 360, dec / DEG];
}
function ring(pts) {   // close a ring and wind it so d3 fills the small side
  const c = pts.concat([pts[0]]);
  let g = { type: "Polygon", coordinates: [c] };
  if (d3.geoArea(g) > 2 * Math.PI) g = { type: "Polygon", coordinates: [c.slice().reverse()] };
  return g;
}
function drawSky(D, id, shadow, south) {   // shadow: pick out the nine plates the shadow crossed; south: add the 22 southern ones
  const s = svg(id), W = 760, H = 760, cx = W / 2, cy = H / 2, Rd = W / 2 - 44;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const plates = D.plates.map(p => Object.assign({}, p, { g: ring(p.poly) }));
  const sp = south ? D.south.plates.map(p => Object.assign({}, p, { g: ring(p.poly), south: true })) : [];
  const clip = d3.max(plates.concat(sp).flatMap(p => p.poly), v => 90 - v[1]) + 2.5;
  const proj = d3.geoAzimuthalEqualArea().reflectX(true).rotate([0, -90]).clipAngle(clip)
    .scale(Rd / (2 * Math.sin(clip * DEG / 2))).translate([cx, cy]).precision(0.3);
  const path = d3.geoPath(proj);
  const col = d3.scaleSequentialLog(d3.interpolateViridis).domain([5, 1500]).clamp(true);
  const plane = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => galToEq(l, 0)) };

  s.append("path").datum({ type: "Sphere" }).attr("d", path).attr("fill", V("--paper")).attr("stroke", V("--rule"));
  s.append("path").datum(d3.geoGraticule().step([30, 30]).extentMajor([[-180, -80], [180, 80.01]])()).attr("d", path)
    .attr("fill", "none").attr("stroke", V("--rule")).attr("stroke-width", 0.7);
  s.append("g").selectAll("path").data(D.plates_out.map(p => Object.assign({}, p, { g: ring(p.poly) }))).join("path").attr("class", "hov").attr("d", p => path(p.g))
    .attr("fill", V("--ink-3")).attr("fill-opacity", shadow ? 0.15 : 0.35).attr("stroke", V("--surface")).attr("stroke-width", 0.5)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients in the catalogue<br>Left out of VASCO’s published tests: its archive date is not on VASCO’s list of observing nights`));
  s.append("g").selectAll("path").data(plates.filter(p => !(shadow && p.circle))).join("path").attr("class", "hov").attr("d", p => path(p.g))
    .attr("fill", p => p.N > 0 ? col(Math.max(p.N, 5)) : V("--surface")).attr("fill-opacity", shadow ? 0.3 : 1).attr("stroke", V("--surface")).attr("stroke-width", 0.5)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients`));
  // plane of the Milky Way, drawn over the tiles: a soft band with a dashed centre line
  const mw = s.append("g").attr("pointer-events", "none");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-opacity", 0.38).attr("stroke-width", 16).attr("stroke-linecap", "round");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#111418").attr("stroke-opacity", 0.85).attr("stroke-width", 1.5).attr("stroke-dasharray", "5 4");
  const sh = shadow ? plates.filter(p => p.circle).concat(sp) : [];
  // the nine shadow plates outlined in black on a white halo; the shadow itself as a ring in a hue viridis never uses
  s.append("g").selectAll("path").data(sh).join("path").attr("pointer-events", "none").attr("d", p => path(p.g)).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-width", 5);
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "hov").attr("d", p => path(p.g))
    .attr("fill", p => p.N !== undefined ? col(Math.max(p.N, 5)) : V("--surface")).attr("stroke", "#111418").attr("stroke-width", 1.8)
    .call(sel => hov(sel, p => p.south
      ? `<b>Plate ${p.id}</b> · ${p.date}<br>` + (p.N !== undefined
        ? `${fmt(p.N)} transients found by our extraction<br>${Math.round(p.frac * 100)}% of the plate inside the shadow<br>${p.obs} in the shadowed part, against ${Math.round(p.exp)} expected`
        : `${Math.round(p.frac * 100)}% of the plate inside the shadow<br>still being processed`)
      : `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients on the plate<br>${Math.round(p.frac * 100)}% of the plate inside the shadow<br>${p.obs} transients in the shadowed part; a typical plate would have about ${Math.round(p.exp)}`));

  s.append("g").attr("pointer-events", "none").selectAll("path").data(sh).join("path").attr("d", p => path(ring(p.circle)))
    .attr("fill", "none").attr("stroke", SHADOW).attr("stroke-width", 2.4);
  const darkHalo = sel => sel.attr("class", "ann").style("font-weight", 600).style("fill", "#111418").style("paint-order", "stroke").style("stroke", "#ffffff").style("stroke-width", "3px");
  for (let h = 0; h < 24; h += 2) {   // hours of right ascension round the rim
    const [x, y] = proj([h * 15, 90 - clip + 0.5]), r = Math.hypot(x - cx, y - cy), R = Rd + 16;
    s.append("text").attr("x", cx + (x - cx) / r * R).attr("y", cy + (y - cy) / r * R + 4).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text(h + "h");
  }
  for (const d of (south ? [-30, 0, 30, 60] : [0, 30, 60])) {
    const [x, y] = proj([172, d]);
    s.append("text").attr("x", x).attr("y", y - 3).attr("text-anchor", "middle").attr("fill", V("--ink-3"))
      .style("paint-order", "stroke").style("stroke", V("--paper")).style("stroke-width", "3px").text((d > 0 ? "+" : "") + d + "°");
  }
  const m = proj(galToEq(100, 0)); s.append("text").attr("x", m[0] + 12).attr("y", m[1] - 10).call(darkHalo).text("Milky Way");
  const t = proj(galToEq(168, -42)); if (!shadow) s.append("text").attr("x", t[0]).attr("y", t[1] + 4).attr("text-anchor", "middle").call(darkHalo).text("Taurus–Perseus");
}

/* ---------- how many transients a plate holds ---------- */
function drawSpread(D) {
  const s = svg("c-spread"), W = 760, H = 300, m = { l: 48, r: 16, t: 26, b: 44 };
  const N = D.plates.map(p => p.N), mu = d3.mean(N), sd = Math.sqrt(mu), lo = mu - 2 * sd, hi = mu + 2 * sd;
  const bins = d3.range(0, 1000, 25).map(a => ({ lo: a, hi: a + 25, n: N.filter(v => v >= a && v < a + 25).length }));
  const over = { lo: 1025, hi: 1075, n: N.filter(v => v >= 1000).length, more: true };
  const inBand = N.filter(v => v >= lo && v <= hi).length;
  const x = d3.scaleLinear().domain([0, 1080]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 60]).range([H - m.b, m.t]);
  s.append("rect").attr("x", x(lo)).attr("width", x(hi) - x(lo)).attr("y", m.t - 14).attr("height", H - m.b - m.t + 14).attr("fill", V("--theirs-soft"));
  s.append("text").attr("class", "ann").attr("x", x(hi) + 8).attr("y", m.t - 2).text(`${inBand} of ${N.length} plates fall in this range`);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 200, 400, 600, 800]).tickFormat(d3.format(",")).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", x(1050)).attr("y", H - m.b + 15).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text("Over 1,000");
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Transients per plate");
  s.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Number of plates");
  s.append("g").selectAll("rect.b").data(bins.concat([over])).join("rect").attr("class", "hov").attr("x", b => x(b.lo) + 0.75).attr("width", b => x(b.hi) - x(b.lo) - 1.5)
    .attr("y", b => y(b.n)).attr("height", b => y(0) - y(b.n)).attr("fill", V("--ours"))
    .call(sel => hov(sel, b => b.more ? `${b.n} plates hold more than 1,000 transients, up to ${fmt(d3.max(N))}` : `${b.n} plates hold ${fmt(b.lo)} to ${fmt(b.hi - 1)} transients`));
}
function erf(x) {   // Abramowitz and Stegun 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x)), y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}

/* ---------- 2. two ideas of chance ---------- */
function drawChance(D) {
  const P = D.paste, s = svg("c-chance"), W = 760, H = 330, m = { l: 44, r: 16, t: 14, b: 40 };
  const x = d3.scaleLinear().domain([0, 1300]).range([m.l, W - m.r]);
  const binw = P.bins[1] - P.bins[0];
  const dens = P.counts.map((c, i) => ({ x0: P.bins[i], x1: P.bins[i + 1], p: c / P.n / binw }));
  const mu = P.area_expected, sd = P.binom_sd;
  const norm = d3.range(mu - 5 * sd, mu + 5 * sd, 0.5).map(v => [v, Math.exp(-0.5 * ((v - mu) / sd) ** 2) / (sd * Math.sqrt(2 * Math.PI))]);
  const y = d3.scaleLinear().domain([0, d3.max(norm, d => d[1]) * 1.08]).range([H - m.b, m.t]);
  const y2 = d3.scaleLinear().domain([0, d3.max(dens, d => d.p) * 1.15]).range([H - m.b, m.t + 40]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(7).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Number of transients inside the shadow");
  s.append("g").selectAll("rect").data(dens.filter(d => d.x0 < 1300)).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1))
    .attr("y", d => y2(d.p)).attr("height", d => y2(0) - y2(d.p)).attr("rx", 1.5).attr("fill", V("--ours")).attr("fill-opacity", 0.75);
  s.append("path").datum(norm).attr("d", d3.area().x(d => x(d[0])).y0(y(0)).y1(d => y(d[1])).curve(d3.curveBasis))
    .attr("fill", V("--theirs")).attr("fill-opacity", 0.85);
  s.append("line").attr("x1", x(P.observed)).attr("x2", x(P.observed)).attr("y1", m.t).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 2);
  s.append("text").attr("class", "ann").attr("x", x(P.observed) - 8).attr("y", m.t + 10).attr("text-anchor", "end").text(`Observed: ${P.observed}`);
}

function drawPhiSet(D, id, key, vasco) {   // one likelihood curve from shadow_phi_combined.py, styled as drawPhi; vasco: mark 0.32
  const s = svg(id), W = 760, H = 320, m = { l: 62, r: 20, t: 40, b: 46 }, R = D.south.phi[key];
  const x = d3.scaleLinear().domain([0, 1]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 1.08]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(10).tickFormat(d3.format(".1f")).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Fraction of transients that vanish inside the shadow");
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0, 0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".2f")).tickSize(3)).call(axisStyle);
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Likelihood, relative to the best value");
  const c = R.likelihood.phi.map((v, i) => [v, R.likelihood.rel[i]]), up = R.fc90;
  const line = d3.line().x(d => x(d[0])).y(d => y(d[1])).curve(d3.curveMonotoneX);
  const area = d3.area().x(d => x(d[0])).y0(y(0)).y1(d => y(d[1])).curve(d3.curveMonotoneX);
  const k = c.findIndex(d => d[0] > up), yUp = c[k - 1][1] + (c[k][1] - c[k - 1][1]) * (up - c[k - 1][0]) / (c[k][0] - c[k - 1][0]);
  s.append("path").attr("d", area(c.slice(0, k).concat([[up, yUp]]))).attr("fill", V("--ours-soft"));
  s.append("path").attr("d", line(c)).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.6);
  const halo = sel => sel.attr("class", "ann").style("paint-order", "stroke").style("stroke", V("--surface")).style("stroke-width", "4px");
  s.append("line").attr("x1", x(up)).attr("x2", x(up)).attr("y1", y(0)).attr("y2", y(0.55)).attr("stroke", V("--ours")).attr("stroke-width", 1.5);
  s.append("text").attr("x", x(up) + 8).attr("y", y(0.55) + 4).call(halo).style("fill", V("--ours")).text(`90% upper limit: ${up.toFixed(2)}`);
  s.append("circle").attr("cx", x(R.best)).attr("cy", y(1)).attr("r", 6).attr("fill", V("--ours")).attr("stroke", V("--surface")).attr("stroke-width", 2);
  s.append("text").attr("x", x(R.best) + 12).attr("y", y(1) - 8).call(halo).style("fill", V("--ours")).text(`Best estimate: ${R.best.toFixed(2)}`);
  if (!vasco) return;
  const v = 0.32;
  s.append("line").attr("x1", x(v)).attr("x2", x(v)).attr("y1", y(0)).attr("y2", y(0.95)).attr("stroke", V("--theirs")).attr("stroke-width", 1.4);
  s.append("path").attr("d", d3.symbol(d3.symbolDiamond, 80)()).attr("transform", `translate(${x(v)},${y(0.95)})`).attr("fill", V("--theirs"));
  s.append("text").attr("x", x(v) + 9).attr("y", y(0.95) + 4).call(halo).style("fill", V("--theirs")).text("Implied by VASCO’s headline deficit, about 0.32");
}

function drawPhi(D) {
  const s = svg("c-phi"), W = 760, H = 320, m = { l: 62, r: 20, t: 40, b: 46 };
  const x = d3.scaleLinear().domain([0, 1]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 1.08]).range([H - m.b, m.t]);
  const curve = g => { const P = D.phi[g].posterior, mx = d3.max(P.density); return P.phi.map((v, i) => [v, P.density[i] / mx]); };
  const geo = curve("geo"), up = D.phi.geo.fc["0.9"][1];
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(10).tickFormat(d3.format(".1f")).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Fraction of transients that vanish inside the shadow");
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0, 0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".2f")).tickSize(3)).call(axisStyle);
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Likelihood, relative to the best value");
  const area = d3.area().x(d => x(d[0])).y0(y(0)).y1(d => y(d[1])).curve(d3.curveMonotoneX);
  const line = d3.line().x(d => x(d[0])).y(d => y(d[1])).curve(d3.curveMonotoneX);
  const k = geo.findIndex(d => d[0] > up), yUp = geo[k - 1][1] + (geo[k][1] - geo[k - 1][1]) * (up - geo[k - 1][0]) / (geo[k][0] - geo[k - 1][0]);
  s.append("path").attr("d", area(geo.slice(0, k).concat([[up, yUp]]))).attr("fill", V("--ours-soft"));
  s.append("path").attr("d", line(geo)).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.6);
  s.append("line").attr("x1", x(up)).attr("x2", x(up)).attr("y1", y(0)).attr("y2", y(0.55)).attr("stroke", V("--ours")).attr("stroke-width", 1.5);
  s.append("text").attr("class", "ann").attr("x", x(up) + 8).attr("y", y(0.55) + 4).style("fill", V("--ours")).text(`Upper limit at 90% confidence: ${up.toFixed(2)}`);
  const mo = D.phi.geo.mode, my = geo.reduce((b, d) => (Math.abs(d[0] - mo) < Math.abs(b[0] - mo) ? d : b))[1];
  s.append("circle").attr("cx", x(mo)).attr("cy", y(my)).attr("r", 6).attr("fill", V("--ours")).attr("stroke", V("--surface")).attr("stroke-width", 2);
  s.append("text").attr("class", "ann").attr("x", x(mo) + 12).attr("y", y(my) - 8).style("fill", V("--ours")).text(`Best estimate: ${mo.toFixed(2)}`);
  const v = 0.32;
  s.append("line").attr("x1", x(v)).attr("x2", x(v)).attr("y1", y(0)).attr("y2", y(0.95)).attr("stroke", V("--theirs")).attr("stroke-width", 1.4);
  s.append("path").attr("d", d3.symbol(d3.symbolDiamond, 80)()).attr("transform", `translate(${x(v)},${y(0.95)})`).attr("fill", V("--theirs"));
  s.append("text").attr("class", "ann").attr("x", x(v) + 9).attr("y", y(0.95) + 4).style("fill", V("--theirs")).text("Implied by VASCO’s headline deficit, about 0.32");
}



/* ---------- nuclear 3: months ---------- */
function drawMonths(D, id, keys) {
  const all = { n: { label: "Transients per night", c: V("--ours"), max: 450, f: v => Math.round(v) },
                b: { label: "Degrees from the Milky Way", c: V("--teal"), max: 60, f: v => Math.round(v) + "°" },
                win: { label: "Nights in VASCO’s window", c: V("--emulsion"), max: 7, f: v => v || "" } };
  const rows = keys.map(k => ({ key: k, ...all[k] })), W = 760, m = { l: 190, r: 16 }, H = 14 + rows.length * 112 + 8;
  const s = svg(id); s.attr("viewBox", `0 0 ${W} ${H}`);
  const x = d3.scaleBand().domain(D.months.map(d => d.m)).range([m.l, W - m.r]).padding(0.28);
  rows.forEach((r, i) => {
    const top = 14 + i * 112, bot = top + 86;
    const y = d3.scaleLinear().domain([0, r.max]).range([bot, top]);
    s.append("text").attr("x", 0).attr("y", bot - 4).attr("class", "ann").text(r.label);
    s.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", bot).attr("y2", bot).attr("stroke", V("--rule"));
    s.append("g").selectAll("rect").data(D.months).join("rect").attr("class", "hov").attr("x", d => x(d.m)).attr("width", x.bandwidth()).attr("y", d => y(d[r.key])).attr("height", d => bot - y(d[r.key])).attr("rx", 2).attr("fill", r.c)
      .call(sel => hov(sel, d => `<b>${d.m}</b><br>${Math.round(d.n)} transients per night on average<br>Plates typically ${Math.round(d.b)}° from the plane of the Milky Way<br>${d.win} observing night${d.win === 1 ? "" : "s"} within a day of a US test`));
    s.append("g").selectAll("text.v").data(D.months).join("text").attr("x", d => x(d.m) + x.bandwidth() / 2).attr("y", d => y(d[r.key]) - 4).attr("text-anchor", "middle").style("font-size", "10px").text(d => r.f(d[r.key]));
  });
  s.append("g").selectAll("text.m").data(D.months).join("text").attr("x", d => x(d.m) + x.bandwidth() / 2).attr("y", 14 + (rows.length - 1) * 112 + 86 + 18).attr("text-anchor", "middle").text(d => d.m);
}


/* ---------- alignments 1: fan ---------- */
function drawFan(D) {
  const s = svg("c-fan"), W = 760, H = 300, R = 160, cy = 240;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const half = arr => arr.concat(arr.slice(0, 90).reverse());       // 0..180 degrees from east
  const fans = [{ cx: 190, v: half(D.pearson), c: V("--theirs"), title: "VASCO’s selection" },
                { cx: 570, v: half(D.trackfit), c: V("--ours"), title: "Straight-line fit" }];
  fans.forEach(({ cx, v, c, title }) => {
    const pt = (deg, r) => [cx + Math.cos(deg * Math.PI / 180) * r * R, cy - Math.sin(deg * Math.PI / 180) * r * R];
    const wedge = (a0, a1) => d3.arc()({ innerRadius: 0, outerRadius: R * 1.08, startAngle: Math.PI / 2 - a1 * Math.PI / 180, endAngle: Math.PI / 2 - a0 * Math.PI / 180 });
    for (const [a0, a1] of [[0, 8], [172, 180]]) s.append("path").attr("d", wedge(a0, a1)).attr("transform", `translate(${cx},${cy})`).attr("fill", V("--teal")).attr("fill-opacity", 0.35);
    for (const rr of [0.5, 1]) s.append("path").attr("d", d3.arc()({ innerRadius: rr * R, outerRadius: rr * R, startAngle: -Math.PI / 2, endAngle: Math.PI / 2 })).attr("transform", `translate(${cx},${cy})`).attr("fill", "none").attr("stroke", V("--rule"));
    s.append("path").attr("d", "M" + cx + "," + cy + "L" + v.map((x, a) => pt(a, x).join(",")).join("L") + "Z").attr("fill", c).attr("fill-opacity", 0.3).attr("stroke", c).attr("stroke-width", 2);
    s.append("text").attr("class", "ann").attr("x", cx).attr("y", 18).attr("text-anchor", "middle").style("font-weight", 600).text(title);
    s.append("text").attr("x", cx).attr("y", cy - R - 8).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text("North");
    s.append("text").attr("x", cx + R * 1.08 + 4).attr("y", cy + 4).attr("fill", V("--ink-3")).text("E");
    s.append("text").attr("x", cx - R * 1.08 - 4).attr("y", cy + 4).attr("text-anchor", "end").attr("fill", V("--ink-3")).text("W");
    const ew = Math.round((v[0] + v[180]) / 2 * 100);
    s.append("text").attr("class", "ann").attr("x", cx).attr("y", cy + 26).attr("text-anchor", "middle").style("fill", V("--teal"))
      .text(`East–west lines kept: ${ew}%`);
    s.append("g").selectAll("path.h").data(d3.range(0, 181, 3)).join("path")
      .attr("d", a => d3.arc()({ innerRadius: 0, outerRadius: R, startAngle: Math.PI / 2 - (a + 1.5) * Math.PI / 180, endAngle: Math.PI / 2 - (a - 1.5) * Math.PI / 180 })).attr("transform", `translate(${cx},${cy})`).attr("fill", "transparent")
      .call(sel => hov(sel, a => `${title}<br>A line ${a <= 90 ? a : 180 - a}° from east–west: ${Math.round(v[a] * 100)}% of perfect lines kept`));
  });
}


/* ---------- alignments: VASCO's search and score, real against shuffled copies ---------- */
function drawTracks(D) {
  const T = D.trackdist, s = svg("c-counts"), W = 760, m = { l: 52, r: 10, t: 44, b: 52 }, gx = 12;
  const rows = [{ key: "dist", h: 160, lab: "Number of track candidates at each fit quality" }, { key: "all", h: 110, lab: "Difference between real catalogue and random chance (standard deviations)" }];
  const gy = 34, H = m.t + d3.sum(rows, r => r.h) + gy * (rows.length - 1) + m.b;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const pw = (W - m.l - m.r - 3 * gx) / 4, step = 0.2, lo = T.lo;
  const fmtX = v => ({ "-4": "10⁻⁴", "-2": "0.01", "0": "1", "2": "100" })[String(v)];
  let top = m.t;
  rows.forEach((row, ri) => {
    const bot = top + row.h;
    s.append("text").attr("class", "ann").attr("x", m.l).attr("y", top - 8).style("font-weight", 600).text(row.lab);
    [0, 1, 2, 3].forEach(ci => {
      const x0 = m.l + ci * (pw + gx), x = d3.scaleLinear().domain([-4, 3]).range([x0, x0 + pw]);
      const d = T[row.key === "inner" ? "inner" : "all"][ci];
      if (ri === 0) s.append("text").attr("x", x0 + pw / 2).attr("y", 14).attr("text-anchor", "middle").style("font-weight", 600).text(`${d.r}-point candidates`);
      s.append("rect").attr("x", x0).attr("y", top).attr("width", pw).attr("height", row.h).attr("fill", "none").attr("stroke", V("--rule"));
      const bins = lo.map((l, k) => ({ l, h: l + step, n: d.real[k], b: d.bmean[k], sd: d.bsd[k], sg: d.sig[k] }));
      if (row.key === "dist") {
        const y = d3.scaleLog().domain([0.5, d3.max(bins, q => Math.max(q.n, q.b)) * 3]).range([bot, top]).clamp(true);
        const stepPath = acc => d3.line().x(q => x(q[0])).y(q => y(Math.max(0.5, q[1]))).curve(d3.curveStepAfter)(bins.map(q => [q.l, acc(q)]).concat([[bins[bins.length - 1].h, acc(bins[bins.length - 1])]]));
        s.append("path").attr("d", stepPath(q => q.b) + `L${x(3)},${bot}L${x(-4)},${bot}Z`).attr("fill", V("--ink-3")).attr("fill-opacity", 0.35);
        const pts = bins.filter(q => q.n > 0);
        s.append("g").selectAll("line").data(pts).join("line").attr("x1", q => x(q.l + step / 2)).attr("x2", q => x(q.l + step / 2))
          .attr("y1", q => y(Math.max(0.5, q.n - Math.sqrt(q.n)))).attr("y2", q => y(q.n + Math.sqrt(q.n))).attr("stroke", V("--emulsion")).attr("stroke-width", 1.4);
        s.append("g").selectAll("circle").data(pts).join("circle").attr("cx", q => x(q.l + step / 2)).attr("cy", q => y(q.n)).attr("r", 2.6).attr("fill", V("--emulsion"));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([1, 100, 10000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
      } else {
        const y = d3.scaleLinear().domain([-4, 4]).range([bot, top]).clamp(true);
        for (const [k, op] of [[2, 0.18], [1, 0.3]]) s.append("rect").attr("x", x0).attr("width", pw).attr("y", y(k)).attr("height", y(-k) - y(k)).attr("fill", V("--ink-3")).attr("fill-opacity", op);
        s.append("line").attr("x1", x0).attr("x2", x0 + pw).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", V("--ink-3")).attr("stroke-dasharray", "3 3");
        const ok = bins.filter(q => q.sd > 0 && q.b >= 1);   // pulls mean nothing where chance expects under one group
        s.append("g").selectAll("circle").data(ok).join("circle").attr("class", "hov").attr("cx", q => x(q.l + step / 2)).attr("cy", q => y((q.n - q.b) / q.sd)).attr("r", 2.8).attr("fill", V("--emulsion"))
          .call(sel => hov(sel, q => `${d.r}-point candidates, χ²/ndof ${(10 ** q.l).toPrecision(2)} to ${(10 ** q.h).toPrecision(2)}<br>Real catalogue ${d3.format(",")(q.n)}; chance ${d3.format(",.1f")(q.b)} ± ${q.sd.toFixed(1)}<br>${((q.n - q.b) / q.sd).toFixed(1)} standard deviations`));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([-4, -2, 0, 2, 4]).tickFormat(v => (v > 0 ? "+" : "") + v + "σ").tickSize(3)).call(axisStyle);
      }
      if (ri === rows.length - 1) s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).tickValues([-4, -2, 0, 2]).tickFormat(fmtX).tickSize(3).tickSizeOuter(0)).call(axisStyle);
    });
    top = bot + gy;
  });
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 12).attr("text-anchor", "middle").text("Track χ²/ndof, for 1″ hit resolution (lower values are better fits)");
}
function wrap(t, width) {   // split an SVG label over lines to fit a width
  const words = t.text().split(" "), x = t.attr("x"), y = +t.attr("y"); t.text(null);
  let line = [], tsp = t.append("tspan").attr("x", x).attr("y", y), n = 0;
  for (const w of words) {
    line.push(w); tsp.text(line.join(" "));
    if (tsp.node().getComputedTextLength() > width && line.length > 1) { line.pop(); tsp.text(line.join(" ")); line = [w]; tsp = t.append("tspan").attr("x", x).attr("y", y + ++n * 14).text(w); }
  }
}

/* ---------- dates ---------- */
function fillDates(D) {
  const notes = { XE113: "One day early; rest of its run correct", XE389: "Correct; rest of its run one day early" };
  document.querySelector("#t-exc tbody").innerHTML = D.date_summary.exceptions.map(e => `<tr><td class="m">${e.plate}</td><td class="m">E-${e.e}</td><td class="m">${e.header}</td><td class="m">${e.correct}</td><td>${notes[e.plate] || "Out by " + Math.abs(e.days) + " days"}</td></tr>`).join("");
}

/* ---------- catalogue ---------- */
function drawRadial(D) {
  const R = D.radial, s = svg("c-radial"), W = 760, H = 320, m = { l: 44, r: 16, t: 14, b: 42 };
  const x = d3.scaleLinear().domain([0, 3]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 3]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6).tickFormat(d => d + "°").tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickFormat(d => d + "×").tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Distance from the centre of the plate");
  const line = k => d3.line().x((v, i) => x(R.r[i])).y(v => y(v)).curve(d3.curveMonotoneX)(R[k]);
  s.append("path").attr("d", d3.area().x((v, i) => x(R.r[i])).y0((v, i) => y(R.lo[i])).y1((v, i) => y(R.hi[i])).curve(d3.curveMonotoneX)(R.r)).attr("fill", V("--ours-soft"));
  s.append("path").attr("d", line("control")).attr("fill", "none").attr("stroke", V("--ink-3")).attr("stroke-width", 2).attr("stroke-dasharray", "4 3");
  s.append("path").attr("d", line("transients")).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.4);
  s.append("g").selectAll("circle").data(R.transients).join("circle").attr("cx", (v, i) => x(R.r[i])).attr("cy", v => y(v)).attr("r", 3.5).attr("fill", V("--ours"))
    .call(sel => hov(sel, (v, i) => `${v.toFixed(1)} times the density near the centre`));
}
function drawCrowd(D) {
  const s = svg("c-crowd"), W = 760, H = 320, m = { l: 52, r: 16, t: 14, b: 42 };
  const pts = D.crowding;
  const half = 6;   // quartiles over a sliding 12-degree window, stepped every 2 degrees
  const bins = d3.range(2, 79, 2).map(c => {
    const v = pts.filter(p => Math.abs(p.b - c) <= half).map(p => p.N).sort(d3.ascending);
    return { b: c, n: v.length, q25: d3.quantile(v, 0.25), q50: d3.quantile(v, 0.5), q75: d3.quantile(v, 0.75) };
  }).filter(d => d.n >= 20);
  const w = [1, 2, 3, 2, 1];   // light triangular smoothing of the three curves across neighbouring points
  for (const k of ["q25", "q50", "q75"]) {
    const raw = bins.map(d => d[k]);
    bins.forEach((d, i) => { let t = 0, ws = 0; w.forEach((wj, j) => { const v = raw[i + j - 2]; if (v !== undefined) { t += wj * v; ws += wj; } }); d[k] = t / ws; });
  }
  const x = d3.scaleLinear().domain([0, 80]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 400]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(8).tickFormat(d => d + "°").tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Degrees from the plane of the Milky Way");
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients per plate");
  s.append("path").datum(bins).attr("d", d3.area().x(d => x(d.b)).y0(d => y(d.q25)).y1(d => y(d.q75)).curve(d3.curveMonotoneX)).attr("fill", V("--ours-soft"));
  s.append("path").datum(bins).attr("d", d3.line().x(d => x(d.b)).y(d => y(d.q50)).curve(d3.curveMonotoneX)).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.4);
  s.append("g").selectAll("circle").data(bins).join("circle").attr("cx", d => x(d.b)).attr("cy", d => y(d.q50)).attr("r", 7).attr("fill", "transparent")
    .call(sel => hov(sel, d => `${Math.max(0, d.b - half)}° to ${d.b + half}° from the Milky Way: ${d.n} plates<br>Median ${Math.round(d.q50)} transients; middle half ${Math.round(d.q25)} to ${Math.round(d.q75)}`));
}



/* ---------- nuclear: VASCO's random-night test, before and after allowing for where the telescope pointed ---------- */
function drawPerm(D, id, key, lab) {
  const q = D.perm[key], s = svg(id), W = 760, H = 250, m = { l: 52, r: 16, t: 34, b: 44 }, bw = 0.05;
  const x = d3.scaleLinear().domain([0, 3]).range([m.l, W - m.r]);
  const top = m.t, bot = H - m.b, tot = d3.sum(q.hist);
  const y = d3.scaleLinear().domain([0, 0.12]).range([bot, top + 8]);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0, 0.05, 0.1]).tickFormat(d3.format(".0%")).tickSize(3)).call(axisStyle);
  s.append("g").selectAll("rect").data(q.hist).join("rect").attr("x", (c, j) => x(j * bw) + 0.75).attr("width", x(bw) - x(0) - 1.5)
    .attr("y", c => y(c / tot)).attr("height", c => bot - y(c / tot)).attr("fill", V("--ink-3")).attr("fill-opacity", (c, j) => j * bw >= q.real - 1e-9 ? 0.9 : 0.45);
  s.append("line").attr("x1", x(q.real)).attr("x2", x(q.real)).attr("y1", top).attr("y2", bot).attr("stroke", V("--emulsion")).attr("stroke-width", 2.5);
  s.append("text").attr("class", "ann").attr("x", m.l).attr("y", 16).style("font-weight", 600).text(lab);
  const sd = `${q.sigma.toFixed(1)} standard deviation${q.sigma.toFixed(1) === "1.0" ? "" : "s"}`;
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 14).style("fill", V("--emulsion")).text(`Nights near tests: ${q.real.toFixed(2)}×`);
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 30).style("fill", V("--emulsion")).text(`Random sets reach it 1 time in ${Math.round(1 / q.p)} (${sd})`);
  s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).ticks(6).tickFormat(d => d + "×").tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Transients on the chosen nights relative to all other nights");
  s.append("text").attr("transform", `translate(12,${(top + bot) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Share of random sets");
}

let DATA = null;
function render() {
  if (!DATA) return;
  drawSky(DATA, "c-sky0", false); drawSpread(DATA); drawSky(DATA, "c-sky", true); drawChance(DATA); drawPhi(DATA); drawMonths(DATA, "c-months", ["n", "b"]); drawMonths(DATA, "c-window", ["win"]); drawPerm(DATA, "c-perm-any", "raw", "VASCO’s test"); drawPerm(DATA, "c-perm-adj", "adj", "VASCO’s test, with the Milky Way effect removed");
  drawFan(DATA); drawTracks(DATA); fillDates(DATA); drawRadial(DATA); drawCrowd(DATA);
  if (DATA.south) { drawSky(DATA, "c-south-sky", true, true); if (DATA.south.phi) { drawPhiSet(DATA, "c-phi-south", "south", false); drawPhiSet(DATA, "c-phi-comb", "combined", false); } }
}
fetch("data/site.json").then(r => r.json()).then(d => { DATA = d; render(); })
  .catch(() => document.querySelectorAll(".chart svg").forEach(s => s.insertAdjacentHTML("afterend", '<div class="loading">The chart data could not be loaded.</div>')));
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
