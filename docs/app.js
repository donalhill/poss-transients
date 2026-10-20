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
function drawSky(D, id, shadow) {   // shadow: pick out the plates the shadow crossed
  const s = svg(id), W = 760, H = 760, cx = W / 2, cy = H / 2, Rd = W / 2 - 44;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const plates = D.plates.map(p => Object.assign({}, p, { g: ring(p.poly) }));
  const clip = d3.max(plates.flatMap(p => p.poly), v => 90 - v[1]) + 2.5;
  const proj = d3.geoAzimuthalEqualArea().reflectX(true).rotate([0, -90]).clipAngle(clip)
    .scale(Rd / (2 * Math.sin(clip * DEG / 2))).translate([cx, cy]).precision(0.3);
  const path = d3.geoPath(proj);
  const col = d3.scaleSequentialLog(d3.interpolateViridis).domain([5, 1500]).clamp(true);
  const plane = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => galToEq(l, 0)) };

  s.append("path").datum({ type: "Sphere" }).attr("d", path).attr("fill", V("--paper")).attr("stroke", V("--rule"));
  s.append("path").datum(d3.geoGraticule().step([30, 30]).extentMajor([[-180, -80], [180, 80.01]])()).attr("d", path)
    .attr("fill", "none").attr("stroke", V("--rule")).attr("stroke-width", 0.7);
  s.append("g").selectAll("path").data(plates.filter(p => !(shadow && p.circle))).join("path").attr("class", "hov").attr("d", p => path(p.g))
    .attr("fill", p => p.N > 0 ? col(Math.max(p.N, 5)) : V("--surface")).attr("fill-opacity", shadow ? 0.3 : 1).attr("stroke", V("--surface")).attr("stroke-width", 0.5)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients` + (p.vasco ? "" : "<br>Left out of VASCO’s published tests: its archive date is not on VASCO’s list of observing nights")));
  // plates VASCO's published tests leave out, outlined
  if (!shadow) s.append("g").attr("pointer-events", "none").selectAll("path").data(plates.filter(p => !p.vasco)).join("path").attr("d", p => path(p.g))
    .attr("fill", "none").attr("stroke", "#111418").attr("stroke-width", 1.3).attr("stroke-dasharray", "2.5 2");
  // plane of the Milky Way, drawn over the tiles: a soft band with a dashed centre line
  const mw = s.append("g").attr("pointer-events", "none");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-opacity", 0.38).attr("stroke-width", 16).attr("stroke-linecap", "round");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#111418").attr("stroke-opacity", 0.85).attr("stroke-width", 1.5).attr("stroke-dasharray", "5 4");
  const sh = shadow ? plates.filter(p => p.circle) : [];
  // the shadow plates outlined in black on a white halo; the shadow itself as a ring in a hue viridis never uses
  s.append("g").selectAll("path").data(sh).join("path").attr("pointer-events", "none").attr("d", p => path(p.g)).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-width", 5);
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "hov").attr("d", p => path(p.g))
    .attr("fill", p => p.N !== undefined ? col(Math.max(p.N, 5)) : V("--surface")).attr("stroke", "#111418").attr("stroke-width", 1.8)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients on the whole plate<br>${Math.round(p.frac * 100)}% of the whole plate inside the shadow<br>${p.obs} transients in the shadowed part; an average plate would have about ${Math.round(p.exp)} there`));

  s.append("g").attr("pointer-events", "none").selectAll("path").data(sh).join("path").attr("d", p => path(ring(p.circle)))
    .attr("fill", "none").attr("stroke", SHADOW).attr("stroke-width", 2.4);
  const darkHalo = sel => sel.attr("class", "ann").style("font-weight", 600).style("fill", "#111418").style("paint-order", "stroke").style("stroke", "#ffffff").style("stroke-width", "3px");
  for (let h = 0; h < 24; h += 2) {   // hours of right ascension round the rim
    const [x, y] = proj([h * 15, 90 - clip + 0.5]), r = Math.hypot(x - cx, y - cy), R = Rd + 16;
    s.append("text").attr("x", cx + (x - cx) / r * R).attr("y", cy + (y - cy) / r * R + 4).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text(h + "h");
  }
  for (const d of [0, 30, 60]) {
    const [x, y] = proj([172, d]);
    s.append("text").attr("x", x).attr("y", y - 3).attr("text-anchor", "middle").attr("fill", V("--ink-3"))
      .style("paint-order", "stroke").style("stroke", V("--paper")).style("stroke-width", "3px").text((d > 0 ? "+" : "") + d + "°");
  }
  // colour scale: a log bar in the top-left corner, outside the disc
  const gid = `${id}-grad`, bx = 8, by = 26, bw = 200, lx = d3.scaleLog().domain(col.domain()).range([bx, bx + bw]);
  s.append("defs").append("linearGradient").attr("id", gid).selectAll("stop").data(d3.range(0, 1.001, 0.1)).join("stop")
    .attr("offset", t => t).attr("stop-color", t => col(lx.invert(bx + t * bw)));
  s.append("text").attr("x", bx).attr("y", by - 8).attr("fill", V("--ink-2")).text("Transients per plate");
  s.append("rect").attr("x", bx).attr("y", by).attr("width", bw).attr("height", 10).attr("fill", `url(#${gid})`).attr("fill-opacity", shadow ? 0.3 : 1);
  s.append("text").attr("x", bx + bw + 6).attr("y", by + 9).attr("fill", V("--ink-3")).style("font-size", "11px").text("1,500+");
  [5, 10, 100, 1000].forEach(v => {
    s.append("line").attr("x1", lx(v)).attr("x2", lx(v)).attr("y1", by + 10).attr("y2", by + 14).attr("stroke", V("--ink-3"));
    s.append("text").attr("x", lx(v)).attr("y", by + 25).attr("text-anchor", v === 5 ? "start" : "middle").attr("fill", V("--ink-3")).style("font-size", "11px").text(fmt(v));
  });
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
  s.append("text").attr("class", "ann").attr("x", x(hi) + 8).attr("y", m.t + 16).text(`${inBand} of ${N.length} plates fall in this range`);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 200, 400, 600, 800]).tickFormat(d3.format(",")).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", x(1050)).attr("y", H - m.b + 15).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text("Over 1,000");
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text("Transients per plate");
  s.append("line").attr("class", "avg").attr("x1", x(mu)).attr("x2", x(mu)).attr("y1", m.t - 14).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 1.5).attr("stroke-dasharray", "4 3");
  s.append("text").attr("class", "ann").attr("x", x(hi) + 8).attr("y", m.t - 2).text(`Average (dashed line): ${Math.round(mu)}`);
  s.append("text").attr("transform", `translate(12,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Number of plates");
  s.append("g").selectAll("rect.b").data(bins.concat([over])).join("rect").attr("class", "hov").attr("x", b => x(b.lo) + 0.75).attr("width", b => x(b.hi) - x(b.lo) - 1.5)
    .attr("y", b => y(b.n)).attr("height", b => y(0) - y(b.n)).attr("fill", V("--ours"))
    .call(sel => hov(sel, b => b.more ? `${b.n} plates hold more than 1,000 transients, up to ${fmt(d3.max(N))}` : `${b.n} plates hold ${fmt(b.lo)} to ${fmt(b.hi - 1)} transients`));
  s.select(".avg").raise();
}
function erf(x) {   // Abramowitz and Stegun 7.1.26
  const t = 1 / (1 + 0.3275911 * Math.abs(x)), y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}

/* ---------- 2. two ideas of chance ---------- */
function drawChance(D) {
  const P = D.paste, s = svg("c-chance"), W = 760, H = 330, m = { l: 44, r: 16, t: 14, b: 40 };
  const x = d3.scaleLinear().domain([0, 1600]).range([m.l, W - m.r]);
  const binw = P.bins[1] - P.bins[0];
  const dens = P.counts.map((c, i) => ({ x0: P.bins[i], x1: P.bins[i + 1], p: c / P.n / binw }));
  const mu = P.area_expected, sd = P.binom_sd;
  const norm = d3.range(mu - 5 * sd, mu + 5 * sd, 0.5).map(v => [v, Math.exp(-0.5 * ((v - mu) / sd) ** 2) / (sd * Math.sqrt(2 * Math.PI))]);
  const y = d3.scaleLinear().domain([0, d3.max(norm, d => d[1]) * 1.08]).range([H - m.b, m.t]);
  const y2 = d3.scaleLinear().domain([0, d3.max(dens, d => d.p) * 1.15]).range([H - m.b, m.t + 40]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(7).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Number of transients inside the shadow");
  s.append("g").selectAll("rect").data(dens.filter(d => d.x0 < 1600)).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1))
    .attr("y", d => y2(d.p)).attr("height", d => y2(0) - y2(d.p)).attr("rx", 1.5).attr("fill", V("--ours")).attr("fill-opacity", 0.75);
  s.append("path").datum(norm).attr("d", d3.area().x(d => x(d[0])).y0(y(0)).y1(d => y(d[1])).curve(d3.curveBasis))
    .attr("fill", V("--theirs")).attr("fill-opacity", 0.85);
  s.append("line").attr("x1", x(P.observed)).attr("x2", x(P.observed)).attr("y1", m.t).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 2);
  s.append("text").attr("class", "ann").attr("x", x(P.observed) - 8).attr("y", m.t + 10).attr("text-anchor", "end").text(`Observed: ${P.observed}`);
  s.append("text").attr("class", "ann").attr("x", W - m.r).attr("y", m.t + 10).attr("text-anchor", "end").style("fill", V("--ink-3")).style("font-size", "11px")
    .text("Each shape scaled to fit; the orange curve is about 15 times narrower and taller");
}



/* ---------- the plates partly crossed by Earth's shadow, up close: transients as dots over their density (north
   up, east left), coloured by the side of the shadow's edge they lie on, with each plate's numbers ---------- */
function drawPlates(D, id) {
  const Z = D.shadow, s = svg(id);
  if (!Z || !Z.plates) return;
  const cols = Z.plates.length === 9 ? 3 : 4, W = 760, colW = W / cols, S = cols === 3 ? 196 : 156, rowH = S + 64, top = 4, half = Z.half;
  const H = top + Math.ceil(Z.plates.length / cols) * rowH;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const px = v => (half - v) / (2 * half) * S, py = v => (half - v) / (2 * half) * S;
  Z.plates.forEach((p, k) => {
    const x0 = (k % cols) * colW + (colW - S) / 2, y0 = top + Math.floor(k / cols) * rowH;
    const g = s.append("g").attr("transform", `translate(${x0},${y0})`);
    g.append("clipPath").attr("id", `plate-clip-${k}`).append("rect").attr("width", S).attr("height", S);
    const inner = g.append("g").attr("clip-path", `url(#plate-clip-${k})`);
    inner.append("rect").attr("width", S).attr("height", S).attr("fill", V("--paper"));
    const pts = p.xi.map((v, i) => ({ x: px(v), y: py(p.eta[i]), dark: p.dark[i] === "1" }));
    // the shadow's edge: the 8.5 deg circle about the anti-sun point, projected on the plate
    inner.append("path").attr("d", d3.line().curve(d3.curveLinearClosed)(p.edge.map(e => [px(e[0]), py(e[1])])))
      .attr("fill", SHADOW).attr("fill-opacity", 0.12).attr("stroke", SHADOW).attr("stroke-width", 2);
    inner.append("g").selectAll("circle").data(pts).join("circle").attr("cx", d => d.x).attr("cy", d => d.y).attr("r", 1.25)
      .attr("fill", d => d.dark ? SHADOW : V("--ours")).attr("fill-opacity", 0.9);
    g.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", V("--ink-3"));
    g.append("text").attr("x", S / 2).attr("y", S + 18).attr("text-anchor", "middle").style("font-weight", 600).style("font-size", "13px")
      .text(`${p.plate} · ${Math.round(p.shadow * 100)}% in shadow`);
    g.append("text").attr("class", "ann").attr("x", S / 2).attr("y", S + 35).attr("text-anchor", "middle").style("fill", V("--ink-2")).style("font-size", "12px")
      .text(`${p.O} in shadow, ${p.E.toFixed(1)} expected`);
    g.append("text").attr("class", "ann").attr("x", S / 2).attr("y", S + 51).attr("text-anchor", "middle").style("fill", V("--ink-2")).style("font-size", "12px")
      .text(`ratio ${p.r.toFixed(2)}`);
    g.append("rect").attr("class", "hov").attr("width", S).attr("height", S).attr("fill", "transparent")
      .call(sel => hov(sel, () => `<b>Plate ${p.plate}</b><br>${fmt(p.n)} transients, ${Math.round(p.shadow * 100)}% of the plate in shadow<br>` +
        `${fmt(p.sunlit)} in sunlight, ${p.O} in shadow<br>expected in shadow with no effect: ${p.E.toFixed(1)}<br>` +
        `ratio ${p.r.toFixed(2)} (${p.lo.toFixed(2)} to ${p.hi == null ? "–" : p.hi.toFixed(2)})`));
  });
}

/* ---------- the ratio plate by plate (shadow_fit.py): each plate's ratio with its counting range, and the best fit across
   all of them with its counting range (one standard deviation) ---------- */
function drawForest(D, id) {
  const Z = D.shadow, s = svg(id);
  if (!Z || !Z.fit) return;
  const rows = Z.plates, F = Z.fit;
  const W = 760, rowH = 26, m = { l: 150, r: 24, t: 34, b: 52 }, XMAX = 12;
  const H = m.t + (rows.length + 1.6) * rowH + m.b;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  // a log axis, with zero (plates with none in shadow) set apart on its left beyond a break
  const xl = d3.scaleLog().domain([0.2, XMAX]).range([m.l + 34, W - m.r]), x = v => (v <= 0 ? m.l : xl(Math.max(v, 0.2)));
  const y = i => m.t + (i + 0.5) * rowH, yC = m.t + (rows.length + 1.0) * rowH;
  s.append("g").attr("transform", `translate(0,${H - m.b})`)
    .call(d3.axisBottom(xl).tickValues([0.25, 0.5, 1, 2, 4, 8]).tickFormat(d3.format("~g")).tickSizeOuter(0)).call(axisStyle);
  const ax0 = s.append("g").attr("class", "ax").attr("transform", `translate(${m.l},${H - m.b})`);
  ax0.append("line").attr("y2", 6).attr("stroke", "currentColor");
  ax0.append("text").attr("y", 9).attr("dy", "0.71em").attr("text-anchor", "middle").attr("fill", "currentColor").text("0");
  for (const dx of [14, 19]) s.append("line").attr("x1", m.l + dx - 2).attr("x2", m.l + dx + 2).attr("y1", H - m.b + 5).attr("y2", H - m.b - 5)
    .attr("stroke", V("--ink-3")).attr("stroke-width", 1.2);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 12).attr("text-anchor", "middle")
    .text("r = density of transients in shadow ÷ density in sunlight");
  s.append("line").attr("x1", x(1)).attr("x2", x(1)).attr("y1", m.t - 8).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 1.2);
  s.append("text").attr("class", "ann").attr("x", x(1)).attr("y", m.t - 14).attr("text-anchor", "middle").style("fill", V("--ink-2")).text("No effect");
  const bar = (lo, hi, yy, col, w, cap) => {
    const h = hi == null ? XMAX : Math.min(hi, XMAX);
    s.append("line").attr("x1", x(lo)).attr("x2", x(h)).attr("y1", yy).attr("y2", yy).attr("stroke", col).attr("stroke-width", w);
    if (cap) s.append("line").attr("x1", x(lo)).attr("x2", x(lo)).attr("y1", yy - cap).attr("y2", yy + cap).attr("stroke", col).attr("stroke-width", w);
    if (hi == null || hi > XMAX) s.append("path").attr("d", d3.symbol(d3.symbolTriangle, 34)()).attr("transform", `translate(${x(XMAX) - 3},${yy}) rotate(90)`).attr("fill", col);
    else if (cap) s.append("line").attr("x1", x(h)).attr("x2", x(h)).attr("y1", yy - cap).attr("y2", yy + cap).attr("stroke", col).attr("stroke-width", w);
  };
  rows.forEach((r, i) => {
    s.append("text").attr("x", m.l - 10).attr("y", y(i) + 4).attr("text-anchor", "end").style("fill", V("--ink-2"))
      .text(`${r.plate} · ${Math.round(r.shadow * 100)}% in shadow`);
    bar(r.lo, r.hi, y(i), V("--ink-3"), 1.6, 5);
    s.append("circle").attr("class", "hov").attr("cx", x(r.r)).attr("cy", y(i)).attr("r", 4.5).attr("fill", V("--ink-2"))
      .call(sel => hov(sel, () => `<b>Plate ${r.plate}</b><br>${r.O} of ${fmt(r.n)} transients in shadow, ${r.E.toFixed(1)} expected with no effect<br>ratio ${r.r.toFixed(2)} (${r.lo.toFixed(2)} to ${r.hi == null ? "–" : r.hi.toFixed(2)}, counting only)`));
  });
  s.append("line").attr("x1", m.l - 140).attr("x2", W - m.r).attr("y1", yC - rowH * 0.8).attr("y2", yC - rowH * 0.8).attr("stroke", V("--rule"));
  s.append("text").attr("x", m.l - 10).attr("y", yC + 4).attr("text-anchor", "end").style("font-weight", 600).style("fill", V("--ours")).text("Best fit");
  bar(F.count_lo, F.count_hi, yC, V("--ours"), 2.2, 6);
  s.append("circle").attr("cx", x(F.r)).attr("cy", yC).attr("r", 5.5).attr("fill", V("--ours"));
  s.append("text").attr("class", "ann").attr("x", x(F.count_hi) + 10).attr("y", yC + 4).style("fill", V("--ours"))
    .text(`${F.r.toFixed(2)} (${F.count_lo.toFixed(2)} to ${F.count_hi.toFixed(2)})`);
}

/* ---------- the fitted ratio on toys with no shadow effect: the real shadow outlines laid on other plates,
   against the ratio fitted on the real plates ---------- */
function drawRatioNull(D, id) {
  const Z = D.shadow, s = svg(id);
  if (!Z || !Z.toys) return;
  const R = Z.toys, W = 760, H = 320, m = { l: 20, r: 20, t: 40, b: 46 }, bw = 3 / R.hist_null.length, n = d3.sum(R.hist_null);
  const bins = R.hist_null.map((c, i) => ({ x0: i * bw, x1: (i + 1) * bw, p: c / n }));
  const x = d3.scaleLinear().domain([0.2, 2.6]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, d3.max(bins, d => d.p) * 1.1]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(12).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle")
    .text(`Best fit ratio (r) across ${Z.plates.length} plates`);
  s.append("g").selectAll("rect").data(bins.filter(d => d.x0 >= 0.2 && d.x1 <= 2.6)).join("rect")
    .attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1)).attr("y", d => y(d.p)).attr("height", d => y(0) - y(d.p)).attr("rx", 1.5)
    .attr("fill", V("--ink-3")).attr("fill-opacity", 0.55);
  const halo = sel => sel.attr("class", "ann").style("paint-order", "stroke").style("stroke", V("--surface")).style("stroke-width", "4px");
  s.append("line").attr("x1", x(Z.fit.r)).attr("x2", x(Z.fit.r)).attr("y1", m.t - 6).attr("y2", H - m.b).attr("stroke", V("--ours")).attr("stroke-width", 2.4);
  s.append("text").attr("x", x(Z.fit.r) + 8).attr("y", m.t + 4).call(halo).style("fill", V("--ours"))
    .text(`Fit result on real plates (${Z.fit.r.toFixed(2)})`);
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


/* ---------- the track score against VASCO's score: candidate by candidate where both score (left), and the share of
   injected tracks each recovers by direction at the same false-alarm rate (right) ---------- */
function drawCompare(D, id) {
  const C = D.compare, Rc = D.recovery, s = svg(id);
  if (!C || !Rc) return;
  const W = 760, H = 410, top = 44, S = 290, x0 = 58, b = top + S;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const fmt = { "-6": "10⁻⁶", "-4": "10⁻⁴", "-2": "0.01", "0": "1" };
  // left: heatmap
  const x = d3.scaleLinear().domain([-7, 0]).range([x0, x0 + S]), y = d3.scaleLinear().domain([-7, 0]).range([b, top]);
  const cells = []; C.hist.forEach((row, i) => row.forEach((n, j) => { if (n > 0 && C.lo[i] >= -7 && C.lo[j] >= -7) cells.push({ i: C.lo[i], j: C.lo[j], n }); }));
  const col = d3.scaleSequentialLog(d3.interpolateBlues).domain([1, d3.max(cells, c => c.n)]);
  s.append("rect").attr("x", x0).attr("y", top).attr("width", S).attr("height", S).attr("fill", V("--surface")).attr("stroke", V("--rule"));
  s.append("g").selectAll("rect.c").data(cells).join("rect").attr("x", c => x(c.i)).attr("y", c => y(c.j + 0.2)).attr("width", x(0.2) - x(0)).attr("height", y(0) - y(0.2))
    .attr("fill", c => col(Math.max(c.n, 1.5)));
  s.append("line").attr("x1", x(-7)).attr("y1", y(-7)).attr("x2", x(0)).attr("y2", y(0)).attr("stroke", V("--ink-3")).attr("stroke-dasharray", "3 3");
  s.append("g").attr("transform", `translate(0,${b})`).call(d3.axisBottom(x).tickValues([-6, -4, -2, 0]).tickFormat(v => fmt[v]).tickSize(3)).call(axisStyle);
  s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([-6, -4, -2, 0]).tickFormat(v => fmt[v]).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", x0 + S / 2).attr("y", b + 34).attr("text-anchor", "middle").style("font-size", "12px").text("VASCO’s score");
  s.append("text").attr("transform", `translate(${x0 - 44},${top + S / 2}) rotate(-90)`).attr("text-anchor", "middle").style("font-size", "12px").text("Track score");
  s.append("text").attr("x", x0).attr("y", top - 22).style("font-weight", 600).style("font-size", "12.5px").text("Track candidates both scores rate");
  s.append("text").attr("x", x0 + 10).attr("y", top + 18).style("font-size", "11.5px").style("fill", V("--ink-2"))
    .text(`Rank correlation ${C.rho.toFixed(2)}`);
  // right: recovery by direction
  const x1 = 440, w1 = W - x1 - 12, xr = d3.scaleLinear().domain([0, 90]).range([x1, x1 + w1]), yr = d3.scaleLinear().domain([0, 100]).range([b, top]);
  for (const [a0, a1] of [[0, 5], [85, 90]]) s.append("rect").attr("x", xr(a0)).attr("width", xr(a1) - xr(a0)).attr("y", top).attr("height", S).attr("fill", V("--teal")).attr("fill-opacity", 0.12);
  s.append("rect").attr("x", x1).attr("y", top).attr("width", w1).attr("height", S).attr("fill", "none").attr("stroke", V("--rule"));
  s.append("g").attr("transform", `translate(0,${b})`).call(d3.axisBottom(xr).tickValues([0, 15, 30, 45, 60, 75, 90]).tickFormat(v => `${v}°`).tickSize(3)).call(axisStyle);
  s.append("g").attr("transform", `translate(${x1},0)`).call(d3.axisLeft(yr).tickValues([0, 25, 50, 75, 100]).tickFormat(v => `${v}%`).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", x1 + w1 / 2).attr("y", b + 34).attr("text-anchor", "middle").style("font-size", "12px").text("Direction of the artificial track, from east–west");
  s.append("text").attr("x", x1).attr("y", top - 22).style("font-weight", 600).style("font-size", "12.5px").text("Artificial tracks found");
  s.append("text").attr("x", x1).attr("y", top - 8).style("font-size", "11px").style("fill", V("--ink-2")).text(`Pass mark: ${Rc.false_alarms} chance candidates across the survey`);
  [["track", V("--ours"), "Track score"], ["vasco", V("--theirs"), "VASCO’s score"]].forEach(([k, c, lab], i) => {
    const pts = Rc.angles.map((a, j) => [a, Rc[k][j]]);
    s.append("path").attr("d", d3.line().x(d => xr(d[0])).y(d => yr(d[1]))(pts)).attr("fill", "none").attr("stroke", c).attr("stroke-width", 2);
    s.append("g").selectAll("circle").data(pts).join("circle").attr("cx", d => xr(d[0])).attr("cy", d => yr(d[1])).attr("r", 3).attr("fill", c);
    const ly = yr(48) + i * 20, lx = xr(30);                          // legend in the empty middle of the panel
    s.append("line").attr("x1", lx).attr("x2", lx + 24).attr("y1", ly).attr("y2", ly).attr("stroke", c).attr("stroke-width", 2);
    s.append("circle").attr("cx", lx + 12).attr("cy", ly).attr("r", 3).attr("fill", c);
    s.append("text").attr("x", lx + 32).attr("y", ly + 4).style("font-size", "11.5px").style("fill", c).style("font-weight", 600).text(lab);
  });
  s.append("text").attr("x", xr(2.5)).attr("y", top + 14).attr("text-anchor", "middle").style("font-size", "10.5px").style("fill", V("--teal")).text("E–W");
  s.append("text").attr("x", xr(87.5)).attr("y", top + 14).attr("text-anchor", "middle").style("font-size", "10.5px").style("fill", V("--teal")).text("N–S");
}

/* ---------- the strongest tracks against the star model of their plates: one row per track with the model, each
   transient on the same scale (and its width as a multiple of the stars'), the light away from the centre relative to
   the stars (1 for star-shaped images; above 1 for wider), and the call; then the typical transient, from transients
   picked at random ---------- */
function drawWidths(D, id) {
  const Z = D.widths, s = svg(id);
  if (!Z) return;
  const W = 790, St = 44, rowH = 92, top = 40, slot = 50, xM = 84, xT = 162, xP = 488, wP = 112, hP = 48, xV = 636;
  const rows = Z.tracks.length + 1, yKey = top + rows * rowH + 8, H = yKey + 40;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const grey = d3.scaleSequential(d3.interpolateGreys).domain([-0.05, 1]).clamp(true);
  // a 13 x 13 stamp, cut to a circle; pixels nearer another source than the transient come as nulls and are left blank
  const pic = (x, y, arr) => {
    const cell = St / 13, g = s.append("g").attr("transform", `translate(${x},${y})`);
    g.selectAll("rect").data(arr.flatMap((row, i) => row.map((v, j) => ({ i, j, v }))).filter(d => d.v != null && Math.hypot(d.i - 6, d.j - 6) <= 6.3)).join("rect")
      .attr("x", d => d.j * cell).attr("y", d => d.i * cell).attr("width", cell + 0.15).attr("height", cell + 0.15).attr("fill", d => grey(d.v));
    g.append("circle").attr("cx", St / 2).attr("cy", St / 2).attr("r", St / 2).attr("fill", "none").attr("stroke", V("--rule"));
  };
  const txt = (x, y, t, o = {}) => s.append("text").attr("x", x).attr("y", y).attr("text-anchor", o.anchor || "start").style("font-size", o.size || "11px")
    .style("font-weight", o.bold ? 600 : null).style("fill", o.col || null).text(t);
  // light away from the centre relative to a star: the band real stars give, the track's band, and its average
  const xr = d3.scaleLinear().domain([0, d3.max(Z.r)]).range([0, wP]), yr = d3.scaleLog().domain([0.6, 3.2]).range([hP, 0]).clamp(true);
  const prof = (y, P) => {
    const g = s.append("g").attr("transform", `translate(${xP},${y})`);
    const band = (lo, hi) => d3.area().x((_, i) => xr(Z.r[i])).y0((_, i) => yr(lo[i])).y1((_, i) => yr(hi[i]))(Z.r);
    g.append("rect").attr("width", wP).attr("height", hP).attr("fill", "none").attr("stroke", V("--rule"));
    g.append("path").attr("d", band(P.stars_lo, P.stars_hi)).attr("fill", V("--ink-3")).attr("fill-opacity", 0.22);
    g.append("line").attr("x1", 0).attr("x2", wP).attr("y1", yr(1)).attr("y2", yr(1)).attr("stroke", V("--ink-3")).attr("stroke-dasharray", "2 3");
    g.append("path").attr("d", band(P.lo, P.hi)).attr("fill", V("--emulsion")).attr("fill-opacity", 0.25);
    g.append("path").attr("d", d3.line().x((_, i) => xr(Z.r[i])).y(v => yr(v))(P.mid)).attr("fill", "none").attr("stroke", V("--emulsion")).attr("stroke-width", 2);
    [1, 2].forEach(v => g.append("text").attr("x", -4).attr("y", yr(v) + 3.5).attr("text-anchor", "end").style("font-size", "9.5px").style("fill", V("--ink-3")).text(`${v}×`));
    g.append("g").attr("transform", `translate(0,${hP})`).call(d3.axisBottom(xr).tickValues([0, 1, 2, 3]).tickFormat(v => `${v}″`).tickSize(2)).call(axisStyle)
      .selectAll("text").style("font-size", "9.5px");
    g.append("text").attr("x", wP / 2).attr("y", hP + 25).attr("text-anchor", "middle").style("font-size", "9.5px").style("fill", V("--ink-3")).text("Distance from centre");
  };
  const odds = p => p < 1 / 20000 ? "under 1 in 20,000" : `1 in ${d3.format(",")(1 / p < 10 ? Math.round(1 / p) : +d3.format(".2r")(1 / p))}`;
  const colour = { wider: V("--emulsion"), narrower: V("--ours"), "like the stars": V("--ink-2") };
  const verdict = (y0, t) => {
    txt(xV, y0 + 18, t.call === "like the stars" ? "Like the stars" : `${t.call[0].toUpperCase()}${t.call.slice(1)} than the stars`, { bold: true, size: "12px", col: colour[t.call] });
    if (t.call === "like the stars") txt(xV, y0 + 34, "Within the range of real stars", { col: V("--ink-2") });
    else {
      txt(xV, y0 + 34, `Real stars spread this ${t.call === "wider" ? "wide" : "narrow"}:`, { col: V("--ink-2") });
      txt(xV, y0 + 48, odds(t.call === "wider" ? t.p_wider : t.p_narrower), { col: V("--ink-2") });
    }
  };
  const width = (x, y, r) => txt(x + St / 2, y + St + 13, `${r.toFixed(1)}×`, { anchor: "middle", size: "10.5px", col: V("--ink-2") });
  [["Star model", xM], ["VASCO transients (width compared to typical star below)", xT], ["Spread relative to a star", xP], ["Verdict", xV]]
    .forEach(([t, x]) => txt(x, top - 18, t, { bold: true, size: "11.5px" }));
  Z.tracks.forEach((t, k) => {
    const y0 = top + k * rowH;
    txt(0, y0 + St / 2 + 4, t.plate, { bold: true, size: "12.5px" });
    pic(xM, y0, t.model);
    t.transients.forEach((r, i) => { pic(xT + i * slot, y0, r.img); width(xT + i * slot, y0, r.ratio); });
    prof(y0 - 2, t.profile);
    verdict(y0, t);
  });
  // the typical transient: the average of transients picked at random, set out in the same way
  const yc = top + Z.tracks.length * rowH + 4, C = Z.control;
  s.append("line").attr("x1", 0).attr("x2", W).attr("y1", yc - 14).attr("y2", yc - 14).attr("stroke", V("--rule"));
  txt(0, yc + 16, "Typical", { bold: true, size: "12.5px" }); txt(0, yc + 31, "transient", { bold: true, size: "12.5px" });
  pic(xM, yc, C.model);
  pic(xT, yc, C.img); width(xT, yc, C.ratio);
  txt(xT + slot + 6, yc + 18, `The average of ${C.n} transients picked`, { col: V("--ink-2") });
  txt(xT + slot + 6, yc + 33, "at random from the catalogue", { col: V("--ink-2") });
  prof(yc - 2, C.profile);
  verdict(yc, C);
  // key for the small plots
  [[V("--ink-3"), 0.22, false, "Real stars of similar brightness (90% range)"], [V("--emulsion"), 0.25, true, "Average over the transients in a track (68% range)"]].forEach(([c, op, line, lab], i) => {
    const k = s.append("g").attr("transform", `translate(${xP - 60},${yKey + i * 18})`);
    k.append("rect").attr("width", 22).attr("height", 10).attr("fill", c).attr("fill-opacity", op);
    if (line) k.append("line").attr("x1", 0).attr("x2", 22).attr("y1", 5).attr("y2", 5).attr("stroke", c).attr("stroke-width", 2);
    k.append("text").attr("x", 30).attr("y", 9).style("font-size", "10.5px").style("fill", V("--ink-2")).text(lab);
  });
}

/* ---------- the trail on XE492, followed across the plate: the strip along the straight line through the track's
   transients (vertical scale stretched), with that line, the transients and the fitted arc (top); the trail's offset
   from the line along its length, against a straight path (bottom) ---------- */
function drawTrail(D, id) {
  const T = D.trail, s = svg(id);
  if (!T) return;
  const W = 760, m = { l: 46, r: 12 }, top = 18, h1 = 120, gap = 46, h2 = 120, H = top + h1 + gap + h2 + 44;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const a0 = T.start_arcmin, x = d3.scaleLinear().domain([T.along0 - a0, T.along1 - a0]).range([m.l, W - m.r]);
  const y1 = d3.scaleLinear().domain([T.across[0], T.across[1]]).range([top + h1, top]);
  s.append("image").attr("href", T.img).attr("x", m.l).attr("y", top).attr("width", W - m.l - m.r).attr("height", h1).attr("preserveAspectRatio", "none");
  s.append("rect").attr("x", m.l).attr("y", top).attr("width", W - m.l - m.r).attr("height", h1).attr("fill", "none").attr("stroke", V("--ink-3"));
  s.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y1(0)).attr("y2", y1(0)).attr("stroke", V("--ours")).attr("stroke-width", 2.4).attr("stroke-dasharray", "6 4");
  const quad = a => { const sd = (a + a0) / 60 + T.s_lo; return T.quad[0] * sd * sd + T.quad[1] * sd + T.quad[2]; };
  const arc = d3.range(0, T.end_arcmin - a0 + 0.01, 0.2).map(a => [a, quad(a)]);
  s.append("g").selectAll("circle").data(T.members_along.map((a, k) => [a - a0, T.members_across[k]])).join("circle")
    .attr("cx", d => x(d[0])).attr("cy", d => y1(d[1])).attr("r", 5).attr("fill", "none").attr("stroke", V("--emulsion")).attr("stroke-width", 1.4);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y1).tickValues([-100, 0, 100]).tickFormat(v => `${v}″`).tickSize(3)).call(axisStyle);
  // offset from the straight line along the trail
  const yA = top + h1 + gap, ext = d3.extent(T.offset.concat([0, 15])), y2 = d3.scaleLinear().domain([Math.min(ext[0], -60) - 5, Math.max(ext[1], 15)]).range([yA + h2, yA]);
  s.append("line").attr("x1", x(0)).attr("x2", x(T.end_arcmin - a0)).attr("y1", y2(0)).attr("y2", y2(0)).attr("stroke", V("--ours")).attr("stroke-width", 2.4).attr("stroke-dasharray", "6 4");
  s.append("text").attr("x", x(T.end_arcmin - a0) - 4).attr("y", y2(0) - 8).attr("text-anchor", "end").style("font-size", "10.5px").style("fill", V("--ours")).text("Track fit to transients");
  s.append("g").selectAll("circle").data(T.s.map((a, k) => [a - a0, T.offset[k]])).join("circle").attr("cx", d => x(d[0])).attr("cy", d => y2(d[1])).attr("r", 2.6).attr("fill", V("--ink-2"));
  s.append("path").attr("d", d3.line().x(d => x(d[0])).y(d => y2(d[1]))(arc)).attr("fill", "none").attr("stroke", V("--emulsion")).attr("stroke-width", 1.6);
  s.append("text").attr("x", x(arc[Math.floor(arc.length * 0.15)][0])).attr("y", y2(arc[Math.floor(arc.length * 0.15)][1]) + 18).style("font-size", "10.5px").style("fill", V("--emulsion"))
    .text(`Arc fit to full trail (${Math.round(T.radius_deg)} degrees)`);
  s.append("g").selectAll("circle.m").data(T.members_along.map((a, k) => [a - a0, T.members_across[k]])).join("circle")
    .attr("cx", d => x(d[0])).attr("cy", d => y2(d[1])).attr("r", 5).attr("fill", "none").attr("stroke", V("--emulsion")).attr("stroke-width", 1.4);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y2).ticks(4).tickFormat(v => `${v}″`).tickSize(3)).call(axisStyle);
  s.append("g").attr("transform", `translate(0,${yA + h2})`).call(d3.axisBottom(x).ticks(9).tickFormat(v => `${v}′`).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").style("font-size", "11.5px")
    .text(`Distance along the trail from where it starts (our scans of this plate end at ${(T.end_arcmin - a0).toFixed(0)}′)`);
}

/* ---------- the tracks that stand out: each field as a blink comparison between the red plate the transients were found
   on, the blue plate taken the same night and the second survey's red plate (one pixel grid, rings and fitted track
   fixed on top), with a second panel where it shows what made the track; findings from the copy ---------- */
function drawTail(D, id) {
  const T = D.tail_tracks, host = document.getElementById(id);
  if (!T || !host) return;
  const notes = {};
  document.querySelectorAll("#tail-findings p").forEach(p => { const m = p.textContent.match(/^(XE\d+):\s*([^|]*)\|\s*(.*)$/); if (m) notes[m[1]] = [m[2].trim(), m[3].trim()]; });
  host.innerHTML = "";
  const fmtOdds = v => (1 / v >= 1e6 ? `1 in ${(1 / v / 1e6).toFixed(1)} million` : `1 in ${d3.format(",")(+d3.format(".2r")(1 / v))}`), S = 240;
  const overlay = (s, xs, ys, line) => {
    const M = xs.map((f, i) => [f * S, ys[i] * S]).filter(d => d[0] >= 0 && d[0] <= S && d[1] >= 0 && d[1] <= S);
    if (line && M.length > 1) {
      const [a, b] = d3.greatest(d3.cross(M, M), ([p, q]) => Math.hypot(p[0] - q[0], p[1] - q[1]));
      const L = Math.hypot(b[0] - a[0], b[1] - a[1]), ux = (b[0] - a[0]) / L, uy = (b[1] - a[1]) / L;
      s.append("line").attr("x1", a[0] - 12 * ux).attr("y1", a[1] - 12 * uy).attr("x2", b[0] + 12 * ux).attr("y2", b[1] + 12 * uy)
        .attr("stroke", V("--ours")).attr("stroke-width", 2.2).attr("stroke-dasharray", "6 4");
    }
    s.append("g").selectAll("circle").data(M).join("circle").attr("cx", d => d[0]).attr("cy", d => d[1]).attr("r", 7)
      .attr("fill", "none").attr("stroke", V("--emulsion")).attr("stroke-width", 1.5);
  };
  const frame = (el, label) => {
    el.append("div").attr("class", "tail-lab").text(label);
    const s = el.append("svg").attr("viewBox", `0 0 ${S} ${S}`);
    return s;
  };
  T.tracks.forEach(t => {
    const card = d3.select(host).append("figure").attr("class", "tail-card"), row = card.append("div").attr("class", "tail-imgs");
    const len = t.length < 1 ? `${Math.round(t.length * 60)}″` : `${t.length.toFixed(1)}′`;
    // blink panel
    const left = row.append("div"), labs = [`Red plate, ${t.years[0]}`, `Blue plate, same night`, `Red plate, ${t.years[2]}`];
    const lab = left.append("div").attr("class", "tail-lab").text(labs[0]);
    const s = left.append("svg").attr("viewBox", `0 0 ${S} ${S}`).attr("role", "img").attr("aria-label", `${t.plate}: the field on three plates`);
    t.blink.forEach(src => { new Image().src = src; });                // preload, so the blink is instant
    const img = s.append("image").attr("href", t.blink[0]).attr("width", S).attr("height", S).attr("preserveAspectRatio", "none");
    s.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", V("--ink-3"));
    overlay(s, t.mx, t.my, true);
    const btns = left.append("div").attr("class", "blink-btns");
    ["Red", "Blue", `${t.years[2]}`].forEach((b, k) => btns.append("button").attr("type", "button").classed("on", !k).text(b)
      .on("click", function () { img.attr("href", t.blink[k]); lab.text(labs[k]); btns.selectAll("button").classed("on", false); d3.select(this).classed("on", true); }));
    // evidence panel, or the caption beside the blink panel
    const right = row.append("div");
    if (t.evidence) {
      const e = t.evidence, es = frame(right, e.kind === "zoom" ? `Red plate, high contrast, ${e.size}′` : `Red plate, ${e.size}′ across`);
      es.append("image").attr("href", e.img).attr("width", S).attr("height", S).attr("preserveAspectRatio", "none");
      es.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", V("--ink-3"));
      overlay(es, e.x, e.y, e.kind === "zoom");
    }
    const cap = (t.evidence ? card : right).append("figcaption");
    cap.append("div").attr("class", "tail-head").html(`<b>${t.plate}</b>`);
    // the track's measures on one line, the checks on its images on the next
    [[`${t.r} transients`, `${len} long`, `Track score ${fmtOdds(t.score)}`, `χ²/ndof ${t.chi2.toFixed(1)}`],
     t.width ? [`Images ${t.width.toFixed(1)}× star width`] : []]
      .filter(line => line.length).forEach(line => { const st = cap.append("div").attr("class", "tail-stats"); line.forEach(v => st.append("span").text(v)); });
  });
}

/* ---------- alignments: VASCO's search and score, real against shuffled copies ---------- */
/* ---------- a score's distribution by candidate size, real catalogue against shifted copies, with the difference in
   standard deviations: the track fit (step 3) and VASCO's own score (step 4) ---------- */
function drawTracks(D) {
  drawDist(D.trackdist, "c-counts", { dom: [-7, 0], step: 0.2, ticks: [-6, -4, -2, 0], fmt: { "-6": "10⁻⁶", "-4": "10⁻⁴", "-2": "0.01", "0": "1" },
    xlab: "Track score: expected number of such tracks by chance (lower is less likely)", what: "track score",
    rowlab: "Number of track candidates at each score", floor: 0.01, yticks: [0.01, 1, 100] });
}
function drawDist(T, id, o) {
  const s = svg(id), W = 760, m = { l: 52, r: 10, t: 44, b: 52 }, gx = 12;
  const rows = [{ key: "dist", h: 160, lab: o.rowlab }, { key: "all", h: 110, lab: "Difference between real catalogue and random chance (standard deviations)" }];
  const gy = 34, H = m.t + d3.sum(rows, r => r.h) + gy * (rows.length - 1) + m.b;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const pw = (W - m.l - m.r - 3 * gx) / 4, step = o.step, lo = T.lo;
  const fmtX = v => o.fmt[String(v)];
  let top = m.t;
  rows.forEach((row, ri) => {
    const bot = top + row.h;
    s.append("text").attr("class", "ann").attr("x", m.l).attr("y", top - 8).style("font-weight", 600).text(row.lab);
    [0, 1, 2, 3].forEach(ci => {
      const x0 = m.l + ci * (pw + gx), x = d3.scaleLinear().domain(o.dom).range([x0, x0 + pw]).clamp(true);
      const d = T.all[ci];
      if (ri === 0) s.append("text").attr("x", x0 + pw / 2).attr("y", 14).attr("text-anchor", "middle").style("font-weight", 600).text(`${d.r}-point candidates`);
      s.append("rect").attr("x", x0).attr("y", top).attr("width", pw).attr("height", row.h).attr("fill", "none").attr("stroke", V("--rule"));
      const bins = lo.map((l, k) => ({ l, h: l + step, n: d.real[k], b: d.bmean[k], sd: d.bsd[k] })).filter(q => q.l >= o.dom[0] - 1e-9 && q.h <= o.dom[1] + 1e-9);
      if (row.key === "dist") {
        const fl = o.floor || 0.5, y = d3.scaleLog().domain([fl, d3.max(bins, q => Math.max(q.n, q.b)) * 3]).range([bot, top]).clamp(true);
        const stepPath = acc => d3.line().x(q => x(q[0])).y(q => y(Math.max(fl, q[1]))).curve(d3.curveStepAfter)(bins.map(q => [q.l, acc(q)]).concat([[bins[bins.length - 1].h, acc(bins[bins.length - 1])]]));
        s.append("path").attr("d", stepPath(q => q.b) + `L${x(bins[bins.length - 1].h)},${bot}L${x(bins[0].l)},${bot}Z`).attr("fill", V("--ink-3")).attr("fill-opacity", 0.35);
        if (o.mark) {
          s.append("line").attr("x1", x(o.mark.at)).attr("x2", x(o.mark.at)).attr("y1", top).attr("y2", bot).attr("stroke", V("--theirs")).attr("stroke-dasharray", "3 3");
          if (ci === 0) s.append("text").attr("x", x(o.mark.at) + 4).attr("y", top + 12).style("font-size", "10.5px").style("fill", V("--theirs")).text(o.mark.lab);
        }
        const pts = bins.filter(q => q.n > 0);
        // error bars: Poisson 68% interval (Gehrels 1986)
        s.append("g").selectAll("line").data(pts).join("line").attr("x1", q => x(q.l + step / 2)).attr("x2", q => x(q.l + step / 2))
          .attr("y1", q => y(Math.max(fl, q.n * (1 - 1 / (9 * q.n) - 1 / (3 * Math.sqrt(q.n))) ** 3))).attr("y2", q => y(q.n + 1 + Math.sqrt(q.n + 0.75))).attr("stroke", V("--emulsion")).attr("stroke-width", 1.4);
        s.append("g").selectAll("circle").data(pts).join("circle").attr("cx", q => x(q.l + step / 2)).attr("cy", q => y(q.n)).attr("r", 2.6).attr("fill", V("--emulsion"));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues(o.yticks || [1, 100, 10000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
      } else {
        const y = d3.scaleLinear().domain([-4, 4]).range([bot, top]).clamp(true);
        for (const [k, op] of [[2, 0.18], [1, 0.3]]) s.append("rect").attr("x", x0).attr("width", pw).attr("y", y(k)).attr("height", y(-k) - y(k)).attr("fill", V("--ink-3")).attr("fill-opacity", op);
        s.append("line").attr("x1", x0).attr("x2", x0 + pw).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", V("--ink-3")).attr("stroke-dasharray", "3 3");
        const ok = bins.filter(q => q.sd > 0 && q.b >= 1);   // pulls mean nothing where chance expects under one group
        s.append("g").selectAll("circle").data(ok).join("circle").attr("class", "hov").attr("cx", q => x(q.l + step / 2)).attr("cy", q => y((q.n - q.b) / q.sd)).attr("r", 2.8).attr("fill", V("--emulsion"))
          .call(sel => hov(sel, q => `${d.r}-point candidates, ${o.what} ${(10 ** q.l).toPrecision(2)} to ${(10 ** q.h).toPrecision(2)}<br>Real catalogue ${d3.format(",")(q.n)}; chance ${d3.format(",.1f")(q.b)} ± ${q.sd.toFixed(1)}<br>${((q.n - q.b) / q.sd).toFixed(1)} standard deviations`));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([-2, 0, 2, 4]).tickFormat(v => (v > 0 ? "+" : "") + v + "σ").tickSize(3)).call(axisStyle);
      }
      if (ri === rows.length - 1) s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).tickValues(o.ticks).tickFormat(fmtX).tickSize(3).tickSizeOuter(0)).call(axisStyle);
    });
    top = bot + gy;
  });
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 12).attr("text-anchor", "middle").text(o.xlab);
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
function drawRadialMap(D) {
  const M = D.radial.map, s = svg("c-radial-map"), S = 300, top = 36, xs = [44, 380], n = M.transients.length, c = S / n;
  const col = d3.scaleSequential(d3.interpolateMagma).domain([0, 3.5]).clamp(true);
  const ax = d3.scaleLinear().domain([M.edges[0], M.edges[n]]).range([0, S]);
  [["transients", "Transients"], ["control", "VASCO’s random positions"]].forEach(([k, lab], j) => {
    const g = s.append("g").attr("transform", `translate(${xs[j]},${top})`);
    g.append("text").attr("x", S / 2).attr("y", -12).attr("text-anchor", "middle").style("font-weight", 600).text(lab);
    // rows run south to north, columns west to east in the data; drawn north up, east left
    g.append("g").selectAll("rect").data(M[k].flatMap((row, i) => row.map((v, jj) => ({ v, i, j: jj })))).join("rect").attr("class", "hov")
      .attr("x", d => S - (d.j + 1) * c).attr("y", d => S - (d.i + 1) * c).attr("width", c + 0.3).attr("height", c + 0.3).attr("fill", d => col(d.v))
      .call(sel => hov(sel, d => `${d.v.toFixed(1)} times the density near the centre`));
    g.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", V("--rule"));
    g.append("g").attr("transform", `translate(0,${S})`).call(d3.axisBottom(ax).tickValues([-3, 0, 3]).tickFormat(d => `${d3.format("")(d)}°`).tickSize(3)).call(axisStyle);
    g.append("text").attr("x", -6).attr("y", 12).attr("text-anchor", "end").attr("fill", V("--ink-3")).text("N");
    g.append("text").attr("x", -6).attr("y", S / 2 + 4).attr("text-anchor", "end").attr("fill", V("--ink-3")).text("E");
  });
  // colour bar
  const bx = 716, bh = S, gid = "radial-map-grad", ly = d3.scaleLinear().domain([0, 3.5]).range([top + bh, top]);
  s.append("defs").append("linearGradient").attr("id", gid).attr("x1", 0).attr("x2", 0).attr("y1", 1).attr("y2", 0)
    .selectAll("stop").data(d3.range(0, 1.001, 0.1)).join("stop").attr("offset", t => t).attr("stop-color", t => col(t * 3.5));
  s.append("rect").attr("x", bx).attr("y", top).attr("width", 12).attr("height", bh).attr("fill", `url(#${gid})`);
  s.append("g").attr("transform", `translate(${bx + 12},0)`).call(d3.axisRight(ly).tickValues([0, 1, 2, 3]).tickFormat(d => `${d}×`).tickSize(3)).call(axisStyle);
  s.append("text").attr("transform", `translate(${bx - 8},${top + bh / 2}) rotate(-90)`).attr("text-anchor", "middle").attr("fill", V("--ink-2")).text("Density relative to the centre");
}
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
  s.append("g").selectAll("circle").data(bins).join("circle").attr("cx", d => x(d.b)).attr("cy", d => y(d.q50)).attr("r", 3.5).attr("fill", V("--ours"))
    .attr("stroke", V("--surface")).attr("stroke-width", 6).attr("stroke-opacity", 0.01)
    .call(sel => hov(sel, d => `${Math.max(0, d.b - half)}° to ${d.b + half}° from the Milky Way: ${d.n} plates<br>Median ${Math.round(d.q50)} transients; middle half ${Math.round(d.q25)} to ${Math.round(d.q75)}`));
}

/* ---------- nuclear: VASCO's random-night test, the nights near tests against random sets of nights ---------- */
function drawPerm(D, id, key, lab, o = {}) {
  const q = D.perm[key], s = svg(id), W = 760, H = 250, m = { l: 52, r: 16, t: 34, b: 44 }, bw = 0.05;
  const line = o.line || "Nights near tests", trials = o.trials || "Random sets", xlab = o.xlab || "Transients on the chosen nights relative to all other nights";
  const x = d3.scaleLinear().domain([0, 3]).range([m.l, W - m.r]);
  const top = m.t, bot = H - m.b, tot = d3.sum(q.hist);
  const y = d3.scaleLinear().domain([0, Math.max(0.12, d3.max(q.hist) / tot * 1.05)]).range([bot, top + 8]);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0, 0.05, 0.1]).tickFormat(d3.format(".0%")).tickSize(3)).call(axisStyle);
  // the bars, and the part of them at or beyond the real value drawn darker, split exactly at the red line
  const bars = (g, op) => g.selectAll("rect").data(q.hist).join("rect").attr("x", (c, j) => x(j * bw) + 0.75).attr("width", x(bw) - x(0) - 1.5)
    .attr("y", c => y(c / tot)).attr("height", c => bot - y(c / tot)).attr("fill", V("--ink-3")).attr("fill-opacity", op);
  bars(s.append("g"), 0.45);
  s.append("clipPath").attr("id", `${id}-tail`).append("rect").attr("x", x(q.real)).attr("width", W - x(q.real)).attr("y", 0).attr("height", H);
  bars(s.append("g").attr("clip-path", `url(#${id}-tail)`), 0.9);
  s.append("line").attr("x1", x(q.real)).attr("x2", x(q.real)).attr("y1", top).attr("y2", bot).attr("stroke", V("--emulsion")).attr("stroke-width", 2.5);
  s.append("text").attr("class", "ann").attr("x", m.l).attr("y", 16).style("font-weight", 600).text(lab);
  const sd = `${q.sigma.toFixed(1)} standard deviation${q.sigma.toFixed(1) === "1.0" ? "" : "s"}`;
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 14).style("fill", V("--emulsion")).text(`${line}: ${q.real.toFixed(2)}×`);
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 30).style("fill", V("--emulsion")).text(`${trials} reach it 1 time in ${Math.round(1 / q.p)} (${sd})`);
  s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).ticks(6).tickFormat(d => d + "×").tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text(xlab);
  s.append("text").attr("transform", `translate(12,${(top + bot) / 2}) rotate(-90)`).attr("text-anchor", "middle").text(o.ylab || "Share of random sets");
}

/* ---------- nuclear: how far each further check moves the ratio per photograph allowing for pointing, as bars either
   side of it ---------- */
function drawNucChecks(D, id, key = "nuc_checks") {
  const C = D[key], s = svg(id); if (!C) return;
  const W = 760, rowH = 26, m = { l: 270, r: 90, t: 40, b: 50 }, gap = 14, base = C.base;
  let yy = m.t; const pos = [];
  C.rows.forEach((r, i) => { if (i === 0 || C.rows[i - 1].group !== r.group) { yy += i ? gap : 0; pos.push({ head: r.group, y: yy }); yy += 20; } pos.push({ row: r, y: yy }); yy += rowH; });
  const H = yy + m.b; s.attr("viewBox", `0 0 ${W} ${H}`);
  const x = d3.scaleLinear().domain([-C.span, C.span]).range([m.l, W - m.r]).clamp(true);
  const sg = v => { v = Math.round(v * 100) / 100; return (v > 0 ? "+" : v < 0 ? "−" : "±") + Math.abs(v).toFixed(2); };
  s.append("g").attr("transform", `translate(0,${H - m.b + 8})`).call(d3.axisBottom(x).tickValues(d3.range(-C.span, C.span + 1e-9, 0.1)).tickFormat(sg).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text(`Change from ${base.toFixed(2)}×`);
  s.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", m.t - 10).attr("y2", H - m.b + 8).attr("stroke", V("--emulsion")).attr("stroke-width", 1.5);
  s.append("text").attr("class", "ann").attr("x", x(0)).attr("y", m.t - 16).attr("text-anchor", "middle").style("fill", V("--emulsion"))
    .text(`${C.title}: ${base.toFixed(2)}×`);
  pos.forEach(p_ => {
    if (p_.head) { s.append("text").attr("x", 10).attr("y", p_.y + 6).style("font-weight", 600).text(p_.head); return; }
    const r = p_.row, d = r.r - base;
    s.append("text").attr("x", m.l - 12).attr("y", p_.y + 4).attr("text-anchor", "end").style("fill", V("--ink-2")).text(r.label);
    s.append("rect").attr("x", Math.min(x(0), x(d))).attr("width", Math.max(1.5, Math.abs(x(d) - x(0)))).attr("y", p_.y - 7).attr("height", 14)
      .attr("fill", V("--emulsion")).attr("fill-opacity", 0.8);
    s.append("text").attr("class", "ann").attr("x", W - m.r + 12).attr("y", p_.y + 4).style("fill", V("--ink-2"))
      .text(sg(d));
  });
}

/* ---------- nuclear: the result for every window that includes the test day (up to a week each side), as a histogram;
   VASCO's window marked ---------- */
function drawNucWindowHist(D, id) {
  const F = D.nuc_windows, s = svg(id); if (!F) return;
  const R = F.cells.map(d => d.r), V_ = F.cells.find(d => d.before === 1 && d.after === 1);
  const W = 760, H = 290, m = { l: 52, r: 16, t: 40, b: 50 }, bw = 0.025;
  const lo = Math.floor(d3.min(R) / bw) * bw, hi = Math.ceil(d3.max(R) / bw) * bw;
  const bins = d3.bin().domain([lo, hi]).thresholds(d3.range(lo, hi + 1e-9, bw))(R);
  const x = d3.scaleLinear().domain([Math.min(0.85, lo), Math.max(1.45, hi)]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, d3.max(bins, b => b.length) + 1]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(7).tickFormat(d => `${d}×`).tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickFormat(d3.format("d")).tickSize(3)).call(axisStyle);
  s.append("g").selectAll("rect").data(bins).join("rect").attr("x", b => x(b.x0) + 1).attr("width", b => Math.max(0, x(b.x1) - x(b.x0) - 2))
    .attr("y", b => y(b.length)).attr("height", b => y(0) - y(b.length)).attr("rx", 1.5).attr("fill", V("--ink-3")).attr("fill-opacity", 0.6);
  s.append("line").attr("x1", x(1)).attr("x2", x(1)).attr("y1", m.t - 8).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 1.2);
  s.append("text").attr("class", "ann").attr("x", x(1)).attr("y", m.t - 14).attr("text-anchor", "middle").style("fill", V("--ink-2")).text("No difference");
  s.append("line").attr("x1", x(V_.r)).attr("x2", x(V_.r)).attr("y1", m.t - 8).attr("y2", H - m.b).attr("stroke", V("--emulsion")).attr("stroke-width", 2.4);
  s.append("text").attr("class", "ann").attr("x", x(V_.r) + 8).attr("y", m.t + 4).style("fill", V("--emulsion")).style("font-weight", 600).text(`VASCO’s window: ${V_.r.toFixed(2)}×`);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 10).attr("text-anchor", "middle").text("Transients near tests relative to the rest, for each window");
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Windows");
}

let DATA = null;
function render() {
  if (!DATA) return;
  drawSky(DATA, "c-sky0", false); drawSpread(DATA); drawSky(DATA, "c-sky", true); drawChance(DATA); drawPerm(DATA, "c-perm-any", "raw", "VASCO’s test");
  drawPerm(DATA, "c-nuc-corr", "corrected", "Per photograph, allowing for pointing", { line: "Plates near tests", trials: "Shifted calendars",
    xlab: "Transients per photograph near the chosen dates, against pointing, relative to the rest", ylab: "Share of shifted calendars" }); drawNucChecks(DATA, "c-nuc-checks"); drawNucWindowHist(DATA, "c-nuc-windows");
  drawFan(DATA); drawTracks(DATA); drawCompare(DATA, "c-compare"); drawTail(DATA, "tail-gallery"); drawWidths(DATA, "c-widths"); drawTrail(DATA, "c-trail"); fillDates(DATA); drawRadialMap(DATA); drawRadial(DATA); drawCrowd(DATA);
  drawPlates(DATA, "c-plates"); drawForest(DATA, "c-ratio-vasco"); drawRatioNull(DATA, "c-ratio-null");
}
fetch("data/site.json").then(r => r.json()).then(d => { DATA = d; render(); })
  .catch(() => document.querySelectorAll(".chart svg").forEach(s => s.insertAdjacentHTML("afterend", '<div class="loading">The chart data could not be loaded.</div>')));
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
