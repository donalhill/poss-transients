const tip = document.getElementById("tip");
const sdText = v => v.toFixed(1) === "1.0" ? "1 standard deviation" : `${v.toFixed(1)} standard deviations`;   // exactly one reads as such
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
  if (shadow) {                                                                  // the ecliptic, along which the shadow moves through the year
    const ecl = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => { const e = 23.44 * DEG, L = l * DEG;
      return [((Math.atan2(Math.sin(L) * Math.cos(e), Math.cos(L)) / DEG) + 360) % 360, Math.asin(Math.sin(e) * Math.sin(L)) / DEG]; }) };
    s.append("path").datum(ecl).attr("d", path).attr("fill", "none").attr("stroke", SHADOW).attr("stroke-opacity", 0.55).attr("stroke-width", 1.2)
      .attr("stroke-dasharray", "1 4").attr("stroke-linecap", "round").attr("pointer-events", "none");
  }
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "shp").attr("pointer-events", "none").attr("d", p => path(p.g)).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-width", 5);
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "hov shp").attr("d", p => path(p.g))
    .attr("fill", p => p.N !== undefined ? col(Math.max(p.N, 5)) : V("--surface")).attr("stroke", "#111418").attr("stroke-width", 1.8)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients on the whole plate<br>${Math.round(p.frac * 100)}% of the central 6 × 6° inside the shadow<br>${p.obs} transients in the shadowed part; an average plate would have about ${Math.round(p.exp)} there`));

  s.append("g").attr("pointer-events", "none").selectAll("path").data(sh).join("path").attr("class", "shp").attr("d", p => path(ring(p.circle)))
    .attr("fill", "none").attr("stroke", SHADOW).attr("stroke-width", 2.4);
  if (shadow) playShadow(s, sh, W);
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

/* ---------- the shadow plates one at a time, in the order they were taken: each plate and its shadow light up in turn,
   with a label, then all are shown together ---------- */
function playShadow(s, sh, W) {
  const order = sh.slice().sort((a, b) => a.date.localeCompare(b.date));
  // in the top-right corner, clear of the disc and its hour labels
  const lab = s.append("text").attr("x", W - 4).attr("y", 18).attr("text-anchor", "end").style("font-size", "15px").style("font-weight", 600).style("fill", V("--ink"));
  const sub = s.append("text").attr("x", W - 4).attr("y", 38).attr("text-anchor", "end").style("font-size", "14px").style("fill", V("--ink-2"));
  const day = d3.timeFormat("%b %Y"), all = s.selectAll(".shp");
  d3.select("#shadow-play").on("click", () => {
    const step = matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1300;
    all.interrupt().style("opacity", 0); lab.text(""); sub.text("");
    order.forEach((p, i) => {
      setTimeout(() => {
        all.filter(q => q === p).transition().duration(step ? 400 : 0).style("opacity", 1);
        all.filter(q => q !== p && order.indexOf(q) < i).transition().duration(step ? 400 : 0).style("opacity", 0.35);
        lab.text(`${p.id}, ${day(new Date(p.date))}`);
        sub.text(`${p.obs} in shadow, ${Math.round(p.exp)} expected`);
      }, i * step);
    });
    setTimeout(() => { all.transition().duration(step ? 400 : 0).style("opacity", 1); lab.text(""); sub.text(""); }, order.length * step + (step ? 600 : 0));
  });
}

/* ---------- how many transients a plate holds ---------- */
function drawSpread(D) {
  const s = svg("c-spread"), W = 760, H = 300, m = { l: 48, r: 16, t: 26, b: 44 };
  const N = D.plates.map(p => p.N), mu = d3.mean(N), sd = Math.sqrt(mu), lo = mu - 2 * sd, hi = mu + 2 * sd;
  const bins = d3.range(0, 1000, 25).map(a => ({ lo: a, hi: a + 25, n: N.filter(v => v >= a && (a === 975 ? v <= 1000 : v < a + 25)).length }));
  const over = { lo: 1025, hi: 1075, n: N.filter(v => v > 1000).length, more: true };
  const inBand = N.filter(v => v >= lo && v <= hi).length;
  const x = d3.scaleLinear().domain([0, 1080]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, Math.ceil(d3.max(bins.concat(over), b => b.n) / 10) * 10]).range([H - m.b, m.t]);
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
  // The observed count of transients in shadow against two measures of chance, told as slides the reader steps through:
  // 1 the observed count; 2 VASCO's model (every transient independent), its expected value and standard deviation;
  // 3 the shortfall in those standard deviations; 4 the real spread, built trial by trial from random plates; 5 the trials
  // as low as the observed count, shaded exactly (a clip at the observed count splits the bar it falls in); and a summary,
  // shown on load. Each slide declares its whole end state: what leaves fades out first, then the axis moves, then what
  // arrives fades in, so nothing from one slide shows on another
  const P = D.paste, s = svg("c-chance"), W = 760, H = 402, m = { l: 44, r: 16, t: 100, b: 42 };
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const binw = P.bins[1] - P.bins[0], bins = P.counts.map((c, i) => ({ x0: P.bins[i], x1: P.bins[i + 1], c })).filter(d => d.x0 < 1600);
  const mu = P.area_expected, sd = P.binom_sd, obsN = P.observed, wide = [0, 1600], near = [mu - 10.5 * sd, mu + 4.5 * sd];
  const x = d3.scaleLinear().domain(wide).range([m.l, W - m.r]), yb = d3.scaleLinear().domain([0, d3.max(bins, d => d.c) * 1.05]).range([H - m.b, m.t]);
  const yn = d3.scaleLinear().domain([0, 1.04]).range([H - m.b, m.t]);
  const gauss = d3.range(mu - 4.5 * sd, mu + 4.5 * sd + 0.1, sd / 20).map(v => [v, Math.exp(-0.5 * ((v - mu) / sd) ** 2)]);
  const ink = V("--ink"), ink2 = V("--ink-2"), pur = V("--naive"), blue = V("--ours"), fmtN = d3.format(","), pct = Math.round(P.p * 100);
  const axis = s.append("g").attr("transform", `translate(0,${H - m.b})`);
  const drawAxis = dur => axis.transition().duration(dur).call(d3.axisBottom(x).ticks(8).tickFormat(fmtN).tickSizeOuter(0)).on("end", function () { d3.select(this).call(axisStyle); });
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Number of transients inside the shadow");
  const t1 = s.append("text").attr("x", m.l).attr("y", 20).style("font-size", "13.5px"), t2 = s.append("text").attr("x", m.l).attr("y", 42).style("font-size", "13.5px");
  // the real spread: light bars, and a dark copy clipped at the observed count
  const clip = s.append("defs").append("clipPath").attr("id", "chance-clip").append("rect").attr("x", 0).attr("y", 0).attr("height", H).attr("width", 0);
  const light = s.append("g").attr("fill", blue).attr("fill-opacity", 0.32), dark = s.append("g").attr("fill", blue).attr("fill-opacity", 0.95).attr("clip-path", "url(#chance-clip)");
  const barsAt = (g, h) => g.selectAll("rect").data(bins).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1)).attr("rx", 1.5)
    .attr("y", (d, i) => yb(h[i])).attr("height", (d, i) => yb(0) - yb(h[i]));
  const zero = bins.map(() => 0), full = bins.map(d => d.c);
  // VASCO's model: the curve; its expected value and one standard deviation either side; the shortfall bracket
  const curve = s.append("path").attr("stroke", pur).attr("stroke-width", 1.8).attr("fill", pur);
  const sigma = s.append("path").attr("fill", pur).attr("fill-opacity", 0.32);
  const muLine = s.append("line").attr("stroke", pur).attr("stroke-dasharray", "4 3");
  const muLab = s.append("text").attr("text-anchor", "middle").style("font-size", "12px").style("fill", pur);
  const sdLab = s.append("text").style("font-size", "12px").style("fill", pur);
  const bracket = s.append("g").attr("stroke", pur).attr("fill", "none");
  const brLab = s.append("text").attr("text-anchor", "middle").style("font-size", "12px").style("font-weight", 600).style("fill", pur);
  const placeNaive = () => {                                                      // geometry only, for the current axis
    const area = d3.area().x(d => x(d[0])).y0(yn(0)).y1(d => yn(d[1])).curve(d3.curveBasis);
    curve.attr("d", area(gauss)); sigma.attr("d", area(gauss.filter(d => Math.abs(d[0] - mu) <= sd)));
    muLine.attr("x1", x(mu)).attr("x2", x(mu)).attr("y1", yn(0)).attr("y2", yn(1.12)); muLab.attr("x", x(mu)).attr("y", yn(1.12) - 6).text(`Expected ${Math.round(mu)}`);
    sdLab.attr("x", x(mu + sd) + 6).attr("y", yn(0.62)).text(`±${Math.round(sd)}`);
    const yB = yn(0.72); bracket.selectAll("*").remove();
    bracket.append("path").attr("d", `M${x(obsN)},${yB}H${x(mu)}`).attr("stroke-width", 1.5);
    d3.range(0, Math.floor((mu - obsN) / sd) + 1).forEach(k => bracket.append("path").attr("d", `M${x(mu - k * sd)},${yB - 4}v8`).attr("stroke-width", 1.2));
    brLab.attr("x", (x(obsN) + x(mu)) / 2).attr("y", yB - 10).text(`${sdText(P.binom_z)}`);
  };
  // the observed count, always shown
  const obsLine = s.append("line").attr("stroke", ink).attr("stroke-width", 2), obsLab = s.append("text").attr("text-anchor", "middle").style("font-size", "12px").style("font-weight", 600).text(`Observed: ${obsN}`);
  const placeObs = dur => { obsLine.transition().duration(dur).attr("x1", x(obsN)).attr("x2", x(obsN)).attr("y1", m.t - 10).attr("y2", H - m.b); obsLab.transition().duration(dur).attr("x", x(obsN)).attr("y", m.t - 14); };
  const parts = { curve: [curve], sigma: [sigma, muLine, muLab, sdLab], bracket: [bracket, brLab] };
  Object.values(parts).flat().forEach(e => e.style("opacity", 0));               // every part starts hidden
  const slides = [
    { dom: near, curve: 0, sigma: 0, bracket: 0, bars: "none", shade: "none", text: [[`The shadow holds ${obsN} transients.`, ink, 400]] },
    { dom: near, curve: 1, sigma: 1, bracket: 0, bars: "none", shade: "none",
      text: [[`If every transient were independent, as in VASCO’s model, the shadow should hold ${Math.round(mu)} ± ${Math.round(sd)} transients.`, pur, 400]] },
    { dom: near, curve: 1, sigma: 1, bracket: 1, bars: "none", shade: "none",
      text: [[`The shortfall is then ${sdText(P.binom_z)}: far beyond chance, on this model.`, pur, 400]] },
    { dom: wide, curve: 0, sigma: 0, bracket: 0, bars: "build", shade: "none",
      text: [["But plates differ far more than counting allows. So we lay the real shadows on random plates, trial after trial.", ink, 400]] },
    { dom: wide, curve: 0, sigma: 0, bracket: 0, bars: "full", shade: "sweep",
      text: [[`${pct}% of the trials are as low as observed: ${sdText(P.sig)}.`, blue, 600]] },
    { dom: wide, curve: 0.12, sigma: 0, bracket: 0, bars: "full", shade: "full",
      text: [[`If every transient were independent (purple): ${sdText(P.binom_z)}`, pur, 400],
             [`Laid on random plates (blue): ${pct}% of trials as low, ${sdText(P.sig)}`, blue, 600]] }
  ];
  const last = slides.length - 1, op = { curve: 0, sigma: 0, bracket: 0 };      // what each part currently shows
  let cur = -1, timers = [];
  const at = (ms, f) => timers.push(setTimeout(f, ms));
  const setText = (sl, dur) => [t1, t2].forEach((el, i) => { const tx = sl.text[i];
    el.interrupt().style("opacity", 0).text(tx ? tx[0] : "").style("fill", tx ? tx[1] : ink).style("font-weight", tx ? tx[2] : 400).transition().duration(dur).style("opacity", 1); });
  const tally = k => t2.interrupt().style("opacity", 1).style("fill", ink2).style("font-weight", 400).text(`${fmtN(k)} trials`);
  const show = (k, dur) => {
    timers.forEach(clearTimeout); timers = []; s.selectAll("*").interrupt();
    const sl = slides[k], moved = sl.dom !== x.domain() && (sl.dom[0] !== x.domain()[0] || sl.dom[1] !== x.domain()[1]);
    const out = dur ? 250 : 0, move = moved && dur ? 800 : 0;
    // 1. leave: parts not on the new slide (or on a moving axis) fade out; the shading and texts go at once
    for (const [key, els] of Object.entries(parts)) if (op[key] && (!sl[key] || moved)) { els.forEach(e => e.transition().duration(out).style("opacity", 0)); op[key] = 0; }
    t1.style("opacity", 0); t2.style("opacity", 0);
    if (sl.shade !== "full" && sl.shade !== "sweep") clip.attr("width", 0);
    if (moved || sl.bars === "none") { barsAt(light, sl.bars === "none" || moved ? zero : full); barsAt(dark, sl.bars === "none" || moved ? zero : full); }
    // 2. move the axis
    at(out, () => { x.domain(sl.dom); drawAxis(move); placeObs(move); });
    // 3. arrive: the new parts at their places, then faded in; the bars and shading as the slide asks
    at(out + move, () => {
      placeNaive(); curve.attr("fill-opacity", sl.curve === 1 ? 0.22 : 0.12).attr("stroke-opacity", sl.curve === 1 ? 1 : 0.8);
      for (const [key, els] of Object.entries(parts)) if (sl[key] && op[key] !== sl[key]) { els.forEach(e => e.style("opacity", 0).transition().duration(dur ? 450 : 0).style("opacity", 1)); op[key] = sl[key]; }
      setText(sl, dur ? 400 : 0);
      if (sl.bars === "full") { barsAt(light, full); barsAt(dark, full); }
      if (sl.bars === "build") {
        barsAt(dark, full); clip.attr("width", 0);
        if (!dur) { barsAt(light, full); return tally(P.n); }
        const cum = d3.cumsum(full), tot = cum[cum.length - 1]; let seed = 7; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
        const h = bins.map(() => 0); let done = 0; barsAt(light, h);
        for (let k2 = 1; k2 <= 48; k2++) at(400 + k2 * 60, () => { const target = Math.round(P.n * (k2 / 48) ** 2.2);
          for (; done < target; done++) h[d3.bisectRight(cum, rnd() * tot)] += 1; barsAt(light, h); tally(done); });
        at(400 + 49 * 60, () => { barsAt(light, full); tally(P.n); });
      }
      if (sl.shade === "full") clip.attr("width", x(obsN));
      if (sl.shade === "sweep") clip.attr("width", x(0)).transition().duration(dur ? 1200 : 0).ease(d3.easeCubicInOut).attr("width", x(obsN));
    });
    cur = k;
  };
  // the controls
  const ctl = d3.select("#chance-steps"), back = ctl.select(".back"), next = ctl.select(".next"), pos = ctl.select(".pos");
  const go = (k, instant) => {
    show(k, instant || matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1);
    back.property("disabled", k === 0);
    next.text(k === last ? "Step through ›" : k === last - 1 ? "Summary ›" : "Next ›");
    pos.text(k === last ? "Summary" : `${k + 1} of ${last}`);
  };
  back.on("click", () => go(cur - 1));
  next.on("click", () => go(cur === last ? 0 : cur + 1));
  x.domain(wide); drawAxis(0); placeObs(0); go(last, true);                      // the summary at once on load
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
  // slides the reader steps through: 1 the fit on trials with no shadow effect, built up, with its middle 68%; 2 the fit
  // on the real plates; 3 the same trials with a shortfall as large as VASCO's put in, and the share of them that still
  // fit as high as the real plates; and a summary, shown on load. Each slide declares its end state
  const Z = D.shadow, s = svg(id);
  if (!Z || !Z.toys) return;
  const R = Z.toys, I = Z.injected, W = 760, H = 384, m = { l: 20, r: 20, t: 116, b: 46 }, bw = 3 / R.hist_null.length;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const keep = (h) => h.map((c, i) => ({ x0: i * bw, x1: (i + 1) * bw, c })).filter(d => d.x0 >= 0.2 - 1e-9 && d.x1 <= 2.6 + 1e-9);
  const B0 = keep(R.hist_null), B1 = keep(I.hist), n0 = d3.sum(R.hist_null), r = Z.fit.r;
  const x = d3.scaleLinear().domain([0.2, 2.6]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, d3.max(B0.concat(B1), d => d.c) * 1.08]).range([H - m.b, m.t]);
  const ink = V("--ink"), ink2 = V("--ink-2"), grey = V("--ink-3"), blue = V("--ours"), pur = V("--naive"), fmtN = d3.format(",");
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(12).tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text(`Best fit ratio (r) across ${Z.plates.length} plates`);
  const t1 = s.append("text").attr("x", m.l).attr("y", 20).style("font-size", "13.5px"), t2 = s.append("text").attr("x", m.l).attr("y", 42).style("font-size", "13.5px");
  const bars = (g, B, h) => g.selectAll("rect").data(B).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1)).attr("rx", 1.5)
    .attr("y", (d, i) => y(h[i])).attr("height", (d, i) => y(0) - y(h[i]));
  // the null: grey bars, its middle 68% marked beneath
  const nullG = s.append("g").attr("fill", grey).attr("fill-opacity", 0.55);
  const band = s.append("g"), [lo, hi] = [R.null.q16, R.null.q84], yB = H - m.b + 0;
  band.append("rect").attr("x", x(lo)).attr("width", x(hi) - x(lo)).attr("y", m.t - 6).attr("height", H - m.b - m.t + 6).attr("fill", grey).attr("fill-opacity", 0.1);
  band.append("text").attr("x", (x(lo) + x(hi)) / 2).attr("y", m.t - 10).attr("text-anchor", "middle").style("font-size", "12px").style("fill", ink2).text("Middle 68% of trials");
  // the injected shortfall: purple outline bars, with those as high as the real plates shaded
  const injG = s.append("g").attr("fill", pur).attr("fill-opacity", 0.14).attr("stroke", pur).attr("stroke-opacity", 0.6);
  const cid = `${id}-clip`, clip = s.append("defs").append("clipPath").attr("id", cid).append("rect").attr("y", 0).attr("height", H).attr("x", x(r)).attr("width", 0);
  const injHi = s.append("g").attr("fill", pur).attr("fill-opacity", 0.75).attr("clip-path", `url(#${cid})`);
  const injLab = s.append("text").attr("x", x(I.r_true)).attr("text-anchor", "middle").style("font-size", "12px").style("fill", pur);
  // the real plates
  const obs = s.append("g");
  obs.append("line").attr("x1", x(r)).attr("x2", x(r)).attr("y1", m.t - 26).attr("y2", H - m.b).attr("stroke", blue).attr("stroke-width", 2.4);
  obs.append("text").attr("x", x(r)).attr("y", m.t - 30).attr("text-anchor", "middle").style("font-size", "12px").style("font-weight", 600).style("fill", blue).text(`Real plates: ${r.toFixed(2)}`);
  const parts = { band: [band], obs: [obs], inj: [injG, injHi, injLab] };
  Object.values(parts).flat().forEach(e => e.style("opacity", 0));
  const full0 = B0.map(d => d.c), full1 = B1.map(d => d.c), zero = B0.map(() => 0);
  const cl = Math.round(100 * (1 - I.p_as_high)), pv = Math.round(100 * I.p_as_high);
  const slides = [
    { band: 0, obs: 0, inj: 0, build: true, text: [["We run the fit many times with the real shadow outlines laid on random plates, which the shadow never crossed.", ink, 400],
      [`On average it finds no effect (a ratio of ${R.null.median.toFixed(2)}); the spread is from chance, due to fitting only ${Z.plates.length} plates.`, ink, 400]] },
    { band: 0, obs: 1, inj: 0, text: [[`The real plates give ${r.toFixed(2)}: no sign of fewer transients in shadow.`, blue, 600]] },
    { band: 0, obs: 1, inj: 1, shade: true, text: [[`Had transients been as rare in shadow as VASCO’s count implies (r = ${I.r_true.toFixed(2)}), a fit this high would occur ${pv}% of the time.`, pur, 400],
      [`A shortfall of that size is excluded at ${cl}% confidence (one-sided).`, pur, 600]] },
    { band: 0, obs: 1, inj: 1, shade: true, text: [[[[`Fit result on real plates (${r.toFixed(2)})`, blue], [" falls within random trials with no shadow effect.", ink]], ink, 400],
      [`A shortfall as large as VASCO’s (purple) is excluded at ${cl}% confidence (one-sided).`, pur, 600]] }
  ];
  const last = slides.length - 1, op = { band: 0, obs: 0, inj: 0 };
  let cur = -1, timers = [];
  const at = (ms, f) => timers.push(setTimeout(f, ms));
  const show = (k, dur) => {
    timers.forEach(clearTimeout); timers = []; s.selectAll("*").interrupt();
    const sl = slides[k], out = dur ? 250 : 0;
    for (const [key, els] of Object.entries(parts)) if (op[key] && !sl[key]) { els.forEach(e => e.transition().duration(out).style("opacity", 0)); op[key] = 0; }
    t1.style("opacity", 0); t2.style("opacity", 0);
    if (!sl.shade) clip.attr("width", 0);
    if (sl.build && dur) bars(nullG, B0, zero);
    at(out, () => {
      bars(injG, B1, full1); bars(injHi, B1, full1); injLab.attr("y", y(d3.max(full1)) - 8).text(`VASCO-sized shortfall (r = ${I.r_true.toFixed(2)})`);
      [t1, t2].forEach((el, i) => { const tx = sl.text[i];                      // a line is a string, or [text, colour] pieces
        el.text("").style("fill", tx ? tx[1] : ink).style("font-weight", tx ? tx[2] : 400);
        if (tx && Array.isArray(tx[0])) tx[0].forEach(([str, col]) => el.append("tspan").style("fill", col).text(str)); else el.text(tx ? tx[0] : "");
        el.transition().duration(dur ? 400 : 0).style("opacity", 1); });
      for (const [key, els] of Object.entries(parts)) if (sl[key] && !op[key]) { els.forEach(e => e.style("opacity", 0).transition().duration(dur ? 450 : 0).style("opacity", 1)); op[key] = 1; }
      if (sl.shade) clip.attr("width", dur && cur !== last ? 0 : W - x(r)).transition().delay(dur ? 450 : 0).duration(dur ? 1000 : 0).attr("width", W - x(r));
      if (!(sl.build && dur)) return bars(nullG, B0, full0);
      const cum = d3.cumsum(full0), tot = cum[cum.length - 1]; let seed = 11; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
      const h = zero.slice(); let done = 0; const N = Math.round(n0);
      for (let k2 = 1; k2 <= 40; k2++) at(200 + k2 * 60, () => { const target = Math.round(N * (k2 / 40) ** 2.2);
        for (; done < target; done++) h[d3.bisectRight(cum, rnd() * tot)] += 1; bars(nullG, B0, h); });
      at(200 + 41 * 60, () => bars(nullG, B0, full0));
    });
    cur = k;
  };
  const ctl = d3.select("#ratio-steps"), back = ctl.select(".back"), next = ctl.select(".next"), pos = ctl.select(".pos");
  const go = (k, instant) => {
    show(k, instant || matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1);
    back.property("disabled", k === 0);
    next.text(k === last ? "Step through ›" : k === last - 1 ? "Summary ›" : "Next ›");
    pos.text(k === last ? "Summary" : `${k + 1} of ${last}`);
  };
  back.on("click", () => go(cur - 1));
  next.on("click", () => go(cur === last ? 0 : cur + 1));
  go(last, true);
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
  s.append("text").attr("x", x1).attr("y", top - 22).style("font-weight", 600).style("font-size", "12.5px").text("Track-finding efficiency, ε");
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

/* ---------- the limits on glint tracks of 3, 4, 5, and 6 transients: in each panel, the probability of each number of
   them, given the data (every number equally likely beforehand), with the bars holding 95% of the probability shaded ---------- */
function drawLimit(D, id) {
  const L = D.limit, s = svg(id);
  if (!L) return;
  const W = 760, pw = 306, ph = 150, gx = 72, gy = 78, m = { l: 64, t: 34 }, keys = ["3", "4", "5", "6"];   // 64 + 2 x 306 + 72 = 748
  s.attr("viewBox", `0 0 ${W} ${m.t + 2 * ph + gy + 44}`);
  keys.forEach((f, i) => {
    const P = L[f], x0 = m.l + (i % 2) * (pw + gx), y0 = m.t + Math.floor(i / 2) * (ph + gy);
    const ns = d3.range(P.p.length), g = s.append("g").attr("transform", `translate(${x0},${y0})`);
    const x = d3.scaleBand().domain(ns).range([0, pw]).padding(0.15);
    const y = d3.scaleLinear().domain([0, Math.ceil(d3.max(P.p) * 10) / 10]).range([ph, 0]);
    g.selectAll("rect").data(ns).join("rect").attr("x", n => x(n)).attr("width", x.bandwidth()).attr("y", n => y(P.p[n])).attr("height", n => ph - y(P.p[n]))
      .attr("fill", n => n <= P.s95 ? V("--ours") : V("--ink-3")).attr("fill-opacity", n => n <= P.s95 ? 0.85 : 0.45)
      .call(sel => hov(sel, n => `Glint tracks of ${f} transients: ${n} in the catalogue, probability ${d3.format(".1%")(P.p[n])}`));
    const step = ns.length > 30 ? 10 : ns.length > 12 ? 5 : ns.length > 6 ? 2 : 1;
    g.append("g").attr("transform", `translate(0,${ph})`).call(d3.axisBottom(x).tickValues(ns.filter(n => n % step === 0)).tickSizeOuter(0)).call(axisStyle);
    g.append("g").call(d3.axisLeft(y).ticks(3).tickFormat(d3.format(".0%")).tickSize(3)).call(axisStyle);
    g.append("text").attr("x", 0).attr("y", -12).style("font-weight", 600).style("font-size", "12.5px").text(`${f} transients`);
    g.append("text").attr("x", pw).attr("y", -12).attr("text-anchor", "end").style("font-size", "12px").style("fill", V("--ink-2"))
      .text(`95% probability: at most ${P.s95}`);
    if (i >= 2) g.append("text").attr("x", pw / 2).attr("y", ph + 34).attr("text-anchor", "middle").text("Glint tracks in the catalogue");
    g.append("text").attr("transform", `translate(-46,${ph / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Probability");
  });
}

/* ---------- alignments: VASCO's search and score, real against shuffled copies ---------- */
/* ---------- a score's distribution by candidate size, real catalogue against shifted copies, with the difference in
   standard deviations: the track fit (step 3) and VASCO's own score (step 4) ---------- */
function drawTracks(D) {
  drawDist(D.trackdist, "c-counts", { dom: [-7, 0], step: 0.2, ticks: [-6, -4, -2, 0], fmt: { "-6": "10⁻⁶", "-4": "10⁻⁴", "-2": "0.01", "0": "1" },
    xlab: "Track score: expected number of such tracks by chance (lower is less likely)", what: "track score",
    rowlab: "Number of track candidates at each score", floor: 0.01, yticks: [0.01, 0.1, 1, 10, 100] });
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
      const bins = lo.map((l, k) => ({ l, h: l + step, n: d.real[k], b: d.bmean[k], sd: d.bsd[k], z: d.z ? d.z[k] : 0 })).filter(q => q.l >= o.dom[0] - 1e-9 && q.h <= o.dom[1] + 1e-9);
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
        // each bin's difference as a bar from zero: the real count's place among the copies, in standard deviations
        const y = d3.scaleLinear().domain([-4, 4]).range([bot, top]).clamp(true), bw = x(o.dom[0] + step) - x(o.dom[0]);
        s.append("g").selectAll("rect").data(bins).join("rect").attr("class", "hov")
          .attr("x", q => x(q.l) + 0.5).attr("width", Math.max(1, bw - 1)).attr("y", q => Math.min(y(0), y(q.z))).attr("height", q => Math.abs(y(q.z) - y(0)))
          .attr("fill", V("--emulsion")).attr("fill-opacity", 0.85)
          .call(sel => hov(sel, q => `${d.r}-point candidates, ${o.what} ${(10 ** q.l).toPrecision(2)} to ${(10 ** q.h).toPrecision(2)}<br>Real catalogue ${d3.format(",")(q.n)}; chance ${d3.format(",.2f")(q.b)} on average<br>${q.z >= 0 ? "+" : "−"}${Math.abs(q.z).toFixed(1)} standard deviations`));
        s.append("line").attr("x1", x0).attr("x2", x0 + pw).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", V("--ink-3"));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([-4, -2, 0, 2, 4]).tickFormat(v => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v) + "σ").tickSize(3)).call(axisStyle);
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
  const R = D.radial, s = svg("c-radial"), W = 760, H = 320, m = { l: 62, r: 16, t: 14, b: 42 };
  const x = d3.scaleLinear().domain([0, 3]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, 3]).range([H - m.b, m.t]);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(6).tickFormat(d => d + "°").tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickFormat(d => d + "×").tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Distance from the centre of the plate");
  s.append("text").attr("transform", `translate(16,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Density relative to the centre");
  const line = k => d3.line().x((v, i) => x(R.r[i])).y(v => y(v)).curve(d3.curveMonotoneX)(R[k]);
  s.append("path").attr("d", d3.area().x((v, i) => x(R.r[i])).y0((v, i) => y(R.lo[i])).y1((v, i) => y(R.hi[i])).curve(d3.curveMonotoneX)(R.r)).attr("fill", V("--ours-soft"));
  s.append("path").attr("d", line("control")).attr("fill", "none").attr("stroke", V("--ink-3")).attr("stroke-width", 2).attr("stroke-dasharray", "4 3");
  s.append("path").attr("d", line("transients")).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.4);
  s.append("g").selectAll("circle").data(R.transients).join("circle").attr("cx", (v, i) => x(R.r[i])).attr("cy", v => y(v)).attr("r", 3.5).attr("fill", V("--ours"))
    .call(sel => hov(sel, (v, i) => `${v.toFixed(1)} times the density near the centre`));
}
let crowdState = "graph";
function drawCrowd(D) {
  // every plate as on the sky map; the buttons morph each one between its place on the sky and a tile at its distance from
  // the Milky Way and its count, with a fit for the typical plate
  const s = svg("c-crowd"), W = 760, H = 400, m = { l: 56, r: 16, t: 14, b: 42 }, F = D.crowd_fit;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const poly = new Map(D.plates.map(p => [p.id, p.poly])), pts = D.crowding.filter(p => poly.has(p.id));
  const x = d3.scaleLinear().domain([0, 90]).range([m.l, W - m.r]);
  const y = d3.scaleLog().domain([3, 3000]).range([H - m.b, m.t]).clamp(true);
  const col = d3.scaleSequentialLog(d3.interpolateViridis).domain([5, 1500]).clamp(true);
  // the sky map's projection, scaled to this chart
  const clip = d3.max(D.plates.flatMap(p => p.poly), v => 90 - v[1]) + 2.5, Rd = (H - 20) / 2, cx = W / 2, cy = H / 2;
  const proj = d3.geoAzimuthalEqualArea().reflectX(true).rotate([0, -90]).clipAngle(clip).scale(Rd / (2 * Math.sin(clip * DEG / 2))).translate([cx, cy]);
  const K = 24, xyz = ([l, b]) => [Math.cos(b * DEG) * Math.cos(l * DEG), Math.cos(b * DEG) * Math.sin(l * DEG), Math.sin(b * DEG)];
  const lonlat = v => { const n = Math.hypot(...v); return [Math.atan2(v[1], v[0]) / DEG, Math.asin(v[2] / n) / DEG]; };
  const area = a => d3.sum(a, (p, i) => p[0] * a[(i + 1) % a.length][1] - a[(i + 1) % a.length][0] * p[1]);
  const start = a => { const i = d3.minIndex(a, p => p[0] + p[1]); return a.slice(i).concat(a.slice(0, i)); };   // from the top-left
  const norm = a => start(area(a) < 0 ? a.slice().reverse() : a);
  const resample = a => {                                                       // exactly K points, evenly spaced round the outline
    const L = [0]; a.forEach((p, i) => { const q = a[(i + 1) % a.length]; L.push(L[i] + Math.hypot(q[0] - p[0], q[1] - p[1])); });
    return d3.range(K).map(k => { const s_ = k / K * L[a.length]; let i = d3.bisectRight(L, s_) - 1; i = Math.min(i, a.length - 1);
      const p = a[i], q = a[(i + 1) % a.length], f = (s_ - L[i]) / ((L[i + 1] - L[i]) || 1); return [p[0] + (q[0] - p[0]) * f, p[1] + (q[1] - p[1]) * f]; });
  };
  const skyShape = id => {                                                     // the plate outline, traced finely, then K points
    const c = poly.get(id), out = [];
    c.forEach((v, i) => { const a = xyz(v), b = xyz(c[(i + 1) % c.length]);
      for (let j = 0; j < 12; j++) { const t = j / 12; out.push(proj(lonlat(a.map((q, k) => q + (b[k] - q) * t)))); } });
    return norm(resample(out));
  };
  const tile = 6.5, sqShape = p => {                                            // a small square at the plate's place on the graph
    const X = x(p.b), Y = y(Math.max(p.N, 3)), h = tile / 2, c = [[X - h, Y - h], [X + h, Y - h], [X + h, Y + h], [X - h, Y + h]], out = [];
    c.forEach((a, i) => { const b = c[(i + 1) % 4]; for (let j = 0; j < K / 4; j++) { const t = j / (K / 4); out.push([a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]); } });
    return norm(out);
  };
  const shapes = new Map(pts.map(p => [p.id, { sky: skyShape(p.id), graph: sqShape(p) }]));
  const d = a => "M" + a.map(q => q.map(v => v.toFixed(1)).join(",")).join("L") + "Z";
  // the graph's furniture
  const G = s.append("g");
  G.append("path").attr("d", d3.area().x((_, i) => x(F.b[i])).y0((_, i) => y(F.q25[i])).y1((_, i) => y(F.q75[i])).curve(d3.curveMonotoneX)(F.b)).attr("fill", V("--ours")).attr("fill-opacity", 0.16);
  G.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).ticks(9).tickFormat(v => v + "°").tickSizeOuter(0)).call(axisStyle);
  G.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([3, 10, 30, 100, 300, 1000, 3000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
  G.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Degrees from the plane of the Milky Way");
  G.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients per plate");
  // the sky's furniture
  const path = d3.geoPath(proj), S = s.append("g");
  S.append("path").datum({ type: "Sphere" }).attr("d", path).attr("fill", V("--paper")).attr("stroke", V("--rule"));
  const plates = s.append("g").selectAll("path").data(pts).join("path").attr("fill", p => col(Math.max(p.N, 5))).attr("stroke", V("--surface")).attr("stroke-width", 0.5)
    .call(sel => hov(sel, p => `${p.id}: ${fmt(p.N)} transients, ${p.b}° from the Milky Way`));
  const plane = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => galToEq(l, 0)) };
  const mw = s.append("g").attr("pointer-events", "none");                       // as on the sky map: a soft band with a dashed centre line
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-opacity", 0.38).attr("stroke-width", 10).attr("stroke-linecap", "round");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#111418").attr("stroke-opacity", 0.85).attr("stroke-width", 1.5).attr("stroke-dasharray", "5 4");
  S.lower(); G.lower();
  const fitG = s.append("g").attr("pointer-events", "none"), line = k => d3.line().x((_, i) => x(F.b[i])).y(v => y(v)).curve(d3.curveMonotoneX)(F[k]);
  fitG.append("path").attr("d", line("mid")).attr("fill", "none").attr("stroke", V("--surface")).attr("stroke-width", 6);
  fitG.append("path").attr("d", line("mid")).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 3);
  const layers = { sky: [S, mw], graph: [G, fitG] };
  const show = (key, t) => {
    const from = crowdState; crowdState = key;
    d3.selectAll("#c-crowd").each(function () { const L = d3.select(this.parentNode);
      for (const k of ["sky", "graph"]) L.selectAll(`.crowd-${k}`).style("display", k === key ? null : "none"); });
    // the furniture of the view being left goes first, then the plates move, then the new view's furniture arrives
    const lead = t && from !== "sky" ? 600 : 0, toSky = key === "sky";
    plates.transition().duration(t).delay(p => t ? lead + (toSky ? 90 - p.b : p.b) * 12 : 0).ease(d3.easeCubicInOut)
      .attrTween("d", p => { const i = d3.interpolate(shapes.get(p.id)[from], shapes.get(p.id)[key]); return t_ => d(i(t_)); })
      .attr("stroke-width", toSky ? 0.5 : 0.6);
    if (!t) plates.attr("d", p => d(shapes.get(p.id)[key]));
    const fade = (sel, on, wait) => sel.transition().duration(t ? 450 : 0).delay(t ? wait : 0).style("opacity", on ? 1 : 0);
    for (const [k, ls] of Object.entries(layers)) ls.forEach((l, j) => k === key ? fade(l, true, lead + t + 500 + 300 * j) : fade(l, false, 150 * j));
  };
  show(crowdState, 0);
  d3.selectAll("#crowd-toggle button").on("click", function () {
    if (this.dataset.key === crowdState) return;
    d3.selectAll("#crowd-toggle button").classed("on", false); d3.select(this).classed("on", true);
    show(this.dataset.key, matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1100);
  }).classed("on", function () { return this.dataset.key === crowdState; });
}

/* ---------- the pointing curve's calibration (plate_properties.pointing): every plate's actual count against its predicted
   count, and in eight groups by prediction, each group's average actual count against its average prediction ---------- */
function drawPointing(D) {
  const s = svg("c-pointing"), W = 760, H = 400, m = { l: 56, r: 16, t: 14, b: 42 }, PF = D.pointing_fit, pts = D.crowding;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const [e0, e1] = d3.extent(pts, p => p.E), x = d3.scaleLog().domain([e0 / 1.15, e1 * 1.15]).range([m.l, W - m.r]);   // the predictions' range
  const y = d3.scaleLog().domain([2, 5000]).range([H - m.b, m.t]).clamp(true);
  const col = d3.scaleSequentialLog(d3.interpolateViridis).domain([5, 1500]).clamp(true);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([3, 5, 10, 20, 30, 50, 100, 200, 300, 500].filter(v => v >= x.domain()[0] && v <= x.domain()[1]))
    .tickFormat(d3.format(",")).tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([3, 10, 30, 100, 300, 1000, 3000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Transients the plate’s pointing predicts");
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients on the plate");
  const t = 6;
  s.append("g").selectAll("rect").data(pts).join("rect").attr("x", p => x(p.E) - t / 2).attr("y", p => y(Math.max(p.N, 2)) - t / 2).attr("width", t).attr("height", t).attr("rx", 1)
    .attr("fill", p => col(Math.max(p.N, 5))).attr("fill-opacity", 0.45)
    .call(sel => hov(sel, p => `${p.id}: ${fmt(p.N)} transients, ${Math.round(p.E)} predicted from its pointing`));
  const diag = `M${x(x.domain()[0])},${y(x.domain()[0])}L${x(x.domain()[1])},${y(x.domain()[1])}`;
  s.append("path").attr("d", diag).attr("stroke", V("--surface")).attr("stroke-width", 6).attr("pointer-events", "none");
  s.append("path").attr("d", diag).attr("stroke", V("--ours")).attr("stroke-width", 2.5).attr("pointer-events", "none");
  const G = s.append("g").selectAll("g").data(PF.calibration).join("g");
  G.append("line").attr("x1", c => x(c.E)).attr("x2", c => x(c.E)).attr("y1", c => y(c.lo)).attr("y2", c => y(c.hi)).attr("stroke", V("--ink")).attr("stroke-width", 2);
  G.append("circle").attr("cx", c => x(c.E)).attr("cy", c => y(c.N)).attr("r", 5.5).attr("fill", V("--ink")).attr("stroke", V("--surface")).attr("stroke-width", 1.5)
    .call(sel => hov(sel, c => `${c.n} plates predicted ${Math.round(c.E)} on average<br>Actual average ${Math.round(c.N)} (${Math.round(c.lo)} to ${Math.round(c.hi)})`));
}

/* ---------- the equatorial row: each plate's count over what its pointing predicts without the row term, by the declination
   of its centre, with each row of plates' median ---------- */
function drawRows(D) {
  const R = D.rows, s = svg("c-rows"), W = 760, H = 300, m = { l: 56, r: 16, t: 14, b: 42 };
  if (!R) return;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const x = d3.scaleLinear().domain([-3, 90]).range([m.l, W - m.r]), y = d3.scaleLog().domain([0.01, 20]).range([H - m.b, m.t]).clamp(true);
  s.append("rect").attr("x", x(-3)).attr("width", x(3) - x(-3)).attr("y", m.t).attr("height", H - m.b - m.t).attr("fill", V("--theirs")).attr("fill-opacity", 0.08);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 15, 30, 45, 60, 75, 90]).tickFormat(v => v + "°").tickSizeOuter(0)).call(axisStyle);
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0.01, 0.1, 1, 10]).tickFormat(v => v + "×").tickSize(3)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 6).attr("text-anchor", "middle").text("Declination of the plate’s centre");
  s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients over prediction");
  s.append("line").attr("x1", x(-3)).attr("x2", x(90)).attr("y1", y(1)).attr("y2", y(1)).attr("stroke", V("--ink-3")).attr("stroke-dasharray", "3 3");
  s.append("g").selectAll("circle").data(R.plates).join("circle").attr("cx", p => x(p.dec)).attr("cy", p => y(p.r)).attr("r", 2.6)
    .attr("fill", V("--ink-3")).attr("fill-opacity", 0.45).call(sel => hov(sel, p => `${p.id}: ${p.r.toFixed(2)} times its prediction, centre at ${p.dec}°`));
  const eb = s.append("g").attr("stroke", V("--ours")).attr("stroke-width", 2);
  eb.selectAll("path").data(R.rows).join("path").attr("fill", "none")
    .attr("d", r => `M${x(r.dec)},${y(r.lo)}V${y(r.hi)}M${x(r.dec) - 5},${y(r.lo)}h10M${x(r.dec) - 5},${y(r.hi)}h10`);
  s.append("g").selectAll("circle.m").data(R.rows).join("circle").attr("cx", r => x(r.dec)).attr("cy", r => y(r.median)).attr("r", 3.2).attr("fill", V("--ours")).call(sel => hov(sel, r => `Row at ${r.dec}°: ${r.n} plates, median ${r.median.toFixed(2)} times their prediction (${r.lo.toFixed(2)} to ${r.hi.toFixed(2)})`));
  s.append("text").attr("x", x(3) + 6).attr("y", m.t + 14).style("font-size", "11px").style("fill", V("--theirs")).text("Equatorial row");
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
  const sd = sdText(q.sigma);
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 14).style("fill", V("--emulsion")).text(`${line}: ${q.real.toFixed(2)}×`);
  s.append("text").attr("class", "ann").attr("x", x(q.real) + 8).attr("y", top + 30).style("fill", V("--emulsion")).text(`${trials} reach it 1 time in ${Math.round(1 / q.p)} (${sd})`);
  s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).ticks(6).tickFormat(d => d + "×").tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text(xlab);
  s.append("text").attr("transform", `translate(12,${(top + bot) / 2}) rotate(-90)`).attr("text-anchor", "middle").text(o.ylab || "Share of random sets");
}

let DATA = null;
/* ---------- how the track score rates one real candidate: the candidate on the sky in its box; the straightness test it
   must pass first (with the close-pair rule as a schematic); and the four numbers that set the score ---------- */
function drawMethod(D, id) {
  // slides the reader steps through, in the order the method runs: 1 a real candidate; 2 the close-pair rule that forms
  // it; 3 the straightness test, which throws away candidates under the cut; 4 the four numbers measured on one that
  // passes; 5 the track score; and a summary, shown on load. Each slide declares how strongly each part shows
  const M = D.method, s = svg(id);
  if (!M) return;
  const W = 760, H = 398, top = 92, ours = V("--ours"), grey = V("--ink-3"), ink = V("--ink"), ink2 = V("--ink-2"), rule = V("--rule");
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const strip = Math.max(M.w, 1), N = M.N, r = M.r, one = d3.format(",")(Math.round(1 / M.score / 100) * 100), pc = Math.round(100 * M.p_chi2);
  const head = (g, x, t) => g.append("text").attr("x", x).attr("y", top - 14).style("font-weight", 600).style("font-size", "13px").style("fill", ink).text(t);
  const lab = (g, x, y, t, o = {}) => g.append("text").attr("x", x).attr("y", y).attr("text-anchor", o.anchor || "start")
    .style("font-size", o.size || "11px").style("font-weight", o.bold ? 600 : null).style("fill", o.col || ink2).text(t);
  const part = () => s.append("g").style("opacity", 0);
  const sky = part(), pair = part(), fit = part(), meas = part(), score = part();

  // the candidate on the sky, to scale
  const S = 196, half = 0.8 * M.side, b = M.side / 2, g1 = sky.append("g").attr("transform", `translate(0,${top})`), m1 = meas.append("g").attr("transform", `translate(0,${top})`);
  const sx = d3.scaleLinear().domain([-half, half]).range([0, S]), sy = d3.scaleLinear().domain([-half, half]).range([S, 0]);
  head(sky, 0, `The candidate (${M.plate})`);
  g1.append("rect").attr("width", S).attr("height", S).attr("rx", 8).attr("fill", V("--surface")).attr("stroke", rule);
  m1.append("rect").attr("x", sx(-b)).attr("y", sy(b)).attr("width", sx(b) - sx(-b)).attr("height", sy(-b) - sy(b)).attr("fill", "none").attr("stroke", ink2).attr("stroke-dasharray", "4 3");
  lab(m1, sx(-b), sy(b) - 6, `${Math.round(M.side)}′ box`);
  const m = M.members, ex = d3.extent(m, d => d[0]), e0 = m.find(d => d[0] === ex[0]), e1 = m.find(d => d[0] === ex[1]);
  const p0 = [sx(e0[0]), sy(e0[1])], p1 = [sx(e1[0]), sy(e1[1])], dl = Math.hypot(p1[0] - p0[0], p1[1] - p0[1]);
  const ux = (p1[0] - p0[0]) / dl, uy = (p1[1] - p0[1]) / dl, side = -uy > 0 ? -1 : 1;
  const at = (t, d) => [p0[0] + ux * t - uy * d * side, p0[1] + uy * t + ux * d * side];
  g1.append("line").attr("x1", at(-10, 0)[0]).attr("y1", at(-10, 0)[1]).attr("x2", at(dl + 10, 0)[0]).attr("y2", at(dl + 10, 0)[1]).attr("stroke", ours).attr("stroke-opacity", 0.6);
  g1.selectAll("circle.o").data(M.others.filter(o => !m.some(q => Math.hypot(q[0] - o[0], q[1] - o[1]) < 0.01))).join("circle")
    .attr("cx", d => sx(d[0])).attr("cy", d => sy(d[1])).attr("r", 3.2).attr("fill", grey).attr("fill-opacity", d => Math.max(Math.abs(d[0]), Math.abs(d[1])) <= b ? 0.95 : 0.4);
  g1.selectAll("circle.m").data(m).join("circle").attr("cx", d => sx(d[0])).attr("cy", d => sy(d[1])).attr("r", 4.5).attr("fill", ours);
  const q0 = at(0, 15), q1 = at(dl, 15), lp = at(dl / 2, 40), rp = at(dl * 0.2, -16);
  m1.append("path").attr("d", `M${q0}L${q1}M${at(0, 11)}L${at(0, 19)}M${at(dl, 11)}L${at(dl, 19)}`).attr("stroke", ink2).attr("fill", "none");
  lab(m1, lp[0], lp[1] + 4, `L = ${M.L.toFixed(1)}′`, { anchor: "middle", bold: true, col: ours });
  lab(m1, rp[0], rp[1] + 4, `r = ${r}`, { anchor: "middle", bold: true, col: ours });
  lab(m1, sx(0), sy(-b) + 16, `N = ${N} in the box`, { anchor: "middle", bold: true, col: ours });

  // the close-pair rule, as a schematic, under the candidate
  const y3 = top + S + 30, g3 = pair.append("g").attr("transform", `translate(0,${y3})`), xs3 = d3.scaleLinear().domain([0, 160]).range([12, S - 12]), yl = 22;
  lab(g3, 0, 0, "Close pairs count once (schematic)", { bold: true, col: ink, size: "12px" });
  g3.append("line").attr("x1", 0).attr("x2", S).attr("y1", yl).attr("y2", yl).attr("stroke", ours);
  [[15, 3, 1], [92, -3, 1], [104, 8, 0], [150, -3, 1]].forEach(([a, o, keep]) => g3.append("circle").attr("cx", xs3(a)).attr("cy", yl - o).attr("r", 4.5)
    .attr("fill", keep ? ours : V("--surface")).attr("stroke", ours).attr("stroke-width", keep ? 0 : 1.5));
  [`Two transients fall within ${M.occupied}″ of each`, "other along the line; the transient", "closer to the line is kept"].forEach((l, i) => lab(g3, 0, yl + 24 + 14 * i, l));

  // the straightness test: the track laid flat, distances across it enlarged, and the chi-squared probability against the cut
  const x2 = 246, w2 = 248, h2 = 96, g2 = fit.append("g").attr("transform", `translate(${x2},${top})`);
  head(fit, x2, "1. Is it straight?");
  g2.append("rect").attr("width", w2).attr("height", h2).attr("rx", 8).attr("fill", V("--surface")).attr("stroke", rule);
  const k = 20, yc = h2 / 2, xa = d3.scaleLinear().domain([0, M.L]).range([18, w2 - 18]);
  g2.append("rect").attr("x", 10).attr("width", w2 - 20).attr("y", yc - k).attr("height", 2 * k).attr("fill", ours).attr("fill-opacity", 0.13);
  g2.append("line").attr("x1", 10).attr("x2", w2 - 10).attr("y1", yc).attr("y2", yc).attr("stroke", ours);
  M.res.forEach(([a, o]) => {
    g2.append("line").attr("x1", xa(a)).attr("x2", xa(a)).attr("y1", yc).attr("y2", yc - k * o).attr("stroke", ours);
    g2.append("circle").attr("cx", xa(a)).attr("cy", yc - k * o).attr("r", 4.5).attr("fill", ours);
  });
  lab(g2, w2 - 12, yc + k + 13, "±1″", { anchor: "end" });
  [[at(-6, 0), [x2, top + 8]], [at(dl + 6, 0), [x2, top + h2 - 8]]].forEach(([pa, pb]) =>
    fit.append("line").attr("x1", pa[0]).attr("y1", top + pa[1]).attr("x2", pb[0]).attr("y2", pb[1]).attr("stroke", ours).attr("stroke-opacity", 0.3).attr("stroke-dasharray", "2 3"));
  // the chi-squared probability on a log scale: below the cut a candidate is thrown away
  const yg = h2 + 34, hg = 16, xg = d3.scaleLog().domain([0.001, 1]).range([0, w2]);
  lab(g2, 0, yg - 8, "χ² probability of a straight line");
  g2.append("rect").attr("x", 0).attr("width", xg(M.p_min)).attr("y", yg).attr("height", hg).attr("fill", grey).attr("fill-opacity", 0.35);
  g2.append("rect").attr("x", xg(M.p_min)).attr("width", w2 - xg(M.p_min)).attr("y", yg).attr("height", hg).attr("fill", ours).attr("fill-opacity", 0.18);
  lab(g2, xg(M.p_min) / 2, yg + 11.5, "thrown away", { anchor: "middle", size: "10.5px" });
  lab(g2, (xg(M.p_min) + xg(M.p_chi2)) / 2, yg + 11.5, "kept", { anchor: "middle", size: "10.5px", col: ours });
  g2.append("line").attr("x1", xg(M.p_chi2)).attr("x2", xg(M.p_chi2)).attr("y1", yg - 3).attr("y2", yg + hg + 3).attr("stroke", ours).attr("stroke-width", 2.5);
  lab(g2, xg(M.p_chi2), yg - 8, `${pc}%`, { anchor: "middle", bold: true, col: ours });
  [0.001, 0.01, 0.1, 1].forEach(v => lab(g2, xg(v), yg + hg + 13, `${+(100 * v).toPrecision(1)}%`, { anchor: v === 0.001 ? "start" : v === 1 ? "end" : "middle", size: "10px" }));

  // the score: four numbers in, one out
  const x4 = 538, g4 = meas.append("g").attr("transform", `translate(${x4},${top})`), g5 = score.append("g").attr("transform", `translate(${x4},${top})`);
  head(meas, x4, "2. How often would chance do this?");
  [[`r = ${r}`, "transients"], [`L = ${M.L.toFixed(1)}′`, "long"], [`w = ±${strip.toFixed(0)}″`, "wide"], [`N = ${N}`, "in the box"]].forEach(([v, t], i) => {
    const cx = (i % 2) * 112, cy = Math.floor(i / 2) * 50;
    g4.append("rect").attr("x", cx).attr("y", cy).attr("width", 104).attr("height", 40).attr("rx", 8).attr("fill", ours).attr("fill-opacity", 0.1).attr("stroke", ours).attr("stroke-opacity", 0.4);
    lab(g4, cx + 52, cy + 18, v, { anchor: "middle", bold: true, col: ours, size: "12.5px" });
    lab(g4, cx + 52, cy + 32, t, { anchor: "middle", size: "10.5px" });
  });
  g5.append("path").attr("d", "M108,108 l0,18 m-6,-7 l6,7 l6,-7").attr("stroke", ink2).attr("fill", "none").attr("stroke-width", 1.5);
  g5.append("text").attr("x", 108).attr("y", 162).attr("text-anchor", "middle").style("font-size", "26px").style("font-weight", 600).style("fill", ours).text(`1 in ${one}`);
  lab(g5, 108, 182, "expected by chance in the box", { anchor: "middle" });

  const t1 = s.append("text").attr("x", 0).attr("y", 20).style("font-size", "13.5px"), t2 = s.append("text").attr("x", 0).attr("y", 42).style("font-size", "13.5px");
  const parts = { sky, pair, fit, meas, score }, dim = 0.3;
  const slides = [
    { sky: 1, text: [`A track candidate on plate ${M.plate}: two end points up to ${M.lmax}′ apart, and transients between them`,
      `lying within ${M.hmax}″ of the line joining them.`] },
    { sky: 1, pair: 1, text: ["Transients close together along the line count once, so a tight group cannot make a line look stronger."] },
    { sky: dim, pair: dim, fit: 1, text: [`We fit a straight line, allowing ${strip.toFixed(0)}″ per transient, the plates’ precision.`,
      `Candidates with a χ² probability under ${Math.round(100 * M.p_min)}% are thrown away. This one has ${pc}% and is kept.`] },
    { sky: 1, pair: dim, fit: dim, meas: 1, text: ["A kept candidate is measured: the number of transients in it (r), its length (L), its width (w, at least 1″),",
      "and the number of transients in a box around it (N)."] },
    { sky: dim, pair: dim, fit: dim, meas: 1, score: 1, text: ["Together they give the track score: the number of tracks like this expected by chance in the box.",
      `For this candidate, about 1 in ${one}.`] },
    { sky: 1, pair: 1, fit: 1, meas: 1, score: 1, text: ["Each candidate is first tested for straightness. Those that pass get a track score",
      "from four measured numbers."] }
  ];
  const last = slides.length - 1;
  let cur = -1, timers = [];
  const at2 = (ms, f) => timers.push(setTimeout(f, ms));
  const show = (k, dur) => {
    timers.forEach(clearTimeout); timers = []; s.selectAll("*").interrupt();
    const sl = slides[k], out = dur ? 250 : 0, lv = key => sl[key] || 0;
    for (const [key, g] of Object.entries(parts)) if (lv(key) < +g.style("opacity")) g.transition().duration(out).style("opacity", lv(key));
    [t1, t2].forEach(el => el.transition().duration(out).style("opacity", 0));
    at2(out, () => {
      [t1, t2].forEach((el, i) => el.text(sl.text[i] || "").style("fill", ink).transition().duration(dur ? 400 : 0).style("opacity", 1));
      for (const [key, g] of Object.entries(parts)) if (lv(key) > +g.style("opacity")) g.transition().duration(dur ? 450 : 0).style("opacity", lv(key));
    });
    cur = k;
  };
  const ctl = d3.select("#method-steps"), back = ctl.select(".back"), next = ctl.select(".next"), pos = ctl.select(".pos");
  const go = (k, instant) => {
    show(k, instant || matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1);
    back.property("disabled", k === 0);
    next.text(k === last ? "Step through ›" : k === last - 1 ? "Summary ›" : "Next ›");
    pos.text(k === last ? "Summary" : `${k + 1} of ${last}`);
  };
  back.on("click", () => go(cur - 1));
  next.on("click", () => go(cur === last ? 0 : cur + 1));
  go(last, true);
}

function render() {
  if (!DATA) return;
  drawSky(DATA, "c-sky0", false); drawSpread(DATA); drawSky(DATA, "c-sky", true); drawChance(DATA); drawPerm(DATA, "c-perm-any", "raw", "VASCO’s test");
  drawPerm(DATA, "c-nuc-corr", "corrected", "Per photograph, allowing for pointing", { line: "Plates near tests", trials: "Shifted calendars",
    xlab: "Transients per photograph near the chosen dates, against pointing, relative to the rest", ylab: "Share of shifted calendars" });
  drawFan(DATA); drawMethod(DATA, "c-method"); drawTracks(DATA); drawCompare(DATA, "c-compare"); drawLimit(DATA, "c-limit"); fillDates(DATA); drawRadialMap(DATA); drawRadial(DATA); drawCrowd(DATA); drawPointing(DATA); drawRows(DATA);
  drawPlates(DATA, "c-plates"); drawForest(DATA, "c-ratio-vasco"); drawRatioNull(DATA, "c-ratio-null");
}
fetch("data/site.json").then(r => r.json()).then(d => { DATA = d; render(); })
  .catch(() => document.querySelectorAll(".chart svg").forEach(s => s.insertAdjacentHTML("afterend", '<div class="loading">The chart data could not be loaded.</div>')));
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
new MutationObserver(render).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
