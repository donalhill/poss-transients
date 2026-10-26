// The page's charts, each drawn from figures.json (report/figures.py) and nothing else. Every
// <svg class="chart-svg" data-fig="..."> of the page is drawn by the function of that id,
// which also sets the chart's height.
const tip = document.getElementById("tip");
const V = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt = d3.format(",");
function hov(sel, html) {
  sel.on("mousemove", (ev, d) => {
    tip.innerHTML = html(d); tip.hidden = false;
    tip.style.left = Math.min(ev.clientX + 14, innerWidth - tip.offsetWidth - 10) + "px";
    tip.style.top = Math.min(ev.clientY + 14, innerHeight - tip.offsetHeight - 10) + "px";
  }).on("mouseleave", () => { tip.hidden = true; });
}
function axisStyle(g) { g.attr("class", "ax"); g.select(".domain").remove(); }
// a ratio axis on a log scale: ticks at powers of 2 within [lo, hi], so every gap is the same factor, with fractions below 1
const pow2Ticks = (lo, hi) => d3.range(Math.ceil(Math.log2(lo) - 1e-9), Math.floor(Math.log2(hi) + 1e-9) + 1).map(k => 2 ** k);
const pow2Label = v => v >= 1 ? d3.format("~g")(v) : `1/${Math.round(1 / v)}`;
const ann = (s, x, y, t, o = {}) => s.append("text").attr("class", "ann").attr("x", x).attr("y", y).attr("text-anchor", o.anchor || "start")
  .style("font-weight", o.bold ? 600 : null).style("fill", o.fill || null).text(t);
const size = (s, h) => s.attr("viewBox", `0 0 760 ${h}`);
const pctText = v => v < 0.01 && v > 0 ? (100 * v).toFixed(1) : v > 0.995 && v < 1 ? ">99" : String(Math.round(100 * v));
const DEG = Math.PI / 180, SHADOW = "#35d0f0";   // cyan: magma, the count map's scale, has no blues
function galToEq(l, b) {
  const L = l * DEG, B = b * DEG, aG = 192.85948 * DEG, dG = 27.12825 * DEG, lN = 122.93192 * DEG;
  const dec = Math.asin(Math.sin(dG) * Math.sin(B) + Math.cos(dG) * Math.cos(B) * Math.cos(lN - L));
  const y = Math.cos(B) * Math.sin(lN - L), x = Math.cos(dG) * Math.sin(B) - Math.sin(dG) * Math.cos(B) * Math.cos(lN - L);
  return [((aG + Math.atan2(y, x)) / DEG + 360) % 360, dec / DEG];
}
function ring(pts) {   // close a ring and wind it so d3 fills the small side
  const c = pts.concat([pts[0]]); let g = { type: "Polygon", coordinates: [c] };
  if (d3.geoArea(g) > 2 * Math.PI) g = { type: "Polygon", coordinates: [c.slice().reverse()] };
  return g;
}

// the sky: VASCO's plates, coloured by their transients, and, for the shadow map,
// the plates the shadow partly covered with its edge at mid-exposure, played in date order
// a plate's outline shrunk about its centre to the central 6 by 6 degrees (the plates are 6.5 degrees across on a
// 6-degree grid, so full outlines overlap and a missing plate's gap shows its neighbours' edges): gnomonic about the centre
function central(poly, side = 6.0) {
  const R = Math.PI / 180, v = poly.map(([a, d]) => [Math.cos(d * R) * Math.cos(a * R), Math.cos(d * R) * Math.sin(a * R), Math.sin(d * R)]);
  const m = [0, 1, 2].map(k => d3.mean(v, u => u[k])), a0 = Math.atan2(m[1], m[0]), d0 = Math.atan2(m[2], Math.hypot(m[0], m[1]));
  const fwd = ([a, d]) => { a *= R; d *= R; const c = Math.sin(d0) * Math.sin(d) + Math.cos(d0) * Math.cos(d) * Math.cos(a - a0);
    return [Math.cos(d) * Math.sin(a - a0) / c, (Math.cos(d0) * Math.sin(d) - Math.sin(d0) * Math.cos(d) * Math.cos(a - a0)) / c]; };
  const P = poly.map(fwd), half = d3.max(P, q => Math.max(Math.abs(q[0]), Math.abs(q[1]))), k = Math.tan(side / 2 * R) / half;
  return P.map(([x, y]) => { x *= k; y *= k; const r = Math.hypot(x, y), c = Math.atan(r);
    if (r === 0) return [a0 / R, d0 / R];
    const d = Math.asin(Math.cos(c) * Math.sin(d0) + y * Math.sin(c) * Math.cos(d0) / r), a = a0 + Math.atan2(x * Math.sin(c), r * Math.cos(d0) * Math.cos(c) - y * Math.sin(d0) * Math.sin(c));
    return [((a / R) % 360 + 360) % 360, d / R]; });
}
function skyMap(s, plates0, shadow, fid, o = {}) {
  const W = 760, H = 760, cx = W / 2, cy = H / 2, Rd = W / 2 - 44; size(s, H + (shadow ? 70 : 50));   // the colour bar (and the shadow's key) in a strip below the disc
  const plates = plates0.map(p => Object.assign({}, p, { g: ring(central(p.poly)) }));   // each plate's central 6 degrees: the grid's cells, so neighbours do not overlap
  const clip = d3.max(plates.flatMap(p => p.poly), v => 90 - v[1]) + 2.5;
  const proj = d3.geoAzimuthalEqualArea().reflectX(true).rotate([0, -90]).clipAngle(clip).scale(Rd / (2 * Math.sin(clip * DEG / 2))).translate([cx, cy]).precision(0.3);
  const path = d3.geoPath(proj), col = o.col || d3.scaleSequentialLog(t => d3.interpolateMagma(0.15 + 0.85 * t)).domain([10, Math.ceil(d3.max(plates0, p => p.N) / 100) * 100]).clamp(true);   // counts: log, 10 or fewer (zero too) at the floor
  const val = o.value || (p => p.N), has = p => o.value ? isFinite(val(p)) : true, fill = p => has(p) ? col(o.col ? val(p) : Math.max(val(p), 10)) : V("--surface");   // zero transients: the scale's floor
  const plane = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => galToEq(l, 0)) };
  s.append("path").datum({ type: "Sphere" }).attr("d", path).attr("fill", V("--paper")).attr("stroke", V("--rule"));
  s.append("path").datum(d3.geoGraticule().step([30, 30]).extentMajor([[-180, -80], [180, 80.01]])()).attr("d", path).attr("fill", "none").attr("stroke", V("--rule")).attr("stroke-width", 0.7);
  s.append("g").selectAll("path").data(plates.filter(p => !(shadow && p.circle))).join("path").attr("d", p => path(p.g))
    .attr("fill", fill).attr("fill-opacity", shadow ? 0.3 : 1).attr("stroke", V("--surface")).attr("stroke-width", 0.5)
    .call(sel => hov(sel, o.tip || (p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients`)));
  const mw = s.append("g").attr("pointer-events", "none");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-opacity", 0.38).attr("stroke-width", 16).attr("stroke-linecap", "round");
  mw.append("path").datum(plane).attr("d", path).attr("fill", "none").attr("stroke", "#111418").attr("stroke-opacity", 0.85).attr("stroke-width", 1.5).attr("stroke-dasharray", "5 4");
  const sh = shadow ? plates.filter(p => p.circle) : [];
  if (shadow) {
    const ecl = { type: "LineString", coordinates: d3.range(0, 361, 2).map(l => { const e = 23.44 * DEG, L = l * DEG;
      return [((Math.atan2(Math.sin(L) * Math.cos(e), Math.cos(L)) / DEG) + 360) % 360, Math.asin(Math.sin(e) * Math.sin(L)) / DEG]; }) };
    s.append("path").datum(ecl).attr("d", path).attr("fill", "none").attr("stroke", "#111418").attr("stroke-opacity", 0.5).attr("stroke-width", 4).attr("pointer-events", "none");   // a dark under-line, then the cyan
    s.append("path").datum(ecl).attr("d", path).attr("fill", "none").attr("stroke", SHADOW).attr("stroke-width", 2).attr("pointer-events", "none");
  }
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "shp").attr("pointer-events", "none").attr("d", p => path(p.g)).attr("fill", "none").attr("stroke", "#ffffff").attr("stroke-width", 5);
  s.append("g").selectAll("path").data(sh).join("path").attr("class", "shp").attr("d", p => path(p.g))
    .attr("fill", p => col(Math.max(p.N, 10))).attr("stroke", "#111418").attr("stroke-width", 1.8)
    .call(sel => hov(sel, p => `<b>Plate ${p.id}</b> · ${p.date}<br>${fmt(p.N)} transients in the area the analysis keeps<br>${Math.round(p.frac * 100)}% of the plate in shadow over the exposure<br>counted in shadow ${p.obs}; expected with no effect ${p.exp}`));
  s.append("g").attr("pointer-events", "none").selectAll("path").data(sh).join("path").attr("class", "shp").attr("d", p => path(ring(p.circle))).attr("fill", SHADOW).attr("fill-opacity", 0.14).attr("stroke", SHADOW).attr("stroke-width", 2.4);
  if (shadow) {
    const order = sh.slice().sort((a, b) => a.date.localeCompare(b.date)), all = s.selectAll(".shp"), dayf = d3.timeFormat("%b %Y");
    const lab = s.append("text").attr("x", W - 4).attr("y", 18).attr("text-anchor", "end").style("font-size", "15px").style("font-weight", 600).style("fill", V("--ink"));
    const sub = s.append("text").attr("x", W - 4).attr("y", 38).attr("text-anchor", "end").style("font-size", "14px").style("fill", V("--ink-2"));
    d3.select("#shadow-play").on("click", () => {
      const step = matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1300;
      all.interrupt().style("opacity", 0); lab.text(""); sub.text("");
      order.forEach((p, i) => setTimeout(() => {
        all.filter(q => q === p).transition().duration(step ? 400 : 0).style("opacity", 1);
        all.filter(q => q !== p && order.indexOf(q) < i).transition().duration(step ? 400 : 0).style("opacity", 0.35);
        lab.text(`${p.id}, ${dayf(new Date(p.date))}`); sub.text(`${p.obs} in shadow, ${p.exp} expected`);
      }, i * step));
      setTimeout(() => { all.transition().duration(step ? 400 : 0).style("opacity", 1); lab.text(""); sub.text(""); }, order.length * step + (step ? 600 : 0));
    });
  }
  const darkHalo = sel => sel.attr("class", "ann").style("font-weight", 600).style("fill", "#111418").style("paint-order", "stroke").style("stroke", "#ffffff").style("stroke-width", "3px");
  for (let h = 0; h < 24; h += 2) { const [x, y] = proj([h * 15, 90 - clip + 0.5]), r = Math.hypot(x - cx, y - cy), R = Rd + 16;
    s.append("text").attr("x", cx + (x - cx) / r * R).attr("y", cy + (y - cy) / r * R + 4).attr("text-anchor", "middle").attr("fill", V("--ink-3")).text(h + "h"); }
  for (const d of [0, 30, 60]) { const [x, y] = proj([172, d]);
    s.append("text").attr("x", x).attr("y", y - 3).attr("text-anchor", "middle").attr("fill", V("--ink-3")).style("paint-order", "stroke").style("stroke", V("--paper")).style("stroke-width", "3px").text((d > 0 ? "+" : "") + d + "°"); }
  const gid = `${fid}-grad`.replace(/\./g, "-"), bx = 8, by = H + 2, bw = 260, bh = 14, lx = (o.col ? d3.scaleLinear() : d3.scaleLog()).domain(col.domain()).range([bx, bx + bw]);
  s.append("defs").append("linearGradient").attr("id", gid).selectAll("stop").data(d3.range(0, 1.001, 0.1)).join("stop").attr("offset", t => t).attr("stop-color", t => col(lx.invert(bx + t * bw)));
  s.append("text").attr("x", bx).attr("y", by - 10).attr("fill", V("--ink-2")).style("font-size", "17px").text(o.legend || "Transients per plate");
  s.append("rect").attr("x", bx).attr("y", by).attr("width", bw).attr("height", bh).attr("rx", 2).attr("fill", `url(#${gid})`).attr("fill-opacity", shadow ? 0.3 : 1);
  (o.ticks || [10, 100, 1000]).forEach((v, i) => { s.append("line").attr("x1", lx(v)).attr("x2", lx(v)).attr("y1", by + bh).attr("y2", by + bh + 5).attr("stroke", V("--ink-3"));
    s.append("text").attr("x", lx(v)).attr("y", by + bh + 22).attr("text-anchor", i === 0 ? "start" : i === (o.ticks || [10, 100, 1000]).length - 1 ? "end" : "middle").attr("fill", V("--ink-2")).style("font-size", "15px").text(o.tickFormat ? o.tickFormat(v) : (i === 0 ? "0–10" : fmt(v))); });
  if (shadow) {                                                             // what the outline and the cyan mean, beside the colour bar
    const kx = W - 440, ky = H + 2;
    s.append("rect").attr("x", kx).attr("y", ky - 2).attr("width", 18).attr("height", 14).attr("fill", "none").attr("stroke", "#111418").attr("stroke-width", 1.8);
    s.append("text").attr("x", kx + 26).attr("y", ky + 10).attr("fill", V("--ink-2")).style("font-size", "15px").text("Plates in the measurement");
    s.append("ellipse").attr("cx", kx + 9).attr("cy", ky + 31).attr("rx", 9).attr("ry", 6).attr("fill", SHADOW).attr("fill-opacity", 0.14).attr("stroke", SHADOW).attr("stroke-width", 2.4);
    s.append("text").attr("x", kx + 26).attr("y", ky + 36).attr("fill", V("--ink-2")).style("font-size", "15px").text("The shadow at mid-exposure");
    s.append("line").attr("x1", kx).attr("x2", kx + 18).attr("y1", ky + 57).attr("y2", ky + 57).attr("stroke", SHADOW).attr("stroke-width", 2.5);
    s.append("text").attr("x", kx + 26).attr("y", ky + 62).attr("fill", V("--ink-2")).style("font-size", "15px").text("The shadow's path through the year, from the Earth's centre");
  }
  const m = proj(galToEq(100, 0)); s.append("text").attr("x", m[0] + 12).attr("y", m[1] - 10).call(darkHalo).text("Milky Way");
}
const mid = e => e.slice(0, -1).map((v, i) => (v + e[i + 1]) / 2);
// a legend row inside the chart: [label, colour token, "dot" | "line" | "band"]
function legend(s, x, y, items) {
  let cx = x;
  items.forEach(([lab, col, kind]) => {
    if (kind === "band") s.append("rect").attr("x", cx).attr("y", y - 8).attr("width", 16).attr("height", 10).attr("fill", V(col));
    else if (kind === "line") s.append("line").attr("x1", cx).attr("x2", cx + 16).attr("y1", y - 3).attr("y2", y - 3).attr("stroke", V(col)).attr("stroke-width", 2.5);
    else s.append("circle").attr("cx", cx + 6).attr("cy", y - 3).attr("r", 4.5).attr("fill", V(col));
    const t = ann(s, cx + 22, y + 1, lab); cx += 22 + t.node().getComputedTextLength() + 18;
  });
}



const shareText = p => {                                                        // a chance in words, as "1 time in 9"
  if (p < 0.25) return `1 time in ${Math.round(1 / p)}`;
  let best = [1, 2];
  for (let n = 2; n <= 5; n++) for (let k = 1; k < n; k++) if (Math.abs(k / n - p) < Math.abs(best[0] / best[1] - p)) best = [k, n];
  return `${best[0]} ${best[0] === 1 ? "time" : "times"} in ${best[1]}`;
};
const sdText = v => v.toFixed(1) === "1.0" ? "1 standard deviation" : `${v.toFixed(1)} standard deviations`;
const fid_ = fid => fid.replace(/\./g, "-");
// the count in shadow expected with no effect, given the plate's total: its count times the fit's no-effect share in shadow
const E0 = p => p.nw * p.weighted;
const shown_ = Z => Z.plates;                                                   // every plate in the measurement, however little it holds in shadow
const listed = q => q.length < 3 ? q.join(" and ") : `${q.slice(0, -1).join(", ")}, and ${q[q.length - 1]}`;
const notShown = Z => { const L = Z.left || [], T = Z.thin || [];                // the plates left out, one line per reason
  return [L.length ? `Left out: ${listed(L)}, with no part sunlit throughout the exposure.` : "", T.length ? `Left out: ${listed(T)}, with less than 1% of their area in shadow.` : ""].filter(Boolean); };

// each plate the shadow partly covers, as a day-night map: the share of its exposure each place spent in shadow, the edge at
// the start (dashed) and the end (solid), and every transient; the plate that weighs most in the fit drawn large
function platesChart(s, Z, fid) {
  const half = Z.half, day = V("--plate-day"), night = V("--plate-night"), edgeC = SHADOW, ink = V("--ink"), ink2 = V("--ink-2"), ink3 = V("--ink-3");
  const shade = d3.interpolateLab(day, night), narrow = s.node().getBoundingClientRect().width < 560, id = fid_(fid);
  const shown = shown_(Z), weight = p => p.nw * p.shadow * (1 - p.shadow);
  const all = shown.slice().sort((a, b) => weight(b) - weight(a));
  const played = [];                                                         // each panel's animation, for the play button
  // every plate the same size, five to a row (three when there are nine or fewer; two on a phone), the most informative first
  const W = narrow ? 400 : 760, cols = narrow ? 2 : all.length <= 9 ? 3 : 5, gap = narrow ? 16 : 20, lab = 52, top = 10, x0 = narrow ? 10 : 0;
  const SS = (W - 2 * x0 - (cols - 1) * gap) / cols, step = SS + lab + 6;
  const keyY = top + Math.ceil(all.length / cols) * step + 34;
  s.attr("viewBox", `0 0 ${W} ${keyY + (narrow ? 110 : 64)}`);
  function panel(p, x0, y0, S, big) {
    const k = S / (2 * half), px = v => (half - v) * k, py = v => (half - v) * k;   // east to the left, north up
    const g = s.append("g").attr("transform", `translate(${x0},${y0})`), cid = `${id}-clip-${p.plate}`;
    g.append("clipPath").attr("id", cid).append("rect").attr("width", S).attr("height", S);
    const inner = g.append("g").attr("clip-path", `url(#${cid})`);
    const nc = Math.round(Math.sqrt(p.cells.length)), cv = document.createElement("canvas"); cv.width = cv.height = nc;
    const ctx = cv.getContext("2d");
    for (let r = 0; r < nc; r++) for (let c = 0; c < nc; c++) {
      ctx.fillStyle = shade(parseInt(p.cells[(nc - 1 - r) * nc + (nc - 1 - c)], 36) / 35); ctx.fillRect(c, r, 1, 1);
    }
    const still = inner.append("g");                                         // the exposure as a whole: its share in shadow, both edges
    still.append("image").attr("href", cv.toDataURL()).attr("width", S).attr("height", S).attr("preserveAspectRatio", "none");
    const line = d3.line().x(v => px(v[0])).y(v => py(v[1]));
    p.edges.forEach((e, i) => {                                             // a dark outline under each edge, so the cyan shows on grey and black alike
      const dash = i === 0 ? (big ? "6 4" : "4 3") : null;
      still.append("path").attr("d", line(e)).attr("fill", "none").attr("stroke", "#0d1013").attr("stroke-opacity", 0.75).attr("stroke-width", big ? 5 : 3.8).attr("stroke-dasharray", dash);
      still.append("path").attr("d", line(e)).attr("fill", "none").attr("stroke", edgeC).attr("stroke-width", big ? 2.6 : 2).attr("stroke-dasharray", dash);
    });
    // the exposure played: the shadow (inside the edge's ring) and its edge, moving from the start of the exposure to its end
    const moving = inner.append("g").style("display", "none");
    moving.append("rect").attr("width", S).attr("height", S).attr("fill", day);
    const ring = f => line(p.sweep[f]) + "Z";
    const fill = moving.append("path").attr("fill", night).attr("d", ring(0));
    const under = moving.append("path").attr("fill", "none").attr("stroke", "#0d1013").attr("stroke-opacity", 0.75).attr("stroke-width", big ? 5 : 3.8).attr("d", ring(0));
    const edge = moving.append("path").attr("fill", "none").attr("stroke", edgeC).attr("stroke-width", big ? 2.6 : 2).attr("d", ring(0));
    played.push({ still, moving, frames: p.sweep.length, set: x => {      // x from 0 (start) to 1 (end), between the stored frames
      const a = Math.min(Math.floor(x * (p.sweep.length - 1)), p.sweep.length - 2), f = x * (p.sweep.length - 1) - a;
      const pts = p.sweep[a].map((q, i) => [q[0] + (p.sweep[a + 1][i][0] - q[0]) * f, q[1] + (p.sweep[a + 1][i][1] - q[1]) * f]);
      const d = line(pts) + "Z"; fill.attr("d", d); under.attr("d", d); edge.attr("d", d);
    } });
    inner.append("g").selectAll("circle").data(p.xi.map((v, i) => [px(v), py(p.eta[i])])).join("circle")
      .attr("cx", d => d[0]).attr("cy", d => d[1]).attr("r", big ? 1.9 : 1.15)
      .attr("fill", "#ffffff").attr("stroke", "#0d1013").attr("stroke-width", big ? 0.6 : 0.4).attr("stroke-opacity", 0.8);
    g.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", ink3).attr("stroke-width", 0.8);
    g.append("text").attr("y", S + (big ? 20 : 17)).style("font-weight", 600).style("font-size", big ? "14px" : "12.5px").style("fill", ink)
      .text(big ? `${p.plate} · ${pctText(p.shadow)}% in shadow` : `${p.plate} · ${pctText(p.shadow)}%`);
    (big ? [`${p.Ow.toFixed(1)} in shadow · ${E0(p).toFixed(1)} expected`] : [`${p.Ow.toFixed(1)} in shadow`, `${E0(p).toFixed(1)} expected`])
      .forEach((t, i) => g.append("text").attr("y", S + (big ? 38 : 32 + 14 * i)).style("font-size", big ? "12.5px" : "11px").style("fill", ink2).text(t));
    g.append("rect").attr("width", S).attr("height", S).attr("fill", "transparent")
      .call(sel => hov(sel, () => `<b>Plate ${p.plate}</b><br>${fmt(p.n)} transients, ${pctText(p.shadow)}% of the plate in shadow over the exposure<br>` +
        `counted by each transient's share of the exposure in shadow: ${p.Ow.toFixed(1)} in shadow, ${(p.nw - p.Ow).toFixed(1)} in sunlight<br>` +
        `expected in shadow with no effect: ${E0(p).toFixed(1)}<br>ratio ${p.r.toFixed(2)} (${p.lo.toFixed(2)} to ${p.hi == null ? "–" : p.hi.toFixed(2)}, counting only)`));
    return { g, px, py };
  }
  all.forEach((p, i) => panel(p, x0 + (i % cols) * (SS + gap), top + Math.floor(i / cols) * step, SS, false));
  const kw = 220, kx = x0, key = s.append("g").attr("transform", `translate(${kx},${keyY})`), gid = `${id}-ramp`;
  const grad = s.append("defs").append("linearGradient").attr("id", gid);
  d3.range(0, 1.001, 0.1).forEach(t => grad.append("stop").attr("offset", t).attr("stop-color", shade(t)));
  key.append("rect").attr("width", kw).attr("height", 8).attr("fill", `url(#${gid})`).attr("stroke", ink3).attr("stroke-width", 0.5);
  key.append("text").attr("y", 22).style("font-size", "11px").style("fill", ink2).text("Sunlit throughout");
  key.append("text").attr("x", kw).attr("y", 22).attr("text-anchor", "end").style("font-size", "11px").style("fill", ink2).text("In shadow throughout");
  key.append("text").attr("y", -7).style("font-size", "11px").style("fill", ink2).text("Share of the exposure in shadow");
  [["Edge at the start of the exposure", "4 3"], ["Edge at the end", null]].forEach(([t, dash], i) => {     // the two edges, as each plate draws them
    const ex = narrow ? 0 : kw + 30, ey = narrow ? 40 + 16 * i : 16 * i - 2, eg = key.append("g").attr("transform", `translate(${ex},${ey})`);
    eg.append("line").attr("x2", 26).attr("stroke", edgeC).attr("stroke-width", 2).attr("stroke-dasharray", dash);
    eg.append("text").attr("x", 32).attr("y", 4).style("font-size", "11px").style("fill", ink2).text(t);
  });
  key.append("text").attr("x", narrow ? 0 : kw + 30).attr("y", narrow ? 80 : 36).style("font-size", "11px").style("fill", ink3).text("North up, east left");
  // the play button: every panel's exposure from start to end together, then back to the exposure as a whole
  d3.select("#plates-play").on("click", function () {
    const btn = d3.select(this).property("disabled", true), still = matchMedia("(prefers-reduced-motion: reduce)").matches;
    const fadeIn = still ? 0 : 500, hold = still ? 1500 : 500, fadeOut = still ? 0 : 1200;
    // the moving layer (above the still one) fades in over the whole exposure, plays, then fades out to reveal it again
    played.forEach(q => { q.set(still ? 1 : 0); q.moving.interrupt().style("display", null).style("opacity", 0).transition().duration(fadeIn).style("opacity", 1); });
    const done = () => {
      played.forEach(q => q.moving.transition().duration(fadeOut).style("opacity", 0).on("end", () => q.moving.style("display", "none")));
      setTimeout(() => btn.property("disabled", false), fadeOut);
    };
    if (still) { setTimeout(done, hold); return; }
    const t0 = performance.now() + fadeIn, dur = 4000;
    const tick = now => { const x = Math.min(Math.max((now - t0) / dur, 0), 1); played.forEach(q => q.set(x)); x < 1 ? requestAnimationFrame(tick) : setTimeout(done, hold); };
    requestAnimationFrame(tick);
  });
  const notes = notShown(Z);
  notes.forEach((t, i) => s.append("text").attr("x", x0).attr("y", (narrow ? keyY + 100 : keyY + 54) + 15 * i).style("font-size", "11px").style("fill", ink3).text(t));
  if (notes.length > 1) s.attr("viewBox", `0 0 ${W} ${keyY + (narrow ? 110 : 64) + 15 * (notes.length - 1)}`);
}

// the ratio plate by plate, each with its range from counting alone, and the common fit
function forestChart(s, Z) {
  // each plate's own ratio (counting range only), the plates in date order along the bottom, and the common fit across
  const rows = shown_(Z).slice().sort((a, b) => (a.date || "").localeCompare(b.date || "")), F = Z.fit;
  const lows = rows.flatMap(r => [r.r, r.lo].filter(v => v > 0)), highs = rows.flatMap(r => [r.r, r.hi].filter(v => v != null));
  const YMIN = Math.min(0.1, 0.8 * d3.min(lows)), YMAX = Math.max(12, 1.25 * d3.max(highs));   // the axis holds every plate
  const W = 760, m = { l: 52, r: 150, t: 30, b: 64 }, H0 = 306, notes = notShown(Z), H = H0 + 8 + 15 * notes.length;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  notes.forEach((t, i) => s.append("text").attr("x", 10).attr("y", H0 + 12 + 15 * i).style("font-size", "11px").style("fill", V("--ink-3")).text(t));
  const x = d3.scaleBand().domain(rows.map(r => r.plate)).range([m.l, W - m.r]).padding(0.3), cx = q => x(q) + x.bandwidth() / 2;
  const yl = d3.scaleLog().domain([YMIN, YMAX]).range([H0 - m.b, m.t]), y = v => yl(Math.min(YMAX, Math.max(YMIN, v)));
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(yl).tickValues(pow2Ticks(YMIN, YMAX)).tickFormat(pow2Label).tickSize(3)).call(axisStyle);
  s.append("text").attr("transform", `translate(14,${(m.t + H0 - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("In shadow, observed ÷ expected");
  const xa = s.append("g").attr("transform", `translate(0,${H0 - m.b})`).call(d3.axisBottom(x).tickSize(0).tickSizeOuter(0)).call(axisStyle);
  xa.selectAll("text").attr("transform", "rotate(-60)").attr("text-anchor", "end").attr("dx", "-0.4em").attr("dy", "0.5em");
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H0 - 6).attr("text-anchor", "middle").style("fill", V("--ink-2")).text("Plates in order of date");
  const band = s.append("rect").attr("x", m.l).attr("width", W - m.r - m.l).attr("y", y(F.count_hi)).attr("height", y(F.count_lo) - y(F.count_hi)).attr("fill", V("--ours")).attr("fill-opacity", 0.18);
  s.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(1)).attr("y2", y(1)).attr("stroke", V("--ink")).attr("stroke-width", 1.2);
  const below = F.r <= 1;                                                  // the two labels on either side of the two lines
  s.append("text").attr("class", "ann").attr("x", W - m.r + 8).attr("y", below ? Math.min(y(1), y(F.r)) - 6 : Math.max(y(1), y(F.r)) + 14).style("fill", V("--ink-2")).text("No effect");
  s.append("line").attr("x1", m.l).attr("x2", W - m.r).attr("y1", y(F.r)).attr("y2", y(F.r)).attr("stroke", V("--ours")).attr("stroke-width", 2);
  const fl = below ? Math.max(y(1), y(F.r)) + 15 : Math.min(y(1), y(F.r)) - 22;
  s.append("text").attr("x", W - m.r + 8).attr("y", fl).style("font-weight", 600).style("fill", V("--ours")).text("All plates");
  s.append("text").attr("class", "ann").attr("x", W - m.r + 8).attr("y", fl + 17).style("fill", V("--ours")).text(`${F.r.toFixed(2)} (${F.count_lo.toFixed(2)} to ${F.count_hi.toFixed(2)})`);
  rows.forEach(r => {
    const X = cx(r.plate), hi = r.hi == null ? YMAX : r.hi, c = V("--ink-2");
    s.append("line").attr("x1", X).attr("x2", X).attr("y1", y(Math.max(r.lo, YMIN))).attr("y2", y(hi)).attr("stroke", V("--ink-3")).attr("stroke-width", 1.6);
    if (r.lo > YMIN) s.append("line").attr("x1", X - 4).attr("x2", X + 4).attr("y1", y(r.lo)).attr("y2", y(r.lo)).attr("stroke", V("--ink-3")).attr("stroke-width", 1.6);
    if (hi < YMAX) s.append("line").attr("x1", X - 4).attr("x2", X + 4).attr("y1", y(hi)).attr("y2", y(hi)).attr("stroke", V("--ink-3")).attr("stroke-width", 1.6);
    else s.append("path").attr("d", d3.symbol(d3.symbolTriangle, 30)()).attr("transform", `translate(${X},${y(YMAX) + 2})`).attr("fill", V("--ink-3"));
    if (r.r > YMAX || r.r < YMIN)                                            // beyond the axis: at its end, its value beside it
      s.append("text").attr("x", X).attr("y", r.r > YMAX ? y(YMAX) - 9 : y(YMIN) - 9).attr("text-anchor", "middle").style("font-size", "10.5px").style("fill", c).text(r.r > YMAX ? `${r.r.toFixed(0)}` : r.r.toFixed(2));
    s.append("circle").attr("cx", X).attr("cy", y(r.r)).attr("r", 4.5).attr("fill", c)
      .call(sel => hov(sel, () => `<b>Plate ${r.plate}</b> · ${r.date || ""}<br>${pctText(r.shadow)}% of the plate in shadow over the exposure<br>${r.Ow.toFixed(1)} of ${fmt(r.n)} transients in shadow, each by its share of the exposure there, with ${E0(r).toFixed(1)} expected with no effect<br>ratio ${r.r.toFixed(2)} (${r.lo.toFixed(2)} to ${r.hi == null ? "–" : r.hi.toFixed(2)}, counting only)`));
  });
}


// the ratio over the shuffles of the plate times, stepped through: the shuffles built up with their middle 68%, the
// measured ratio, the shuffles as low (shaded), and a summary with the 68% range of true ratios, shown on load
// the shadow test's model on an invented plate, as four panels: (a) the artefacts, placed as transients lie on the
// plates the shadow never reaches (the map tinted beneath); (b) glints added, spread evenly; (c) Earth's shadow, the
// glints there lost (faint rings) and the artefacts kept; (d) the test without and with the glints: the transients
// observed in shadow against those expected. A sharp edge (the real shadow of the plate holding most transients, at
// mid-exposure) keeps every count whole; the moving edge is in the plates' chart
function modelChart(s, M) {
  const half = 3, ink = V("--ink"), ink2 = V("--ink-2"), ink3 = V("--ink-3"), glint = V("--glint"), red = V("--emulsion");
  const narrow = s.node().getBoundingClientRect().width < 560;
  const W = 760, gap = 20, lab = 40, S = (W - 2 * gap) / 3, DH = 116, H = S + lab + gap + DH;   // three plates in a row, (d) below
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const k = S / (2 * half), px = v => (half - v) * k, py = v => (half - v) * k;     // east to the left, north up
  const nc = Math.round(Math.sqrt(M.acc.length)), cv = document.createElement("canvas"); cv.width = cv.height = nc;
  const ctx = cv.getContext("2d"), amax = d3.max(M.acc);
  for (let r = 0; r < nc; r++) for (let c = 0; c < nc; c++) {
    const a = M.acc[(nc - 1 - r) * nc + (nc - 1 - c)] / amax; ctx.fillStyle = `rgba(235,238,242,${0.42 * a})`; ctx.fillRect(c, r, 1, 1);
  }
  const tintURL = cv.toDataURL(), ring = d3.line()(M.ring.map(q => [px(q[0]), py(q[1])])) + "Z";
  const P = M.points.map(([x, y, isG, sh]) => ({ x: px(x), y: py(y), isG, sh }));
  const A = P.filter(d => !d.isG), G = P.filter(d => d.isG);
  function plate(i, x0, y0, title, sub, draw) {
    const g = s.append("g").attr("transform", `translate(${x0},${y0})`), cid = `model-clip-${i}`;
    g.append("text").attr("y", 15).style("font-size", "13.5px").style("font-weight", 600).style("fill", ink).text(title);
    g.append("text").attr("y", 32).style("font-size", "11.5px").style("fill", ink2).text(sub);
    const p = g.append("g").attr("transform", `translate(0,${lab})`);
    p.append("clipPath").attr("id", cid).append("rect").attr("width", S).attr("height", S);
    const inner = p.append("g").attr("clip-path", `url(#${cid})`);
    inner.append("rect").attr("width", S).attr("height", S).attr("fill", V("--plate-day"));
    inner.append("image").attr("href", tintURL).attr("width", S).attr("height", S).attr("preserveAspectRatio", "none");
    draw(inner);
    p.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", ink3).attr("stroke-width", 0.8);
  }
  const dots = (g, D, fill, r) => g.append("g").selectAll("circle").data(D).join("circle").attr("cx", d => d.x).attr("cy", d => d.y).attr("r", r)
    .attr("fill", fill).attr("stroke", "#0d1013").attr("stroke-width", 0.5);
  const shadow = g => { g.append("path").attr("d", ring).attr("fill", V("--plate-night")).attr("fill-opacity", 0.75);
    g.append("path").attr("d", ring).attr("fill", "none").attr("stroke", SHADOW).attr("stroke-width", 2); };
  const nG = G.length, nGs = G.filter(d => d.sh).length;
  plate(0, 0, 0, "(a) Artefacts", `${A.length}, denser towards the edges`, g => dots(g, A, "#ffffff", 1.5));
  plate(1, S + gap, 0, "(b) Glints added, evenly", `${nG}, ${Math.round(100 * nG / (A.length + nG))}% of the transients`,
    g => { dots(g, A, "#ffffff", 1.5); dots(g, G, glint, 2.6); });
  plate(2, 2 * (S + gap), 0, "(c) Earth’s shadow", `${nGs} glints never appear (rings)`, g => {
    shadow(g); dots(g, A, "#ffffff", 1.5); dots(g, G.filter(d => !d.sh), glint, 2.6);
    g.append("g").selectAll("circle").data(G.filter(d => d.sh)).join("circle").attr("cx", d => d.x).attr("cy", d => d.y).attr("r", 2.6)
      .attr("fill", "none").attr("stroke", glint).attr("stroke-width", 1).attr("stroke-opacity", 0.8);
  });
  // (d) the test without and with the glints: observed and expected in shadow, and their ratio
  const d0 = { O: A.filter(d => d.sh).length, n: A.length, f: M.F }, d1 = { O: A.filter(d => d.sh).length, n: A.length + nG - nGs, f: M.F_glints };
  [d0, d1].forEach(c => { c.E = c.n * c.f; });
  const gd = s.append("g").attr("transform", `translate(0,${S + lab + gap})`);
  gd.append("text").attr("y", 15).style("font-size", "13.5px").style("font-weight", 600).style("fill", ink).text("(d) Transients in shadow");
  [["observed", red], ["expected", ink3]].forEach(([nm, col], j_) => {                // the key, beside the title
    const kx_ = 190 + j_ * 100; gd.append("rect").attr("x", kx_).attr("y", 5).attr("width", 12).attr("height", 10).attr("rx", 2).attr("fill", col);
    gd.append("text").attr("x", kx_ + 17).attr("y", 14.5).style("font-size", "12px").style("fill", ink2).text(nm);
  });
  const lx = 150, bw = 330, bx = d3.scaleLinear().domain([0, d3.max([d0, d1], c => Math.max(c.O, c.E))]).range([0, bw]);
  [[d0, "Artefacts only (a)", ink], [d1, "With glints (c)", glint]].forEach(([c, label, col], i) => {
    const row = gd.append("g").attr("transform", `translate(0,${36 + i * 40})`);
    row.append("text").attr("y", 14).style("font-size", "13px").style("fill", ink).text(label);
    row.append("rect").attr("x", lx).attr("y", 2).attr("height", 9).attr("rx", 2).attr("width", bx(c.O)).attr("fill", red);
    row.append("rect").attr("x", lx).attr("y", 13).attr("height", 9).attr("rx", 2).attr("width", bx(c.E)).attr("fill", ink3);
    const tx = row.append("text").attr("x", lx + bw + 16).attr("y", 16).style("font-size", "13px").style("fill", ink2);
    tx.append("tspan").text(`${c.O} ÷ ${c.E.toFixed(1)} = `);
    tx.append("tspan").style("font-weight", 700).style("fill", col).text((c.O / c.E).toFixed(2));
  });
  if (narrow) s.attr("viewBox", `0 0 ${W} ${H}`);
}

function ratioNullChart(s, Z, steps) {
  // stepped: the spread with no effect built up; the measured value; the spread as far towards a shortfall (shaded), the
  // summary shown on load. Z.right: a shortfall lies to the right (the shortfall statistic), else to the left (a ratio)
  const W = 760, H = 340, m = { l: 20, r: 20, t: 70, b: 46 };
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const B0 = Z.counts.map((c, i) => ({ x0: Z.edges[i], x1: Z.edges[i + 1], c })), n0 = d3.sum(Z.counts), r = Z.r, pLow = Math.round(100 * Z.p_low), R = !!Z.right;
  const nf = v => R ? v.toFixed(1) : v.toFixed(2);
  const Lm = Z.limit, hL = Lm ? Lm.counts.map(c => c * n0 / d3.sum(Lm.counts)) : [];   // the limit's fits, scaled to the same total
  const dom = [Z.edges[0], Z.edges[Z.edges.length - 1]];
  const x = (Z.log ? d3.scaleLog() : d3.scaleLinear()).domain(dom).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, Math.max(d3.max(B0, d => d.c), d3.max(hL) || 0) * 1.08]).range([H - m.b, m.t]);
  const ink = V("--ink"), grey = V("--ink-3"), red = V("--emulsion");
  const ax = Z.log ? d3.axisBottom(x).tickValues(pow2Ticks(dom[0], dom[1])).tickFormat(pow2Label)
                   : d3.axisBottom(x).ticks(12);
  s.append("g").attr("transform", `translate(0,${H - m.b})`).call(ax.tickSizeOuter(0)).call(axisStyle);
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 8).attr("text-anchor", "middle").text(Z.xlabel || "Fitted ratio, r");
  const t1 = s.append("text").attr("x", m.l).attr("y", 20).style("font-size", "13.5px");
  const bars = (g, B, h) => g.selectAll("rect").data(B).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1)).attr("rx", 1.5)
    .attr("y", (d, i) => y(h[i])).attr("height", (d, i) => y(0) - y(h[i]));
  const nullG = s.append("g").attr("fill", grey).attr("fill-opacity", 0.55);
  // the limit: the shuffles' fits with the limiting shortfall put in, scaled to the same total, the part as high as the
  // measured ratio stronger
  const limG = s.append("g");
  if (Lm) {
    limG.selectAll("rect").data(B0).join("rect").attr("x", d => x(d.x0) + 0.5).attr("width", d => Math.max(0, x(d.x1) - x(d.x0) - 1)).attr("rx", 1.5)
      .attr("y", (d, i) => y(hL[i])).attr("height", (d, i) => y(0) - y(hL[i])).attr("fill", V("--ours")).attr("fill-opacity", d => (R ? d.x1 <= r + 1e-9 : d.x0 >= r - 1e-9) ? 0.9 : 0.4);
    const pk = d3.maxIndex(hL), px = x(Math.sqrt(B0[pk].x0 * B0[pk].x1));   // labelled over their own peak (the fit sits a little low),
    const lab = limG.append("text").attr("y", y(hL[pk]) - 10).style("font-size", "12px").style("font-weight", 600)   // clear of the
      .style("fill", V("--ours")).text(Z.replicas ? `Replicas, ${Math.round(100 * Lm.g0)}% glints` : `Fits with a true ratio of ${Lm.r_true.toFixed(2)}`);   // measured line
    const wl = lab.node().getComputedTextLength(), xr = x(r);
    if (Math.abs(px - xr) < wl / 2 + 8) lab.attr("x", px < xr ? Math.max(m.l + wl, xr - 8) : Math.min(W - m.r - wl, xr + 8)).attr("text-anchor", px < xr ? "end" : "start");
    else lab.attr("x", Math.min(Math.max(px, m.l + wl / 2), W - m.r - wl / 2)).attr("text-anchor", "middle");
  }
  const cid = `${fid_(steps)}-clip`, clip = s.append("defs").append("clipPath").attr("id", cid).append("rect").attr("y", 0).attr("height", H).attr("x", x(r)).attr("width", 0);
  const lowG = s.append("g").attr("fill", red).attr("fill-opacity", 0.5).attr("clip-path", `url(#${cid})`);
  const obs = s.append("g");
  obs.append("line").attr("x1", x(r)).attr("x2", x(r)).attr("y1", m.t - 18).attr("y2", H - m.b).attr("stroke", red).attr("stroke-width", 2.4);
  obs.append("text").attr("x", x(r)).attr("y", m.t - 24).attr("text-anchor", "middle").style("font-size", "12px").style("font-weight", 600).style("fill", red).text(`Measured: ${nf(r)}`);
  const parts = { obs: [obs], low: [lowG], lim: [limG] };
  Object.values(parts).flat().forEach(e => e.style("opacity", 0));
  const full0 = B0.map(d => d.c), zero = B0.map(() => 0);
  const slides = Z.replicas ? [
    { obs: 0, low: 0, build: true, text: [`With no glints, replicas of the measurement on other plates give a typical ratio of ${nf(Z.null.median)}.`, ink] },
    { obs: 1, low: 0, text: [`The measured ratio is ${nf(r)}.`, red] },
    { obs: 1, low: 1, shade: true, text: [`${pLow}% of the replicas give ${nf(r)} or ${R ? "more" : "lower"} (shaded).`, red] },
  ].concat(Lm ? [{ obs: 1, low: 0, lim: 1, dim: true, text: [`If ${Math.round(100 * Lm.g0)}% of transients were glints, only ${Math.round(100 * Lm.share_as_high)}% of replicas would give ${nf(r)} or more.`, V("--ours")] }] : []) : [
    { obs: 0, low: 0, build: true, text: [`With no shadow effect, shuffles of the plate times give a typical ratio of ${Z.null.median.toFixed(2)}.`, ink] },
    { obs: 1, low: 0, text: [`The measured ratio is ${r.toFixed(2)}.`, red] },
    { obs: 1, low: 1, shade: true, text: [`${pLow}% of the shuffles give a ratio of ${r.toFixed(2)} or lower (shaded).`, red] },
  ].concat(Lm ? [{ obs: 1, low: 0, lim: 1, dim: true, text: [`With a true ratio of ${Lm.r_true.toFixed(2)}, only ${Math.round(100 * Lm.share_as_high)}% of shuffles give ${r.toFixed(2)} or more, so lower ratios are excluded at 95%.`, V("--ours")] }] : []);
  const last = slides.length - 1, op = { obs: 0, low: 0, lim: 0 };
  let cur = -1, timers = [];
  const at = (ms, f) => timers.push(setTimeout(f, ms));
  const show = (k, dur) => {
    timers.forEach(clearTimeout); timers = []; s.selectAll("*").interrupt();
    const sl = slides[k], out = dur ? 250 : 0;
    for (const [key, els] of Object.entries(parts)) if (op[key] && !sl[key]) { els.forEach(e => e.transition().duration(out).style("opacity", 0)); op[key] = 0; }
    t1.style("opacity", 0); nullG.transition().duration(dur ? 400 : 0).attr("fill-opacity", sl.dim ? 0.18 : 0.55);
    if (!sl.shade) clip.attr("x", x(r)).attr("width", 0);
    const fx = R ? x(r) : m.l, fw = R ? W - m.r - x(r) : x(r) - m.l;          // the shaded region, as far towards a shortfall
    if (sl.build && dur) bars(nullG, B0, zero);
    at(out, () => {
      bars(lowG, B0, full0);
      t1.text(sl.text[0]).style("fill", sl.text[1]).transition().duration(dur ? 400 : 0).style("opacity", 1);
      for (const [key, els] of Object.entries(parts)) if (sl[key] && !op[key]) { els.forEach(e => e.style("opacity", 0).transition().duration(dur ? 450 : 0).style("opacity", 1)); op[key] = 1; }
      if (sl.shade) {
        const start = dur && cur !== last;
        clip.attr("x", start ? x(r) : fx).attr("width", start ? 0 : fw).transition().delay(dur ? 450 : 0).duration(dur ? 1000 : 0).attr("x", fx).attr("width", fw);
      }
      if (!(sl.build && dur)) return bars(nullG, B0, full0);
      const cum = d3.cumsum(full0), tot = cum[cum.length - 1]; let seed = 11; const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
      const h = zero.slice(); let done = 0; const N = Math.round(n0);
      for (let k2 = 1; k2 <= 40; k2++) at(200 + k2 * 60, () => { const target = Math.round(N * (k2 / 40) ** 2.2);
        for (; done < target; done++) h[d3.bisectRight(cum, rnd() * tot)] += 1; bars(nullG, B0, h); });
      at(200 + 41 * 60, () => bars(nullG, B0, full0));
    });
    cur = k;
  };
  const ctl = d3.select("#" + steps), back = ctl.select(".back"), next = ctl.select(".next"), pos = ctl.select(".pos");
  const go = (k, instant) => {
    show(k, instant || matchMedia("(prefers-reduced-motion: reduce)").matches ? 0 : 1);
    back.property("disabled", k === 0);
    next.text(k === last ? "Step through ›" : "Next ›");
    pos.text(`${k + 1} of ${last + 1}`);
  };
  back.on("click", () => go(cur - 1));
  next.on("click", () => go(cur === last ? 0 : cur + 1));
  go(last, true);
}

// a null distribution as bars, the part at or beyond the measured value darker, split at the red line
function permPanel(s, q, y0, o) {
  const W = 760, m = { l: 52, r: 16 }, H = 200, top = y0 + 34, bot = y0 + H - 44, tot = d3.sum(q.counts);
  const x = d3.scaleLinear().domain([q.edges[0], q.edges[q.edges.length - 1]]).range([m.l, W - m.r]);
  const y = d3.scaleLinear().domain([0, d3.max(q.counts) / tot * 1.1]).range([bot, top + 52]);   // room above the bars for the labels
  s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(3).tickFormat(d3.format(".0%")).tickSize(3)).call(axisStyle);
  const bars = (g, op) => g.selectAll("rect").data(q.counts).join("rect").attr("x", (c, j) => x(q.edges[j]) + 0.75).attr("width", (c, j) => Math.max(0.5, x(q.edges[j + 1]) - x(q.edges[j]) - 1.5))
    .attr("y", c => y(c / tot)).attr("height", c => bot - y(c / tot)).attr("fill", V("--ink-3")).attr("fill-opacity", op);
  bars(s.append("g"), 0.45);
  s.append("clipPath").attr("id", o.id).append("rect").attr("x", x(o.value)).attr("width", W - x(o.value)).attr("y", 0).attr("height", y0 + H);
  bars(s.append("g").attr("clip-path", `url(#${o.id})`), 0.9);
  s.append("line").attr("x1", x(o.value)).attr("x2", x(o.value)).attr("y1", top).attr("y2", bot).attr("stroke", V("--emulsion")).attr("stroke-width", 2.5);
  s.append("text").attr("class", "ann").attr("x", m.l).attr("y", y0 + 16).style("font-weight", 600).text(o.title);
  const right = x(o.value) < W * 0.62, lab2 = (yy, txt) => s.append("text").attr("class", "ann").attr("x", x(o.value) + (right ? 8 : -8)).attr("y", yy)
    .attr("text-anchor", right ? "start" : "end").style("fill", V("--emulsion")).text(txt);
  lab2(top + 14, o.line1);
  const l2 = lab2(top + 30, o.line2 + (o.sd ? ` (${o.sd})` : "")), L = l2.node().getComputedTextLength(), X = x(o.value) + (right ? 8 : -8);
  if (o.sd && (right ? X + L > W - 6 : X - L < 6)) { l2.text(o.line2); lab2(top + 46, `(${o.sd})`); }   // too long: the standard deviations below
  s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).ticks(8).tickFormat(d => d + "×").tickSizeOuter(0)).call(axisStyle);
  if (o.xlab) s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", bot + 34).attr("text-anchor", "middle").text(o.xlab);
}

// a score's distribution for each length of line, catalogue against the copies, with the difference in standard deviations
function distChart(s, T, o) {
  const W = 760, m = { l: 52, r: 10, t: 50, b: 52 }, gx = 12;
  const rows = [{ key: "dist", h: 160, lab: o.rowlab }, { key: "z", h: 110, lab: "How rarely the random catalogues hold as many track candidates, or as few (in standard deviations)" }];
  const gy = 34, H = m.t + d3.sum(rows, r => r.h) + gy * (rows.length - 1) + m.b;
  s.attr("viewBox", `0 0 ${W} ${H}`);
  const pw = (W - m.l - m.r - 3 * gx) / 4, step = T.edges[1] - T.edges[0], ymax = 3 * d3.max(T.points, (np, ci) => d3.max([...T.real[ci], ...T.copies_mean[ci]]));
  let top = m.t;
  rows.forEach((row, ri) => {
    const bot = top + row.h;
    s.append("text").attr("class", "ann").attr("x", m.l).attr("y", ri === 0 ? 14 : top - 8).style("font-weight", 600).text(row.lab);   // the chart's title first, then the panels'
    T.points.forEach((np, ci) => {
      const x0 = m.l + ci * (pw + gx), x = d3.scaleLinear().domain(o.dom).range([x0, x0 + pw]).clamp(true);
      if (ri === 0) s.append("text").attr("x", x0 + pw / 2).attr("y", top - 8).attr("text-anchor", "middle").style("font-weight", 600).text(`Tracks of ${np} transients`);
      s.append("rect").attr("x", x0).attr("y", top).attr("width", pw).attr("height", row.h).attr("fill", "none").attr("stroke", V("--rule"));
      const bins = T.edges.slice(0, -1).map((l, k) => ({ l, h: l + step, n: T.real[ci][k], b: T.copies_mean[ci][k], z: T.z[ci][k] })).filter(q => q.l >= o.dom[0] - 1e-9 && q.h <= o.dom[1] + 1e-9);
      if (row.key === "dist") {
        const fl = o.floor, y = d3.scaleLog().domain([fl, ymax]).range([bot, top]).clamp(true);   // one scale for every panel, labelled on the first
        const stepPath = acc => d3.line().x(q => x(q[0])).y(q => y(Math.max(fl, q[1]))).curve(d3.curveStepAfter)(bins.map(q => [q.l, acc(q)]).concat([[bins[bins.length - 1].h, acc(bins[bins.length - 1])]]));
        s.append("path").attr("d", stepPath(q => q.b) + `L${x(bins[bins.length - 1].h)},${bot}L${x(bins[0].l)},${bot}Z`).attr("fill", V("--ink-3")).attr("fill-opacity", 0.35);
        const pts = bins.filter(q => q.n > 0);
        // Poisson 68% intervals (Gehrels 1986)
        s.append("g").selectAll("line").data(pts).join("line").attr("x1", q => x(q.l + step / 2)).attr("x2", q => x(q.l + step / 2))
          .attr("y1", q => y(Math.max(fl, q.n * (1 - 1 / (9 * q.n) - 1 / (3 * Math.sqrt(q.n))) ** 3))).attr("y2", q => y(q.n + 1 + Math.sqrt(q.n + 0.75))).attr("stroke", V("--emulsion")).attr("stroke-width", 1.4);
        s.append("g").selectAll("circle").data(pts).join("circle").attr("cx", q => x(q.l + step / 2)).attr("cy", q => y(q.n)).attr("r", 2.6).attr("fill", V("--emulsion"));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues(o.yticks.filter(v => v <= ymax)).tickFormat(d3.format(",~g")).tickSize(3)).call(axisStyle);
      } else {
        const y = d3.scaleLinear().domain([-3.5, 3.5]).range([bot, top]).clamp(true), bw = x(o.dom[0] + step) - x(o.dom[0]);
        s.append("g").selectAll("rect").data(bins).join("rect")
          .attr("x", q => x(q.l) + 0.5).attr("width", Math.max(1, bw - 1)).attr("y", q => Math.min(y(0), y(q.z))).attr("height", q => Math.abs(y(q.z) - y(0)))
          .attr("fill", V("--emulsion")).attr("fill-opacity", 0.85)
          .call(sel => hov(sel, q => `Tracks of ${np} transients, ${o.what} ${(10 ** q.l).toPrecision(2)} to ${(10 ** q.h).toPrecision(2)}<br>Catalogue ${fmt(q.n)}; copies ${d3.format(",.2f")(q.b)} on average<br>as unusual among the copies as ${q.z >= 0 ? "+" : "−"}${Math.abs(q.z).toFixed(1)} standard deviations`));
        s.append("line").attr("x1", x0).attr("x2", x0 + pw).attr("y1", y(0)).attr("y2", y(0)).attr("stroke", V("--ink-3"));
        if (ci === 0) s.append("g").attr("transform", `translate(${x0},0)`).call(d3.axisLeft(y).tickValues([-3, 0, 3]).tickFormat(v => (v > 0 ? "+" : v < 0 ? "−" : "") + Math.abs(v)).tickSize(3)).call(axisStyle);
      }
      if (ri === rows.length - 1) s.append("g").attr("transform", `translate(0,${bot})`).call(d3.axisBottom(x).tickValues(o.ticks).tickFormat(v => o.fmt[String(v)]).tickSize(3).tickSizeOuter(0)).call(axisStyle);
    });
    top = bot + gy;
  });
  s.append("text").attr("x", (m.l + W - m.r) / 2).attr("y", H - 12).attr("text-anchor", "middle").text(o.xlab);
}





// the track score's method, stepped through on one candidate (drawn from analysis/alignment_example.py's record)
function methodChart(s, M) {
  // slides the reader steps through, in the order the method runs: 1 a real candidate; 2 the close-pair rule that forms
  // it; 3 the straightness test, which throws away candidates under the cut; 4 the four numbers measured on one that
  // passes; 5 the track score; and a summary, shown on load. Each slide declares how strongly each part shows
  const W = 760, H = 398, top = 92, ours = V("--ours"), grey = V("--ink-3"), ink = V("--ink"), ink2 = V("--ink-2"), rule = V("--rule");
  size(s, H);
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
  lab(m1, rp[0], rp[1] + 4, `k = ${r}`, { anchor: "middle", bold: true, col: ours });
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
  [[`k = ${r}`, "transients"], [`L = ${M.L.toFixed(1)}′`, "long"], [`w = ±${strip.toFixed(0)}″`, "wide"], [`N = ${N}`, "in the box"]].forEach(([v, t], i) => {
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
    { sky: 1, pair: dim, fit: dim, meas: 1, text: ["A kept candidate is measured: the number of transients in it (k), its length (L), its width (w, at least 1″),",
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

const DRAW = {
  // VASCO's catalogue, before any test (report/vasco_stats.py): transients per plate against the range counting noise allows
  "fig.vasco_spread": (s, F) => {
    const W = 760, H = 300, m = { l: 52, r: 16, t: 34, b: 46 }; size(s, H);
    const N = F.plates.map(p => p.N), bins = d3.range(0, 1000, 25).map(a => ({ lo: a, hi: a + 25, n: N.filter(v => v >= a && (a === 975 ? v <= 1000 : v < a + 25)).length }));
    const over = { lo: 1025, hi: 1075, n: N.filter(v => v > 1000).length, more: true };
    const x = d3.scaleLinear().domain([0, 1080]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, Math.ceil(d3.max(bins.concat(over), b => b.n) / 10) * 10]).range([H - m.b, m.t]);
    s.append("rect").attr("x", x(F.lo)).attr("width", x(F.hi) - x(F.lo)).attr("y", m.t - 14).attr("height", H - m.b - m.t + 14).attr("fill", V("--theirs-soft"));
    s.append("g").selectAll("rect.b").data(bins.concat([over])).join("rect").attr("x", b => x(b.lo) + 0.75).attr("width", b => x(b.hi) - x(b.lo) - 1.5).attr("rx", 1.5)
      .attr("y", b => y(b.n)).attr("height", b => y(0) - y(b.n)).attr("fill", V("--ours"))
      .call(sel => hov(sel, b => b.more ? `${b.n} plates hold more than 1,000 transients, up to ${fmt(d3.max(N))}` : `${b.n} plates hold ${fmt(b.lo)} to ${fmt(b.hi - 1)} transients`));
    s.append("line").attr("x1", x(F.mean)).attr("x2", x(F.mean)).attr("y1", m.t - 14).attr("y2", H - m.b).attr("stroke", V("--ink")).attr("stroke-width", 1.5).attr("stroke-dasharray", "4 3");
    ann(s, x(F.hi) + 8, m.t - 2, `Average (dashed line): ${Math.round(F.mean)}`);
    ann(s, x(F.hi) + 8, m.t + 16, `Shaded: the range counting noise allows, ${Math.round(F.lo)} to ${Math.round(F.hi)}`);
    ann(s, x(F.hi) + 8, m.t + 34, `${F.in_band} of ${N.length} plates fall in this range`);
    s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 200, 400, 600, 800]).tickFormat(d3.format(",")).tickSizeOuter(0)).call(axisStyle);
    ann(s, x(1050), H - m.b + 15, "Over 1,000", { anchor: "middle", fill: V("--ink-3") });
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickSize(3)).call(axisStyle);
    ann(s, (m.l + W - m.r) / 2, H - 8, "Transients per plate", { anchor: "middle" });
    s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Number of plates");
  },
  "fig.vasco_edge_map": (s, F) => {
    // the density across the plate, each plate weighted equally, over the density near the centre: transients, and VASCO's
    // control positions; north up, east left
    const S = 300, top = 34, xs = [52, 384], M = F, n = M.transients.length, c = S / n; size(s, top + S + 40);
    const col = d3.scaleSequential(d3.interpolateMagma).domain([0, 3.5]).clamp(true), ax = d3.scaleLinear().domain([M.edges[0], M.edges[n]]).range([0, S]);
    [["transients", "Transients"], ["control", "VASCO’s control positions"]].forEach(([k, lab], j) => {
      const g = s.append("g").attr("transform", `translate(${xs[j]},${top})`);
      ann(g, S / 2, -12, lab, { anchor: "middle", bold: true });
      g.append("g").selectAll("rect").data(M[k].flatMap((row, i) => row.map((v, jj) => ({ v, i, j: jj })))).join("rect")
        .attr("x", d => S - (d.j + 1) * c).attr("y", d => S - (d.i + 1) * c).attr("width", c + 0.3).attr("height", c + 0.3).attr("fill", d => col(d.v))
        .call(sel => hov(sel, d => `${d.v.toFixed(1)} times the density near the centre`));
      g.append("rect").attr("width", S).attr("height", S).attr("fill", "none").attr("stroke", V("--rule"));
      g.append("g").attr("transform", `translate(0,${S})`).call(d3.axisBottom(ax).tickValues([-3, -1.5, 0, 1.5, 3]).tickFormat(d => `${d}°`).tickSize(3)).call(axisStyle);
      ann(g, -6, 12, "N", { anchor: "end", fill: V("--ink-3") }); ann(g, -6, S / 2 + 4, "E", { anchor: "end", fill: V("--ink-3") });
    });
    const bx = 712, gid = "vasco-edge-grad", ly = d3.scaleLinear().domain([0, 3.5]).range([top + S, top]);
    s.append("defs").append("linearGradient").attr("id", gid).attr("x1", 0).attr("x2", 0).attr("y1", 1).attr("y2", 0)
      .selectAll("stop").data(d3.range(0, 1.001, 0.1)).join("stop").attr("offset", q => q).attr("stop-color", q => col(q * 3.5));
    s.append("rect").attr("x", bx).attr("y", top).attr("width", 12).attr("height", S).attr("fill", `url(#${gid})`);
    s.append("g").attr("transform", `translate(${bx + 12},0)`).call(d3.axisRight(ly).tickValues([0, 1, 2, 3]).tickFormat(d => `${d}×`).tickSize(3)).call(axisStyle);
    s.append("text").attr("transform", `translate(${bx - 8},${top + S / 2}) rotate(-90)`).attr("text-anchor", "middle").attr("fill", V("--ink-2")).text("Density relative to the centre");
  },
  "fig.vasco_edge": (s, F) => {
    // the same by distance from the centre: transients with their 68% range from resampling the plates, and the control positions
    const W = 760, H = 300, m = { l: 62, r: 16, t: 34, b: 46 }; size(s, H);
    const x = d3.scaleLinear().domain([0, 2.75]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, 3]).range([H - m.b, m.t]);
    const line = k => d3.line().x((v, i) => x(F.r[i])).y(v => y(v)).curve(d3.curveMonotoneX)(F[k]);
    s.append("path").attr("d", d3.area().x((v, i) => x(F.r[i])).y0((v, i) => y(F.lo[i])).y1((v, i) => y(F.hi[i])).curve(d3.curveMonotoneX)(F.r)).attr("fill", V("--ours-soft"));
    s.append("path").attr("d", line("control")).attr("fill", "none").attr("stroke", V("--ink-3")).attr("stroke-width", 2).attr("stroke-dasharray", "4 3");
    s.append("path").attr("d", line("transients")).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2.2);
    s.append("g").selectAll("circle").data(F.transients).join("circle").attr("cx", (v, i) => x(F.r[i])).attr("cy", v => y(v)).attr("r", 4).attr("fill", V("--ours")).attr("stroke", V("--surface")).attr("stroke-width", 1.5)
      .call(sel => hov(sel, (v, i) => `${(F.r[i] - 0.125).toFixed(2)}° to ${(F.r[i] + 0.125).toFixed(2)}° from the centre<br>Transients: ${v.toFixed(2)} times the density near the centre<br>Control positions: ${F.control[i].toFixed(2)}`));
    legend(s, m.l + 10, 16, [["Transients (band: 68% range over resampled plates)", "--ours", "dot"], ["VASCO’s control positions", "--ink-3", "line"]]);
    s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 0.5, 1, 1.5, 2, 2.5]).tickFormat(d => d + "°").tickSizeOuter(0)).call(axisStyle);
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(6).tickFormat(d => d + "×").tickSize(3)).call(axisStyle);
    ann(s, (m.l + W - m.r) / 2, H - 8, "Distance from the centre of the plate (the larger of the east–west and north–south offsets)", { anchor: "middle" });
    s.append("text").attr("transform", `translate(16,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Density relative to the centre");
  },
  "fig.vasco_crowd": (s, F) => {
    // every plate's count against its distance from the plane of the Milky Way, the typical plate's curve and the middle half of plates
    const W = 760, H = 360, m = { l: 56, r: 16, t: 34, b: 46 }, C = F.crowd_fit; size(s, H);
    const x = d3.scaleLinear().domain([0, 90]).range([m.l, W - m.r]), y = d3.scaleLog().domain([1, 5000]).range([H - m.b, m.t]).clamp(true);
    s.append("path").attr("d", d3.area().x((_, i) => x(C.b[i])).y0((_, i) => y(C.q25[i])).y1((_, i) => y(C.q75[i])).curve(d3.curveMonotoneX)(C.b)).attr("fill", V("--ink-3")).attr("fill-opacity", 0.2);
    s.append("g").selectAll("circle").data(F.plates).join("circle").attr("cx", p => x(p.b)).attr("cy", p => y(Math.max(p.N, 1))).attr("r", 3)
      .attr("fill", V("--ours")).attr("fill-opacity", 0.5).call(sel => hov(sel, p => `<b>${p.id}</b><br>${fmt(p.N)} transients, ${p.b}° from the Milky Way`));
    const ln = d3.line().x((_, i) => x(C.b[i])).y(v => y(v)).curve(d3.curveMonotoneX)(C.mid);
    s.append("path").attr("d", ln).attr("fill", "none").attr("stroke", V("--surface")).attr("stroke-width", 5).attr("pointer-events", "none");
    s.append("path").attr("d", ln).attr("fill", "none").attr("stroke", V("--ink")).attr("stroke-width", 2.2).attr("pointer-events", "none");
    legend(s, m.l + 10, 16, [["A plate", "--ours", "dot"], ["The typical plate", "--ink", "line"], ["The middle half of plates", "--ink-3", "band"]]);
    s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 15, 30, 45, 60, 75, 90]).tickFormat(v => v + "°").tickSizeOuter(0)).call(axisStyle);
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([1, 10, 100, 1000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
    ann(s, (m.l + W - m.r) / 2, H - 8, "Distance from the plane of the Milky Way", { anchor: "middle" });
    s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients on the plate");
  },
  "fig.vasco_pointing": (s, F) => {
    // every plate's count against the count the pointing model expects of it; the plates in eight groups by expected count, each
    // group's mean count with its 68% range; the line where the two are equal
    const W = 760, H = 380, m = { l: 56, r: 16, t: 34, b: 46 }; size(s, H);
    const [e0, e1] = d3.extent(F.plates, p => p.E), x = d3.scaleLog().domain([e0 / 1.15, e1 * 1.15]).range([m.l, W - m.r]), y = d3.scaleLog().domain([1, 5000]).range([H - m.b, m.t]).clamp(true);
    s.append("g").selectAll("circle").data(F.plates).join("circle").attr("cx", p => x(p.E)).attr("cy", p => y(Math.max(p.N, 1))).attr("r", 3)
      .attr("fill", V("--ink-3")).attr("fill-opacity", 0.45).call(sel => hov(sel, p => `<b>${p.id}</b><br>${fmt(p.N)} transients, ${fmt(Math.round(p.E))} expected from its pointing`));
    const d = x.domain(), diag = `M${x(d[0])},${y(d[0])}L${x(d[1])},${y(d[1])}`;
    s.append("path").attr("d", diag).attr("stroke", V("--surface")).attr("stroke-width", 5).attr("pointer-events", "none");
    s.append("path").attr("d", diag).attr("stroke", V("--ours")).attr("stroke-width", 2).attr("pointer-events", "none");
    const G = s.append("g").selectAll("g").data(F.calibration).join("g");
    G.append("line").attr("x1", c => x(c.E)).attr("x2", c => x(c.E)).attr("y1", c => y(c.lo)).attr("y2", c => y(c.hi)).attr("stroke", V("--ink")).attr("stroke-width", 2);
    G.append("circle").attr("cx", c => x(c.E)).attr("cy", c => y(c.N)).attr("r", 5).attr("fill", V("--ink")).attr("stroke", V("--surface")).attr("stroke-width", 1.5)
      .call(sel => hov(sel, c => `${c.n} plates expected ${fmt(Math.round(c.E))} on average<br>Actual average ${fmt(Math.round(c.N))} (68% range ${fmt(Math.round(c.lo))} to ${fmt(Math.round(c.hi))})`));
    legend(s, m.l + 10, 16, [["A plate", "--ink-3", "dot"], ["A group of about 80 plates, its average (bar: 68% range)", "--ink", "dot"], ["Actual equals expected", "--ours", "line"]]);
    s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([10, 20, 50, 100, 200, 500, 1000].filter(v => v >= d[0] && v <= d[1])).tickFormat(d3.format(",")).tickSizeOuter(0)).call(axisStyle);
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([1, 10, 100, 1000]).tickFormat(d3.format(",")).tickSize(3)).call(axisStyle);
    ann(s, (m.l + W - m.r) / 2, H - 8, "Transients the pointing model expects on the plate", { anchor: "middle" });
    s.append("text").attr("transform", `translate(14,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Transients on the plate");
  },
  "fig.sky": (s, F) => skyMap(s, F.plates, false, "fig.sky"),
  "fig.sky_shadow": (s, F) => skyMap(s, F.plates, true, "fig.sky_shadow"),
  "fig.shadow_plates": (s, F) => platesChart(s, F, "fig.shadow_plates"),
  "fig.shadow_forest": (s, F) => forestChart(s, F),
  "fig.glint_model": (s, F) => modelChart(s, F),
  "fig.shadow_null": (s, F) => ratioNullChart(s, F, "ratio-steps"),
  "fig.nuclear": (s, F) => {
    // the shifted calendars only, the stricter null (it keeps the tests' bunching); random nights' chance is quoted in the text
    const d = F; size(s, 220);
    permPanel(s, d.calendars, 0, { id: "nuc-cal", value: d.ratio, title: `The real calendar of tests shifted ${fmt(d.n_calendars)} times`,
      line1: `Plates near tests: ${d.ratio.toFixed(2)}×`, line2: `Shifted calendars match or exceed it ${shareText(d.p_calendar)}`, sd: d.p_calendar < 0.5 ? sdText(d.sigma) : null,   // a one-sided excess: no significance when chance does better
      xlab: "Transients on plates near tests relative to the other plates, each against its expected count" });
  },
  "fig.alignments": (s, F) => distChart(s, F, { dom: [-7, 0], ticks: [-6, -4, -2, 0], fmt: { "-6": "10⁻⁶", "-4": "10⁻⁴", "-2": "10⁻²", "0": "1" },
    floor: 0.1, yticks: [0.1, 1, 10, 100, 1000], what: "track score", rowlab: "Straight track candidates at each track score: the catalogue (red) and the random catalogues’ average (grey)",
    xlab: "Track score: the number of such tracks expected by chance (lower is rarer)" }),
  "fig.align_method": (s, F) => methodChart(s, F),
  "fig.align_direction": (s, F) => {
    // the share of artificial glint tracks the search finds at each angle from east-west, with its counting range
    const W = 480, H = 280, m = { l: 62, r: 24, t: 20, b: 50 };
    s.attr("viewBox", `0 0 ${W} ${H}`).style("width", `${100 * W / 760}%`).style("min-width", `${Math.round(620 * W / 760)}px`).style("display", "block").style("margin", "0 auto");   // narrower than the other charts, at the same text size
    const x = d3.scaleLinear().domain([-2, 92]).range([m.l, W - m.r]), y = d3.scaleLinear().domain([0, 1]).range([H - m.b, m.t]);
    const pts = F.angle.map((a, i) => [a, F.found[i], F.err[i]]);
    s.append("path").datum(pts).attr("d", d3.line().x(d => x(d[0])).y(d => y(d[1]))).attr("fill", "none").attr("stroke", V("--ours")).attr("stroke-width", 2);
    s.append("g").selectAll("line").data(pts).join("line").attr("x1", d => x(d[0])).attr("x2", d => x(d[0])).attr("y1", d => y(d[1] - d[2])).attr("y2", d => y(d[1] + d[2])).attr("stroke", V("--ours")).attr("stroke-width", 1.5);
    s.append("g").selectAll("circle").data(pts).join("circle").attr("cx", d => x(d[0])).attr("cy", d => y(d[1])).attr("r", 4).attr("fill", V("--ours"))
      .call(sel => hov(sel, d => `${d[0]}° from east–west<br>${(100 * d[1]).toFixed(0)}% of ${fmt(F.tracks)} tracks found`));
    s.append("g").attr("transform", `translate(0,${H - m.b})`).call(d3.axisBottom(x).tickValues([0, 30, 60, 90]).tickFormat(d => d + "°").tickSizeOuter(0)).call(axisStyle);
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).tickValues([0, 0.25, 0.5, 0.75, 1]).tickFormat(d3.format(".0%")).tickSize(3)).call(axisStyle);
    ann(s, x(0), H - m.b - 8, "east–west", { anchor: "start", fill: V("--ink-3") }); ann(s, x(90), H - m.b - 8, "north–south", { anchor: "end", fill: V("--ink-3") });
    ann(s, (m.l + W - m.r) / 2, H - 8, "Direction of the track, from east–west", { anchor: "middle" });
    s.append("text").attr("transform", `translate(16,${(m.t + H - m.b) / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Tracks found");
  },
  "fig.clock": (s, F) => {
    // every asteroid detection's offset along its track from the written time, over the survey; the survey's offset marked
    const W = 760, m = { l: 62, r: 20 }, h1 = 230, H = 40 + h1 + 50; size(s, H);
    const x = d3.scaleLinear().domain([F.edges[0], F.edges[F.edges.length - 1]]).range([m.l, W - m.r]);
    const y = d3.scaleLinear().domain([0, d3.max(F.counts) * 1.1]).range([40 + h1, 40]);
    s.append("g").selectAll("rect").data(F.counts.map((c, j) => ({ c, lo: F.edges[j], hi: F.edges[j + 1] }))).join("rect").attr("x", d => x(d.lo) + 0.4).attr("width", d => x(d.hi) - x(d.lo) - 0.8).attr("rx", 1)
      .attr("y", d => y(d.c)).attr("height", d => 40 + h1 - y(d.c)).attr("fill", V("--ink-3")).attr("fill-opacity", 0.6)
      .call(sel => hov(sel, d => `${d.lo} to ${d.hi} minutes: ${fmt(d.c)} detections`));
    s.append("line").attr("x1", x(0)).attr("x2", x(0)).attr("y1", 40).attr("y2", 40 + h1).attr("stroke", V("--rule")).attr("stroke-dasharray", "3 3");
    s.append("line").attr("x1", x(F.survey)).attr("x2", x(F.survey)).attr("y1", 30).attr("y2", 40 + h1).attr("stroke", V("--ours")).attr("stroke-width", 2);
    ann(s, x(F.survey) + 6, 36, `Survey offset: ${F.survey.toFixed(1).replace("-", "−")} minutes (median over plates)`, { fill: V("--ours") });
    s.append("g").attr("transform", `translate(${m.l},0)`).call(d3.axisLeft(y).ticks(5).tickFormat(d => d === 0 ? "" : d3.format(",")(d)).tickSize(3)).call(axisStyle);
    ann(s, m.l, 18, "When the asteroids show each plate was taken, against its written time", { bold: true });
    s.append("text").attr("transform", `translate(16,${40 + h1 / 2}) rotate(-90)`).attr("text-anchor", "middle").text("Asteroids");
    s.append("g").attr("transform", `translate(0,${40 + h1})`).call(d3.axisBottom(x).ticks(12).tickSizeOuter(0)).call(axisStyle);
    ann(s, (m.l + W - m.r) / 2, H - 8, "Minutes from the written time (negative: the plate was taken earlier than written)", { anchor: "middle" });
    (F.beyond || []).forEach((b, i) => ann(s, W - m.r, 64 + 18 * i, `Beyond the axis: ${b.plate}, ${b.dt.map(v => (v > 0 ? "+" : "−") + Math.abs(v).toFixed(1)).join(" and ")} minutes`, { anchor: "end", fill: V("--ink-2") }));
  },
};

// tables, drawn into <div class="chart-table" data-fig="...">
const TABLE = {
  "fig.dates_table": F => `<div class="tables"><div class="tablewrap"><table><thead><tr><th>Plate serials</th><th>Taken</th><th>Archive date</th></tr></thead><tbody>` +
    F.runs.map(r => `<tr><td class="m">${r.serials}</td><td>${r.taken}</td><td>${r.archive}</td></tr>`).join("") + "</tbody></table></div>" +
    `<div class="tablewrap"><table><thead><tr><th>Plate</th><th>Serial</th><th>Archive date</th><th>Correct date (UT)</th><th>Note</th></tr></thead><tbody>` +
    F.exceptions.map(e => `<tr><td class="m">${e.plate}</td><td class="m">${e.serial}</td><td class="m">${e.archive}</td><td class="m">${e.correct}</td><td>${e.note}</td></tr>`).join("") + "</tbody></table></div></div>",
};

// a plain table: F.head (column names) and F.rows (cells as text), the first column a label
const PLAIN = F => `<div class="tablewrap"><table><thead><tr>${F.head.map(h => `<th>${h}</th>`).join("")}</tr></thead><tbody>` +
  F.rows.map(r => `<tr>${r.map((c, i) => `<td${i ? ' class="m"' : ""}>${c}</td>`).join("")}</tr>`).join("") + "</tbody></table></div>";

let DATA = null;
function render() {
  if (!DATA) return;
  document.querySelectorAll("div.chart-table").forEach(el => { const F = DATA.figures[el.dataset.fig]; el.innerHTML = F ? (TABLE[el.dataset.fig] || PLAIN)(F) : "<p>Not yet produced.</p>"; });
  document.querySelectorAll("svg.chart-svg").forEach(el => {
    const id = el.dataset.fig, s = d3.select(el); s.selectAll("*").remove();
    if (!DATA.figures[id]) { size(s, 60); ann(s, 20, 30, `No data for ${id}.`); return; }
    DRAW[id](s, DATA.figures[id]);
  });
}
fetch("figures.json").then(r => r.json()).then(d => { DATA = d; render(); });
matchMedia("(prefers-color-scheme: dark)").addEventListener("change", render);
