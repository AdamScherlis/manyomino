"""Turn analysis/convergence.json into what the page and the README need:
chart data, text tokens, and analysis/convergence.md."""
import json, math, os, re

ROOT = os.path.join(os.path.dirname(__file__), "..")


def _stats(zs):
    a = [abs(z) for z in zs]
    return {"n": len(a), "rms": math.sqrt(sum(z * z for z in a) / len(a)), "max": max(a), "over2": sum(z > 2 for z in a), "over3": sum(z > 3 for z in a)}


def conv_report(kept_moves, hero_n, burn_n):
    fn = os.path.join(ROOT, "analysis", "convergence.json")
    if not os.path.exists(fn):
        return None, {}, []
    C = json.load(open(fn))
    T = {}; tests = []; md = []
    zs = lambda rows, key="z": [[r["name"], round(r[key], 2)] for r in rows]
    pct = lambda x: f"{100 * (math.exp(x) - 1):+.1f}%"

    # ---- the tests, in the order they are drawn
    for e in C["enum"]:
        tests.append({"label": ("every shape" if not e["p"] else f"weighted {e['p']:.3g}") + f"|{e['n']} cells", "z": zs(e["rows"]), "group": "enum",
                      "md": (f"all {e['count']:,} shapes of {e['n']} cells" + (f", weighted (1-p)^perimeter, p = {e['p']:.6g}" if e["p"] else ", uniform"),
                             f"{e['samples']:,} chain samples; precision {e['rel_precision']:.1e} sd")})
    for e in C["exact"]:
        n = int(re.search(r"n=(\d+)", e["tag"]).group(1)); p = float(re.search(r"p=([\d.]+)", e["tag"]).group(1))
        tests.append({"label": f"exact clusters|{n:,} cells", "z": zs(e["rows"]), "group": "exact", "n": n, "n_exact": e["n_exact"], "witness": e["auc_z"], "q": (e["q_n"], e["q_over2"], e["q_max"]),
                      "md": (f"{e['n_exact']:,} exact percolation clusters of {n:,} cells at p = {p:.6g}",
                             f"{e['chains']} chain{'s' if e['chains'] > 1 else ''}" + (f" (chi2/dof between them {e['between_chains']:.2f})" if e.get('between_chains') else "") + f"; learned test z = {e['auc_z']:+.2f}; quartiles of each distribution: {e['q_over2']} of {e['q_n']} beyond 2 s.e., largest {e['q_max']:.2f}")})
    for g in C["growth"]:
        ok = abs(g["count_dev"]) < 3 * g["count_dev_se"]
        tests.append({"label": f"grown population|{g['n']:,} cells" + ("" if ok else " *"), "z": zs(g["rows"], "t"), "kind": "chain" if ok else "strained", "group": "growth",
                      "n": g["n"], "count_dev": g["count_dev"], "count_se": g["count_dev_se"], "precision": g["median_precision"], "anc": sum(g["ancestors_m2"]) / len(g["ancestors_m2"]),
                      "md": (f"{g['replicas']} grown populations of {round(g['population'], -3):,} at {g['n']:,} cells, uniform" + ("" if ok else " (fails its own count check)"),
                             f"count check {pct(g['count_dev'])} ± {100 * g['count_dev_se']:.1f}%; {sum(g['ancestors_m2']) / len(g['ancestors_m2']):,.0f} distinct ancestors at n/2; precision {g['median_precision']:.3f} sd")})
    if C.get("symmetry"):
        sy = C["symmetry"]
        tests.append({"label": f"mirror symmetry|{sy[0]['n']:,} to {sy[-1]['n']:,}", "z": [[f"{r['name']} at {s['n']:,} cells", round(r["z"], 2)] for s in sy for r in s["rows"]], "group": "symmetry",
                      "md": (f"{sy[0]['features']} observables that are exactly zero by symmetry, at {len(sy)} sizes from {sy[0]['n']:,} to {sy[-1]['n']:,} cells", "two chains per size")})
    for s in C.get("starts", []):
        tests.append({"label": f"bar, square starts|{s['n']:,} cells", "z": zs(s["rows"]), "group": "starts", "n": s["n"], "witness": s["witness_z"],
                      "md": (f"chains started from a straight bar and a solid square of {s['n']:,} cells, after {s['burn']:,} moves per cell, against the ladder chains",
                             f"learned test z = {s['witness_z']:+.2f}")})
    if "control" in C:
        e = C["control"]
        tests.append({"label": "wrong on purpose|250 cells", "z": zs(e["rows"]), "kind": "control", "group": "control", "witness": e["auc_z"],
                      "md": ("control: chain at p = 0.51 against exact clusters at p = 0.50, 250 cells", f"learned test z = {e['auc_z']:+.2f}")})

    slow = []
    for s in C["slow"]:
        rg = next(x for x in s["single"] if x["name"] == "rg2"); s1 = s["slowest_single"][0]; mo = s["slowest_mode"]
        slow.append({"n": s["n"], "rg2": round(rg["tau"], 2), "rg2_err": round(rg["err"], 2), "single": round(s1["tau"], 2), "single_err": round(s1["err"], 2),
                     "single_name": s1["name"], "mode": round(mo["tau"], 2), "mode_err": round(mo["err"], 2), "samples": sum(s["samples"])})
    timed = [s for s in slow if s["single_err"] < 0.5 * s["single"]]         # runs long enough to time the slowest observable
    conv = {"tests": [{k: t[k] for k in ("label", "z", "kind") if k in t} for t in tests], "slow": timed}

    # ---- numbers for the text
    trusted = [t for t in tests if t.get("kind", "chain") == "chain"]
    allz = [z[1] for t in trusted for z in t["z"]]; st = _stats(allz)
    worst = max(((abs(z[1]), z[0], t["label"].replace("|", ", ")) for t in trusted for z in t["z"]))
    T["CV_NOBS"] = str(max(len(t["z"]) for t in tests if t["group"] == "enum"))
    T["CV_TESTS"] = str(len(trusted)); T["CV_COMPARISONS"] = f"{st['n']:,}"
    T["CV_OVER2"] = str(st["over2"]); T["CV_OVER2_PCT"] = f"{100 * st['over2'] / st['n']:.1f}"; T["CV_OVER3"] = str(st["over3"]); T["CV_EXPECT3"] = f"{0.0027 * st['n']:.0f}"
    T["CV_MAX"] = f"{worst[0]:.1f}"; T["CV_MAX_WHERE"] = f"{worst[1]}, {worst[2]}"
    en = {(e["n"], e["p"]): e for e in C["enum"]}
    for n in (14, 16):
        if (n, 0.0) in en: T[f"CV_COUNT{n}"] = f"{en[(n, 0.0)]['count']:,}"
    T["CV_ENUM_MAX"] = f"{max(e['z_max'] for e in C['enum']):.1f}"; T["CV_ENUM_N"] = f"{sum(e['features'] for e in C['enum']):,}"
    T["CV_ENUM_PREC"] = f"{100 * max(e['rel_precision'] for e in C['enum']):.2f}".rstrip("0")
    for t in tests:
        if t["group"] == "exact": T[f"CV_EXACT{t['n']}"] = f"{t['n_exact']:,}"
        if t["group"] == "growth":
            dev = 100 * (math.exp(t["count_dev"]) - 1)
            T[f"CV_GCOUNT{t['n']}"] = (f"{abs(dev):.0f}% {'low' if dev < 0 else 'high'} (&plusmn;&nbsp;{100 * t['count_se']:.0f}%)" if t.get("kind") == "strained"
                                       else f"{abs(dev):.1f}%&nbsp;&plusmn;&nbsp;{100 * t['count_se']:.1f}%")
            T[f"CV_GPREC{t['n']}"] = f"{100 * t['precision']:.1f}"
            T[f"CV_GMAX{t['n']}"] = f"{max(abs(z[1]) for z in t['z']):.1f}"
            T[f"CV_GANC{t['n']}"] = f"{round(t['anc'], -2 if t['anc'] > 1000 else -1):,.0f}"
    if "control" in C:
        T["CV_CTL_MAX"] = f"{C['control']['z_max']:.0f}"; T["CV_CTL_OVER2"] = str(C["control"]["z_over2"]); T["CV_CTL_WITNESS"] = f"{C['control']['auc_z']:.1f}"
        r = max(C["control"]["rows"], key=lambda r: abs(r["z"])); T["CV_CTL_SHIFT"] = f"{abs(r['chain'] - r['exact']) / r['sd']:.2f}"
        e250 = next((e for e in C["exact"] if "n=250" in e["tag"]), None)
        if e250: T["CV_EX250_PREC"] = f"{100 * sorted(math.hypot(r['chain_se'], r['exact_se']) / r['sd'] for r in e250['rows'])[len(e250['rows']) // 2]:.1f}"
    T["CV_WITNESS_MAX"] = f"{max(abs(t['witness']) for t in tests if 'witness' in t and t.get('kind') != 'control'):.1f}"
    ab = C.get("ab", [])
    if ab: T["CV_AB_MAX"] = f"{max(abs(a['z']) for a in ab):.1f}"; T["CV_AB_N"] = str(len(ab))
    big = [s for s in timed if s["n"] >= 16000]
    wmean = lambda key: sum(s[key] / s[key + "_err"] ** 2 for s in big) / sum(1 / s[key + "_err"] ** 2 for s in big)
    tau_rg, tau_single, tau_mode = wmean("rg2"), wmean("single"), wmean("mode")
    T["CV_TAU_RG"] = f"{tau_rg:.0f}"; T["CV_TAU_SINGLE"] = f"{tau_single:.0f}"; T["CV_TAU_MODE"] = f"{tau_mode:.0f}"
    T["CV_SLOW_RATIO"] = f"{tau_single / tau_rg:.1f}"; T["CV_SLOW_NMAX"] = f"{timed[-1]['n']:,}"; T["CV_SLOW_NMIN"] = f"{big[0]['n']:,}"
    sc = C.get("scaling", {})
    if sc:
        T["CV_SLOW_EXP"] = f"{1 + sc['single']['exponent']:.2f}"; T["CV_SLOW_EXP_ERR"] = f"{sc['single']['err']:.2f}"
        T["CV_RG_EXP"] = f"{1 + sc['rg2']['exponent']:.2f}"; T["CV_RG_EXP_ERR"] = f"{sc['rg2']['err']:.2f}"
    # the largest size, counted in the slowest time found (plateau value; and the fitted growth, if any, as the pessimistic end)
    slow_hi = max(tau_mode, sc["single"]["at_1M"] if sc else 0, sc["mode"]["at_1M"] if sc else 0)
    T["NEFF_HERO_SLOW"] = f"{kept_moves / (tau_mode * hero_n):.0f}"; T["NEFF_HERO_SLOW_LO"] = f"{kept_moves / (slow_hi * hero_n):.0f}"
    T["CV_TAU_1M"] = f"{slow_hi:.0f}"
    T["CV_BURN_SLOW"] = f"{burn_n / tau_mode:.0f}"
    if sc:
        x, xe = sc["single"]["exponent"], sc["single"]["err"]
        T["CV_SLOW_TREND"] = (f"Unlike the radius&rsquo;s, this time is still creeping up with size, roughly as n<sup>{1 + x:.2f}</sup> (&plusmn;&nbsp;{xe:.2f}) from 1,000 cells up, though the two largest sizes "
                              f"would also fit a time that has levelled off. If the growth "
                              f"continued to {hero_n:,} cells it would reach about {slow_hi:.0f} moves per cell there, and the largest runs would be about {kept_moves / (slow_hi * hero_n):.0f} "
                              f"relaxation times long instead." if x > 2 * xe else
                              f"A fit from 1,000 cells up gives n<sup>{1 + x:.2f}</sup> (&plusmn;&nbsp;{xe:.2f}): no significant growth beyond proportionality to size.")
    else:
        T["CV_SLOW_TREND"] = ""
    for s in C.get("starts", []):
        T[f"CV_START{s['n']}_MAX"] = f"{s['z_max']:.1f}"; T[f"CV_START{s['n']}_OVER2"] = str(s["z_over2"])
    if C.get("symmetry"):
        z = [r["z"] for s in C["symmetry"] for r in s["rows"]]; ss = _stats(z)
        T["CV_SYM_N"] = str(ss["n"]); T["CV_SYM_MAX"] = f"{ss['max']:.1f}"; T["CV_SYM_OVER2"] = str(ss["over2"])

    # ---- sentences that depend on how the tests came out
    refs = [t for t in trusted if t["group"] in ("enum", "exact", "growth")]
    rz = [z[1] for t in refs for z in t["z"]]; rs = _stats(rz)
    rworst = max(((abs(z[1]), z[0], t["label"].replace("|", ", ")) for t in refs for z in t["z"]))
    T["CV_REF_TESTS"] = str(len(refs)); T["CV_REF_N"] = f"{rs['n']:,}"
    repeat = {}
    for t in refs:
        for nm, z in t["z"]:
            if abs(z) > 2: repeat[nm] = repeat.get(nm, 0) + 1
    most = max(repeat.values()) if repeat else 0
    wit = max(t["witness"] for t in refs if "witness" in t)
    qn = sum(t["q"][0] for t in refs if "q" in t); qo = sum(t["q"][1] for t in refs if "q" in t); qm = max(t["q"][2] for t in refs if "q" in t)
    T["CV_VERDICT"] = (f"Leaving out the starred rows, the three kinds of reference give {len(refs)} tests and {rs['n']:,} comparisons. {rs['over2']} of them ({100 * rs['over2'] / rs['n']:.1f}%) fall outside two standard "
                       f"errors, where chance alone would put about 4.6%, and the largest is {rworst[0]:.1f}. The measurements overlap, so these are not {rs['n']:,} independent "
                       f"trials, but nothing repeats: no measurement is outside two standard errors in more than {most} of the {len(refs)} tests. Against the exact clusters the comparison "
                       f"also covers the shape of each distribution, at its quartiles: {qo} of {qn} outside two standard errors, the largest {qm:.1f}. A last kind of test looks for differences "
                       f"that no single measurement shows. A small neural network is trained to tell the walk&rsquo;s shapes from the reference&rsquo;s, then scored on shapes it has "
                       f"not seen. Against the exact clusters it does no better than chance at any size: "
                       + (f"its best margin over a coin flip is {wit:.1f} standard errors." if wit > 0.05 else "it never beats a coin flip.")
                       + " (The 250-cell cluster row was first run with five chains, and one measurement stood 3.6 standard errors out. Eight more chains, "
                         "with more than twice the data, did not reproduce it; the row shows all thirteen.)")
    other = []
    st16 = next((s for s in C.get("starts", []) if s["n"] == 16000), None) or (C.get("starts") or [None])[-1]
    if st16:
        other.append(f"The third attack starts the walk somewhere absurd: a straight bar of {st16['n']:,} cells, and a solid square. After {st16['burn']:,} moves per cell both have become "
                     f"indistinguishable from the runs that came up the size ladder: of {st16['features']} measurements, {st16['z_over2']} differ by more than two standard errors, the largest by "
                     f"{st16['z_max']:.1f}, and the network "
                     + ("cannot tell them apart." if st16["witness_z"] < 2 else f"tells them apart at {st16['witness_z']:.1f} standard errors, which is unexplained."))
    if C.get("symmetry"):
        z = [r["z"] for s in C["symmetry"] for r in s["rows"]]; ss = _stats(z); sy = C["symmetry"]
        other.append(f"There is also a reference that costs nothing. A polyomino and its mirror image are equally likely, so {sy[0]['features']} of the measurements must average exactly zero "
                     f"at every size. Across {len(sy)} sizes from {sy[0]['n']:,} to {sy[-1]['n']:,} cells, {ss['over2']} of the {ss['n']} averages lie outside two standard errors and the largest is {ss['max']:.1f}.")
        rt = C.get("retest")
        if rt and ss["max"] > 3:
            other.append(f"That largest one, at {rt['n']:,} cells, was put to a second test: two fresh runs at the same size, not used to find it, give the same average as "
                         f"{rt['second']:+.2f}&nbsp;&plusmn;&nbsp;{rt['second_se']:.2f} of its spread, against {rt['first']:+.2f}&nbsp;&plusmn;&nbsp;{rt['first_se']:.2f} the first time. "
                         + ("It did not repeat." if abs(rt["z_second"]) < 2 else "It repeated, which is unexplained."))
    mu = C.get("multi")
    if mu:
        other.append(f"Last, {mu['chains']} runs at {mu['n']} cells were started from one and the same shape with different random numbers, {mu['moves'] / 1e9:.0f} billion moves in all. "
                     f"If some hidden property were stuck, their averages would differ by more than their error bars allow. The scatter between runs, relative to those error "
                     f"bars, is {mu['chi2_mean']:.2f} averaged over the {mu['features']} measurements, where 1 is what correct error bars give.")
    T["CV_OTHER"] = " ".join(other).replace("+-", "+").replace("-0.", "&minus;0.")

    # ---- convergence.md
    md += ["| reference (no Markov chain) | details | comparisons | rms z | largest | beyond 2 s.e. |", "|---|---|---:|---:|---:|---:|"]
    for t in tests:
        s = _stats([z[1] for z in t["z"]])
        md.append(f"| {t['md'][0]} | {t['md'][1]} | {s['n']} | {s['rms']:.2f} | {s['max']:.2f} | {s['over2']} |")
    md += ["", f"z is (chain − reference) / standard error for each observable in `src/features.h`; chain errors are batch means. Over the {len(trusted)} tests "
           f"with a trustworthy reference: {st['n']:,} comparisons, {st['over2']} beyond 2 standard errors ({100 * st['over2'] / st['n']:.1f}%; 4.6% expected if every one "
           f"were independent and exact), {st['over3']} beyond 3, largest {worst[0]:.2f} ({worst[1]}, {worst[2]}). Observables are correlated, so these are not {st['n']:,} independent tests.",
           "", "The growth sampler estimates the number of polyominoes as it goes, and that number is known. Its estimate is the reference's own health check:", "",
           "| cells | population | populations | count estimate vs known | distinct ancestors at n/2 | rms t | largest |", "|---:|---:|---:|---:|---:|---:|---:|"]
    for g in C["growth"] + C.get("growth_small", []):
        md.append(f"| {g['n']:,} | {round(g['population'], -2):,} | {g['replicas']} | {100 * (math.exp(g['count_dev']) - 1):+.1f}% ± {100 * g['count_dev_se']:.1f}% | "
                  f"{sum(g['ancestors_m2']) / len(g['ancestors_m2']):,.0f} | {g['t_rms']:.2f} | {g['t_max']:.2f} |")
    if C.get("multi"):
        mu = C["multi"]
        md += ["", f"Eight independent chains at {mu['n']} cells from one starting shape ({mu['moves']:.2g} moves in all): chi-square of the chain means about their mean, per degree "
               f"of freedom, averaged over {mu['features']} observables: {mu['chi2_mean']:.2f} (largest {mu['chi2_max']:.2f}; 1 if the batch-means errors are right)."]
    if C.get("retest"):
        rt = C["retest"]
        md += ["", f"The largest symmetry violation ({rt['name']} at {rt['n']:,} cells, {rt['first']:+.3f} ± {rt['first_se']:.3f} of its standard deviation) was measured again on two "
               f"fresh chains at that size: {rt['second']:+.3f} ± {rt['second_se']:.3f}."]
    md += ["", "Autocorrelation times in moves per cell (integrated, Sokal window), two independent chains per size:", "",
           "| cells | samples (one per n moves) | Rg² | slowest single observable | slowest combination found (held out) | A-vs-B learned test z |", "|---:|---:|---:|---:|---:|---:|"]
    abz = {a["n"]: a["z"] for a in ab}
    for s in slow:
        md.append(f"| {s['n']:,} | {s['samples']:,} | {s['rg2']:.1f} ± {s['rg2_err']:.1f} | {s['single']:.1f} ± {s['single_err']:.1f} ({s['single_name']}) | {s['mode']:.1f} ± {s['mode_err']:.1f} | {abz.get(s['n'], float('nan')):+.2f} |")
    short = [s for s in slow if s not in timed]
    if short:
        md += ["", "At " + " and ".join(f"{s['n']:,}" for s in short) + " cells the runs are too short to time the slow observables (errors above 50%); that size is left out of the fits and the chart."]
    if sc:
        md += ["", f"Fits from 1,000 cells up, tau/n ∝ n^x: Rg² x = {sc['rg2']['exponent']:+.3f} ± {sc['rg2']['err']:.3f}; slowest single observable x = {sc['single']['exponent']:+.3f} ± {sc['single']['err']:.3f}; "
               f"slowest combination x = {sc['mode']['exponent']:+.3f} ± {sc['mode']['err']:.3f}."]
    open(os.path.join(ROOT, "analysis", "convergence.md"), "w").write("\n".join(md) + "\n")
    return conv, T, md
