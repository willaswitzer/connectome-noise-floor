"""Arm-2 LIF sweep (PREREG A8.3): case-study set, FlyWire only. Writes spike counts ONLY.

Graphs, in priority order (resumable; finished runs are skipped):
  1. noise floor: the unperturbed v783 graph, re-run 3 times ("base", rerun r0-r2)
  2. each error type at m = 1, injection seeds k = 0-2; then m = 3; then m = 10
     (single-type grid points above a cap are dropped: NT flips at m = 10)
  The combined condition is not run here (A8.3: only if time allows; cut first).
Every graph runs the 43 distinct case-study seed sets of the Amendment-7 FlyWire runs
(data/derived/lif/flywire/seeds.json), 10 trials each (A8.3), with the Amendment-7 model:
Shiu's model.py unchanged except weights (count x A3.4 sign x 0.275 mV), 150 Hz Poisson input
on every seed neuron through one PoissonGroup + Synapses, stimulus weight 68.75 mV, one trial per
subprocess. NT flips are seed-exempt in LIF (A8.2): for each run, the flipped signs of its own seed
neurons and their fragments are restored for that run only.

Output: data/derived/arm2/lif/<condition>/<k or rN>/<run>.parquet (node_id, trial, n_spikes) and
<run>.json. Injected graphs are written to a scratch folder and deleted when the graph is done.

Usage: python src/run_arm2_lif.py [--trials 10] [--max-graphs N]"""
import argparse, json, shutil, subprocess, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "Drosophila_brain_model"))

L = ROOT / "data/derived/lif/flywire"
OUT = ROOT / "data/derived/arm2/lif"
SCR = ROOT / "data/derived/arm2/lif_graphs"
WSYN, RATE = 0.275, 150.0
# exploratory 2026-10-27 (activity-matched control): scale W_syn for network synapses only; the
# stimulus weight w_syn * f_poi stays 68.75 mV. Set through the environment so workers inherit it.
import os
WSYN_SCALE = float(os.environ.get("ARM2_WSYN_SCALE", "1"))
FRAG_ID0 = 2 ** 62
N_BASE = 138639   # v783 scoped-graph nodes (Amendment 7)

ap = argparse.ArgumentParser()
ap.add_argument("--trials", type=int, default=10)
ap.add_argument("--max-graphs", type=int, default=0)
ap.add_argument("--jobs", type=int, default=2, help="trials in parallel (~2.7 GB each on v783)")
ap.add_argument("--_worker", nargs=4, default=None, help=argparse.SUPPRESS)
ap.add_argument("--_write", nargs=2, default=None, help=argparse.SUPPRESS)
ap.add_argument("--ntcm", default="", help="exploratory 2026-10-15: comma list of run_nt_realistic.py "
                "conditions to run instead of the A8.3 grid (k = 0-2)")
ap.add_argument("--glu", default="", help="exploratory 2026-10-16: comma list of run_glu_sign.py "
                "conditions (k = 0-2, condition-major; no seed exemption)")
ap.add_argument("--rungs", default="", help="only runs serving these rungs, e.g. R0,R1 (default: all)")
args = ap.parse_args()


# ----------------------------------------------------------------------------- worker (one trial)
if args._worker is not None:
    gdir, key, tid, fixfile = args._worker
    gdir = Path(gdir)
    from brian2 import Hz, mV, ms, PoissonGroup, Synapses
    import model as shiu

    def poi_group(neu, exc, exc2, params):
        objs = []
        for idx, rate in ((exc, params["r_poi"]), (exc2, params["r_poi2"])):
            if len(idx) == 0:
                continue
            idx = np.asarray(idx)
            P = PoissonGroup(len(idx), rates=rate)
            S = Synapses(P, neu, on_pre="v += w_poi", namespace={"w_poi": params["w_syn"] * params["f_poi"]})
            S.connect(i=np.arange(len(idx)), j=idx)
            neu.rfc[idx] = 0 * ms
            objs += [P, S]
        return objs, neu
    shiu.poi = poi_group

    fix = json.loads(Path(fixfile).read_text()) if fixfile != "-" else []
    fixset = np.asarray(fix, dtype=np.int64)
    _read = pd.read_parquet

    def read_lean(path, *a, **k):
        """model.py uses only these three columns of the connectivity table; reading just them,
        as int32, keeps the largest graphs (spurious m = 10, ~37M edges) within memory.
        Seed-exempt NT flips (A8.2): negate the outgoing weights of the fixed presynaptic indices."""
        df = _read(path, columns=["Presynaptic_Index", "Postsynaptic_Index", "Excitatory x Connectivity"])
        df = df.astype(np.int32)
        if len(fixset):
            m = np.isin(df["Presynaptic_Index"].to_numpy(), fixset)
            df.loc[m, "Excitatory x Connectivity"] *= -1
        return df

    class _PD:  # model.py calls pd.read_parquet / pd.read_csv through its module global
        def __getattr__(self, name):
            return read_lean if name == "read_parquet" else getattr(pd, name)
    shiu.pd = _PD()

    seeds = json.loads((L / "seeds.json").read_text())["runs"][key]
    params = dict(shiu.default_params)
    params.update({"r_poi": RATE * Hz, "n_run": 1, "w_syn": WSYN * WSYN_SCALE * mV,
                   "f_poi": 250 * 0.275 / (WSYN * WSYN_SCALE)})
    shiu.run_exp(exp_name=f"_tmp_{key}_{tid}", neu_exc=seeds, path_res=str(gdir / "tmp"),
                 path_comp=str(gdir / "completeness.csv"), path_con=str(gdir / "connectivity.parquet"),
                 params=params, n_proc=1, force_overwrite=True)
    sys.exit(0)


# ----------------------------------------------------------------------------- driver
import inject_errors as ie


def graph_list():
    if os.environ.get("ARM2_SCALEONLY"):   # exploratory 2026-10-27: unperturbed graph at the scaled W_syn
        return [(f"scaleonly_ws{WSYN_SCALE:g}", "r0", None)]
    if args.glu:
        sfx = "" if WSYN_SCALE == 1 else f"_ws{WSYN_SCALE:g}"
        ks = [int(k) for k in os.environ.get("ARM2_KS", "0,1,2").split(",")]
        return [(c + sfx, f"k{k}", ("glu", c, k)) for c in args.glu.split(",") for k in ks]
    if args.ntcm:
        return [(c, f"k{k}", ("ntcm", c, k)) for k in range(3) for c in args.ntcm.split(",")]
    out = [("base", f"r{i}", None) for i in range(3)]
    for m in (1, 3, 10):
        for t in ie.TYPES:
            rates, _ = ie.condition_rates("single", m, t)
            if rates is None:
                continue
            for k in range(3):
                out.append((f"{t}_m{m:g}", f"k{k}", (rates, k)))
    return out


def write_graph(gdir, spec, base):
    """Shiu-format completeness + connectivity for the base or an injected graph.
    Returns (node ids, flipped node indices, fragment parents)."""
    gdir.mkdir(parents=True, exist_ok=True)
    if spec is None:
        for f in ("completeness.csv", "connectivity.parquet"):
            if not (gdir / f).exists():
                shutil.copy(L / f, gdir / f)
        ids = pd.read_csv(gdir / "completeness.csv", index_col=0).index.to_numpy()
        return ids, np.zeros(0, np.int64), np.zeros(0, np.int64)
    if spec[0] == "glu":   # exploratory 2026-10-16: glutamate sign (run_glu_sign.py draws)
        import run_glu_sign as gs
        _, cond, k = spec
        glu, grp = gs.node_tables(base)
        sign, post_sel, info = gs.draw(base, cond, k, glu, grp)
        pre, post, w = base.pre, base.post, base.w
        esign = sign[pre].astype(int)
        if post_sel is not None:
            m = glu[pre] & post_sel[post]
            esign[m] *= -1
            info["n_edges_flipped"] = int(m.sum())
        ids = base.ids
        con = pd.DataFrame({"Presynaptic_ID": ids[pre], "Postsynaptic_ID": ids[post],
                            "Presynaptic_Index": pre, "Postsynaptic_Index": post,
                            "Connectivity": w.astype(np.int64)})
        con["Excitatory"] = esign
        con["Excitatory x Connectivity"] = con["Excitatory"] * con["Connectivity"]
        con.to_parquet(gdir / "connectivity.parquet")
        pd.DataFrame({"Completed": True}, index=pd.Index(ids, name=None)).to_csv(gdir / "completeness.csv")
        (gdir / "graph.json").write_text(json.dumps({"glu": cond, "k": k, **info}, indent=1))
        return ids, np.zeros(0, np.int64), np.zeros(0, np.int64)
    if spec[0] == "ntcm":   # exploratory 2026-10-15: confusion-matrix flips (run_nt_realistic.py)
        import run_nt_realistic as ntr
        _, cond, k = spec
        flipped, rates, _ = ntr.draw_flips(cond, base, k)
        rates = {"ntcm": cond, "cross_sign_rates": rates}
        pert = ntr.perturbed(base, flipped, k)
    else:
        rates, k = spec
        pert = ie.inject(base, rates, k)
    pre, post, w = pert.edges()
    nfr = len(pert.frag_parent)
    ids = np.concatenate([base.ids, FRAG_ID0 + np.arange(nfr, dtype=np.int64)])
    sign = pert.sign.astype(int)
    con = pd.DataFrame({"Presynaptic_ID": ids[pre], "Postsynaptic_ID": ids[post],
                        "Presynaptic_Index": pre, "Postsynaptic_Index": post,
                        "Connectivity": w.astype(np.int64)})
    con["Excitatory"] = sign[pre]
    con["Excitatory x Connectivity"] = con["Excitatory"] * con["Connectivity"]
    con.to_parquet(gdir / "connectivity.parquet")
    pd.DataFrame({"Completed": True}, index=pd.Index(ids, name=None)).to_csv(gdir / "completeness.csv")
    (gdir / "graph.json").write_text(json.dumps({"rates": rates, "k": k, "edges": int(len(w)),
                                                 "synapses": int(w.sum()), "fragments": int(nfr),
                                                 "flipped": int(len(pert.flipped)),
                                                 "audit": ie.audit(base, pert)}, indent=1))
    return ids, pert.flipped, pert.frag_parent


def write_in_subprocess(cond, rep, gdir):
    """Write the injected graph in a child process, so the driver never holds the ~1.5 GB base
    graph while the two LIF workers run."""
    subprocess.run([sys.executable, __file__, "--_write", cond, rep]
                   + (["--ntcm", args.ntcm] if args.ntcm else [])
                   + (["--glu", args.glu] if args.glu else []), check=True)
    z = np.load(gdir / "flips.npz")
    ids = pd.read_csv(gdir / "completeness.csv", index_col=0).index.to_numpy()
    return ids, z["flipped"], z["frag_parent"]


def main():
    if args._write:
        cond, rep = args._write
        spec = {(c, r): s for c, r, s in graph_list()}[(cond, rep)]
        gdir = SCR / cond / rep
        ids, flipped, fparent = write_graph(gdir, spec, ie.load_flywire_base())
        np.savez(gdir / "flips.npz", flipped=flipped, frag_parent=fparent)
        return
    sj = json.loads((L / "seeds.json").read_text())
    seeds = sj["runs"]
    if args.rungs:
        keep = {x["run"] for x in sj["index"] if x["rung"] in args.rungs.split(",")}
        seeds = {k: v for k, v in seeds.items() if k in keep}
    base = None
    n_graphs = 0
    for cond, rep, spec in graph_list():
        od = OUT / cond / rep
        todo = [k for k in seeds if not (od / f"{k}.parquet").exists()]
        if not todo:
            continue
        if args.max_graphs and n_graphs >= args.max_graphs:
            break
        n_graphs += 1
        t0 = time.time()
        if spec is None:
            gdir = SCR / "base"
            ids, flipped, fparent = write_graph(gdir, None, None)
        else:
            gdir = SCR / cond / rep
            if (gdir / "flips.npz").exists():
                z = np.load(gdir / "flips.npz")
                ids = pd.read_csv(gdir / "completeness.csv", index_col=0).index.to_numpy()
                flipped, fparent = z["flipped"], z["frag_parent"]
            else:
                ids, flipped, fparent = write_in_subprocess(cond, rep, gdir)
        print(f"{cond}/{rep}: graph ready {time.time() - t0:.0f}s", flush=True)
        (gdir / "tmp").mkdir(exist_ok=True)
        od.mkdir(parents=True, exist_ok=True)
        idx_of = pd.Series(np.arange(len(ids)), index=ids)
        jobs = args.jobs   # drops to 1 for the rest of a graph once parallel trials run out of memory
        for key in todo:
            t0 = time.time()
            fixfile = "-"
            n_fix = 0
            if len(flipped):
                sidx = idx_of.reindex(seeds[key]).dropna().astype(np.int64).to_numpy()
                s = np.intersect1d(sidx, flipped)
                frag = N_BASE + np.flatnonzero(np.isin(fparent, s)) if len(s) else np.zeros(0, np.int64)
                fix = np.concatenate([s, frag]).astype(int).tolist()
                n_fix = len(fix)
                if fix:
                    fixfile = str(gdir / f"fix_{key}.json")
                    Path(fixfile).write_text(json.dumps(fix))
            def trial(tid):
                r = subprocess.run([sys.executable, __file__, "--_worker", str(gdir), key, str(tid), fixfile],
                                   stdout=subprocess.DEVNULL, stderr=open(od / "_worker.log", "a"))
                return tid, r.returncode
            from concurrent.futures import ThreadPoolExecutor
            with ThreadPoolExecutor(jobs) as ex:
                failed = [t for t, rc in ex.map(trial, range(args.trials)) if rc != 0]
            for tid in failed:   # e.g. out of memory with 2 in parallel: retry alone
                assert trial(tid)[1] == 0, f"{key} trial {tid} failed twice"
            if len(failed) >= args.trials // 2:
                jobs = 1
            parts = []
            for tid in range(args.trials):
                tmp = gdir / "tmp" / f"_tmp_{key}_{tid}.parquet"
                sp = pd.read_parquet(tmp)
                c = sp.groupby(["flywire_id", "trial"]).size().rename("n_spikes").reset_index()
                c["trial"] = tid
                parts.append(c)
                tmp.unlink()
            res = pd.concat(parts, ignore_index=True).rename(columns={"flywire_id": "node_id"})
            res.to_parquet(od / f"{key}.parquet")
            meta = {"condition": cond, "rep": rep, "run": key, "trials": args.trials, "wsyn_mV": WSYN * WSYN_SCALE,
                    "rate_hz": RATE, "jobs": jobs, "retried": len(failed), "n_seed_neurons": len(seeds[key]), "n_seed_exempt": n_fix,
                    "n_active_neurons": int(res["node_id"].nunique()),
                    "elapsed_s": round(time.time() - t0, 1), "date": time.strftime("%Y-%m-%d %H:%M")}
            (od / f"{key}.json").write_text(json.dumps(meta, indent=1))
            print(f"  {cond}/{rep} {key}: {meta['n_active_neurons']} active, {meta['elapsed_s']:.0f}s", flush=True)
        if spec is not None:
            shutil.rmtree(gdir, ignore_errors=True)


if __name__ == "__main__":
    main()
