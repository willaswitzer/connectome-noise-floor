"""Run the Shiu LIF model on the case-study set (Amendment 7). Resumable: finished runs are skipped.

Inputs from src/lif_prep.py (data/derived/lif/<dataset>/). For every distinct seed set in
seeds.json: 30 trials x 1 s, Poisson input at 150 Hz on every seed neuron, Shiu's model.py
unchanged except:
  - weights: synapse count x A3.4 sign x Wsyn
  - Wsyn: --wsyn primary = structural rescaling (A1.4, results/wsyn.csv)
          --wsyn fixed   = 0.275 mV (A1.4 sensitivity)
  - Poisson input weight is held at Shiu's 0.275 mV x 250 = 68.75 mV in every condition
    (f_poi = 250 x 0.275 / Wsyn). Wsyn rescales synapses between neurons (A1.4); it must not
    change how hard the stimulus drives the seed neurons.
Each trial runs in its own subprocess (--chunk 1): in one process brian2 slows down about 10x and
leaks memory after the first trial with large seed sets. --jobs runs subprocesses in parallel.

Output per run: data/derived/lif/<dataset>/<wsyn>/<run>.parquet with columns node_id, trial,
n_spikes (neurons with >= 1 spike in that trial only), and <run>.json with timing.

Usage: python src/run_lif.py --dataset flywire [--wsyn primary] [--n-run 30] [--chunk 10]"""
import argparse, json, subprocess, sys, time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Drosophila_brain_model"))

ap = argparse.ArgumentParser()
ap.add_argument("--dataset", required=True, choices=["flywire", "male_cns", "hemibrain"])
ap.add_argument("--wsyn", default="primary", choices=["primary", "fixed"])
ap.add_argument("--rate", type=float, default=150.0)
ap.add_argument("--n-run", type=int, default=30)
ap.add_argument("--chunk", type=int, default=1,
                help="trials per subprocess; brian2 slows down and leaks memory after the first trial in a process")
ap.add_argument("--jobs", type=int, default=1, help="subprocesses run at once (each needs ~2.6 GB on FlyWire)")
ap.add_argument("--only", default="", help="comma-separated run keys")
ap.add_argument("--poisson", default="group", choices=["group", "shiu"],
                help="group: one PoissonGroup + Synapses for all seeds (default); shiu: one PoissonInput per seed")
ap.add_argument("--out-suffix", default="", help=argparse.SUPPRESS)
ap.add_argument("--_worker", default=None, help=argparse.SUPPRESS)
args = ap.parse_args()

D = ROOT / f"data/derived/lif/{args.dataset}"
OUT = D / (args.wsyn + args.out_suffix)
OUT.mkdir(parents=True, exist_ok=True)
seeds = json.loads((D / "seeds.json").read_text())
wsyn_mv = 0.275 if args.wsyn == "fixed" else float(
    pd.read_csv(ROOT / "results/wsyn.csv").set_index("dataset").loc[args.dataset, "wsyn_mV"])

if args._worker is not None:  # one chunk of trials for one run, in a fresh process
    from brian2 import Hz, mV
    import model as shiu
    key, chunk_id = args._worker.split(":")
    if args.poisson == "group":
        import numpy as np
        from brian2 import PoissonGroup, Synapses, ms

        def poi_group(neu, exc, exc2, params):
            """Same stimulus as Shiu's poi(): each seed neuron gets an independent Poisson train
            whose events add w_syn * f_poi to v, and seeds have no refractory period. One
            PoissonGroup + Synapses replaces one PoissonInput per neuron, whose per-subgroup code
            brian2 compiles separately (minutes per new seed set with hundreds of seeds)."""
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
    params = dict(shiu.default_params)
    params.update({"r_poi": args.rate * Hz, "n_run": args.n_run, "w_syn": wsyn_mv * mV,
                   "f_poi": 250 * 0.275 / wsyn_mv})
    shiu.run_exp(exp_name=f"_tmp_{key}_{chunk_id}", neu_exc=seeds["runs"][key], path_res=str(OUT),
                 path_comp=str(D / "completeness.csv"), path_con=str(D / "connectivity.parquet"),
                 params=params, n_proc=1, force_overwrite=True)
    sys.exit(0)

keys = [k for k in seeds["runs"] if not args.only or k in args.only.split(",")]
todo = [k for k in keys if not (OUT / f"{k}.parquet").exists()]
print(f"{args.dataset}/{args.wsyn}: Wsyn {wsyn_mv:.4f} mV, {len(keys)} runs, {len(todo)} to do", flush=True)
def run_chunk(key, cid, n):
    subprocess.run([sys.executable, __file__, "--dataset", args.dataset, "--wsyn", args.wsyn,
                    "--rate", str(args.rate), "--n-run", str(n), "--poisson", args.poisson,
                    "--out-suffix", args.out_suffix, "--_worker", f"{key}:{cid}"],
                   check=True, stdout=subprocess.DEVNULL, stderr=open(OUT / "_worker.log", "a"))
    return cid


from concurrent.futures import ThreadPoolExecutor
for i, key in enumerate(todo):
    t0 = time.time()
    sizes = [min(args.chunk, args.n_run - s0) for s0 in range(0, args.n_run, args.chunk)]
    with ThreadPoolExecutor(args.jobs) as ex:
        list(ex.map(lambda c: run_chunk(key, c, sizes[c]), range(len(sizes))))
    parts, done = [], 0
    for cid, n in enumerate(sizes):
        tmp = OUT / f"_tmp_{key}_{cid}.parquet"
        sp = pd.read_parquet(tmp)
        c = sp.groupby(["flywire_id", "trial"]).size().rename("n_spikes").reset_index()
        c["trial"] += done
        parts.append(c)
        tmp.unlink()
        done += n
    res = pd.concat(parts, ignore_index=True).rename(columns={"flywire_id": "node_id"})
    res.to_parquet(OUT / f"{key}.parquet")
    meta = {"run": key, "dataset": args.dataset, "wsyn": args.wsyn, "wsyn_mV": wsyn_mv, "rate_hz": args.rate,
            "n_run": args.n_run, "t_run_s": 1.0, "trials_per_process": args.chunk, "jobs": args.jobs, "poisson": args.poisson, "n_seed_neurons": len(seeds["runs"][key]),
            "n_active_neurons": int(res["node_id"].nunique()), "elapsed_s": time.time() - t0,
            "date": time.strftime("%Y-%m-%d %H:%M")}
    (OUT / f"{key}.json").write_text(json.dumps(meta, indent=1))
    print(f"[{i + 1}/{len(todo)}] {key}: {meta['n_seed_neurons']} seeds, {meta['n_active_neurons']} active, "
          f"{meta['elapsed_s']:.0f}s", flush=True)
