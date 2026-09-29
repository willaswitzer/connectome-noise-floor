"""Week-5 checkpoint rerun on FlyWire v783 (A3.4 signs), with timing for the A6.4 LIF-scope decision.

The Sep 21 checkpoint (PREREG, Week 5 section) ran Shiu's example notebook, which uses the v630
files (2023_03_23_*_630_final) and Shiu's `Excitatory` column. All Arm-1 LIF runs use v783 and
A3.4 signs, so this reruns the same experiment on v783:
  stimulus  Shiu's sugar GRN list (example.ipynb cell 3); the 20 IDs present in v783
            (720575940620900446 is not a v783 root ID). Poisson input 100 Hz.
  readout   MN9, root 720575940660219265 (present in v783). Target 67.0 Hz.
  model     Shiu model.py unchanged (w_syn 0.275 mV, 30 trials x 1 s); only the synapse sign
            differs by --signs:
              shiu          Shiu's v783 `Excitatory` column (Connectivity_783.parquet)
              a34           A3.4 sign per presynaptic neuron (signs_flywire.parquet `sign`);
                            neurons with no sign (tier 4 or absent) pass nothing on (w = 0)
              a34_predonly  A3.4 sensitivity signs (`sign_predonly`)
              shiu630       control: Shiu's v630 files and all 21 IDs, i.e. the Sep 21 checkpoint
                            re-run in this environment, to separate environment from version
Writes results/checkpoints/shiu_783/<signs>/ : spikes parquet, rates.csv, summary.json.

Trials run serially in subprocesses of --chunk trials each (brian2 leaks memory across trials, and
2 parallel workers exceed 8 GB). Usage (from repo root):
    python src/shiu_checkpoint_783.py --signs a34 [--n-run 30] [--chunk 10]"""
import argparse, json, platform, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SHIU = ROOT / "Drosophila_brain_model"
sys.path.insert(0, str(SHIU))
from brian2 import Hz, ms  # noqa: E402
import model as shiu  # noqa: E402

SUGAR = [720575940624963786, 720575940630233916, 720575940637568838, 720575940638202345,
         720575940617000768, 720575940630797113, 720575940632889389, 720575940621754367,
         720575940621502051, 720575940640649691, 720575940639332736, 720575940616885538,
         720575940639198653, 720575940620900446, 720575940617937543, 720575940632425919,
         720575940633143833, 720575940612670570, 720575940628853239, 720575940629176663,
         720575940611875570]
MN9 = 720575940660219265
TARGET_HZ = 67.0

ap = argparse.ArgumentParser()
ap.add_argument("--signs", choices=["shiu", "a34", "a34_predonly", "shiu630"], required=True)
ap.add_argument("--n-run", type=int, default=30)
ap.add_argument("--chunk", type=int, default=10,
                help="trials per subprocess; brian2 networks leak memory across trials in one process")
ap.add_argument("--_worker", type=int, default=None, help=argparse.SUPPRESS)
ap.add_argument("--_con", default=None, help=argparse.SUPPRESS)
args = ap.parse_args()

out = ROOT / "results/checkpoints/shiu_783" / args.signs
out.mkdir(parents=True, exist_ok=True)
V630 = args.signs == "shiu630"
comp_path = SHIU / ("2023_03_23_completeness_630_final.csv" if V630 else "Completeness_783.csv")
if args._worker is not None:  # one chunk of trials in a fresh process
    params = dict(shiu.default_params)
    params["r_poi"] = 100 * Hz
    params["n_run"] = args.n_run
    comp = pd.read_csv(comp_path, index_col=0)
    shiu.run_exp(exp_name=f"chunk{args._worker}", neu_exc=[s for s in SUGAR if s in comp.index],
                 path_res=str(out), path_comp=str(comp_path), path_con=args._con, params=params,
                 n_proc=1, force_overwrite=True)
    sys.exit(0)
CON_PATH = SHIU / ("2023_03_23_connectivity_630_final.parquet" if V630 else "Connectivity_783.parquet")
import pyarrow.parquet as pq
SIGNED = args.signs not in ("shiu", "shiu630")
con = pd.read_parquet(CON_PATH) if SIGNED else None  # loaded only when signs are replaced (memory)
comp = pd.read_csv(comp_path, index_col=0)
stim = [s for s in SUGAR if s in comp.index]
assert MN9 in comp.index and len(stim) == (21 if V630 else 20), (len(stim), MN9 in comp.index)

info = {"signs": args.signs, "n_synapse_rows": pq.ParquetFile(CON_PATH).metadata.num_rows}
if SIGNED:
    col = "sign" if args.signs == "a34" else "sign_predonly"
    s = pd.read_parquet(ROOT / "data/derived/signs_flywire.parquet", columns=["root_id", col]).set_index("root_id")[col]
    new = con["Presynaptic_ID"].map(s).fillna(0).astype(int)
    shiu_sign = np.sign(con["Excitatory"]).astype(int)
    w = con["Connectivity"]
    info.update({
        "synapses_same_sign_as_shiu": float(w[new == shiu_sign].sum() / w.sum()),
        "synapses_flipped_vs_shiu": float(w[(new != shiu_sign) & (new != 0)].sum() / w.sum()),
        "synapses_zeroed_no_sign": float(w[new == 0].sum() / w.sum()),
        "presyn_neurons_flipped": int(con.loc[(new != shiu_sign) & (new != 0), "Presynaptic_ID"].nunique()),
        "presyn_neurons_no_sign": int(con.loc[new == 0, "Presynaptic_ID"].nunique()),
    })
    con["Excitatory"] = new
    con["Excitatory x Connectivity"] = new * w
    con_path = out / "connectivity_783_signed.parquet"
    con.to_parquet(con_path)
else:
    con_path = CON_PATH
del con

params = dict(shiu.default_params)
params["r_poi"] = 100 * Hz
params["n_run"] = args.n_run
t0 = time.time()
import subprocess
chunks, done = [], 0
while done < args.n_run:
    n = min(args.chunk, args.n_run - done)
    subprocess.run([sys.executable, __file__, "--signs", args.signs, "--n-run", str(n),
                    "--_worker", str(len(chunks)), "--_con", str(con_path)], check=True)
    c = pd.read_parquet(out / f"chunk{len(chunks)}.parquet")
    c["trial"] += done
    chunks.append(c)
    (out / f"chunk{len(chunks) - 1}.parquet").unlink()
    done += n
sp_all = pd.concat(chunks, ignore_index=True)
sp_all["exp_name"] = f"sugar_783_{args.signs}_100Hz"
sp_all.to_parquet(out / f"sugar_783_{args.signs}_100Hz.parquet", compression="brotli")
elapsed = time.time() - t0
if con_path != CON_PATH:
    Path(con_path).unlink()

sp = pd.read_parquet(out / f"sugar_783_{args.signs}_100Hz.parquet")
T = float(params["t_run"] / ms) / 1000
per_trial = sp.groupby(["flywire_id", "trial"]).size().unstack(fill_value=0).reindex(columns=range(args.n_run), fill_value=0) / T
rates = pd.DataFrame({"rate_hz": per_trial.mean(axis=1), "sd_hz": per_trial.std(axis=1)}).sort_values("rate_hz", ascending=False)
rates.to_csv(out / "rates.csv")
mn9 = per_trial.loc[MN9] if MN9 in per_trial.index else pd.Series(0.0, index=range(args.n_run))
info.update({
    "stimulus_ids_used": len(stim), "r_poi_hz": 100, "n_run": args.n_run, "t_run_s": T, "w_syn_mV": 0.275,
    "mn9_rate_hz": float(mn9.mean()), "mn9_sem_hz": float(mn9.std() / np.sqrt(args.n_run)),
    "target_hz": TARGET_HZ, "n_active_neurons": int((rates["rate_hz"] > 0).sum()),
    "elapsed_s": elapsed, "n_proc": 1, "trials_per_process": args.chunk, "machine": platform.platform(), "python": platform.python_version(),
    "date": time.strftime("%Y-%m-%d %H:%M"),
})
(out / "summary.json").write_text(json.dumps(info, indent=1))
print(json.dumps(info, indent=1))
