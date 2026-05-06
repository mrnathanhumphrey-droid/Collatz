"""
render_K7_full_spectrum_chord.py
================================
Sustained chord rendering of the full top-10 K_7 measured spectrum.
Each measured eigenvalue contributes one sine voice, simultaneously.

Voice mapping (per task brief):
- frequency: f_i = F0 * 2^(arg(λ_i) / 2π), wrapped to [F0, 2·F0)
- decay rate α_i: linear in log|λ_i| from ALPHA_PERRON (|λ|=1) to ALPHA_FAST
  (smallest |λ|), so Perron sustains gently across full duration and
  fastest voice decays in 1-2 seconds
- amplitude weight: w_i = 1 / (1 + |log10|λ_i||) — log-cushioned so the
  Perron drone does not drown subleading shimmer

Carrier: pure sine. Mix is summed, peak-normalized to -1 dBFS, soft 30ms
attack and 50ms release applied to avoid clicks.

Calibration choices documented in stdout for transparency.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.io import wavfile

sys.stdout.reconfigure(encoding="utf-8")

INPUT_CSV = Path(r"C:\Collatz\audio_data\eigenvalue_spectra_consolidated.csv")
OUTPUT_WAV = Path(r"C:\Collatz\audio_data\measured_K7_full_spectrum_chord.wav")
OUTPUT_LOG = Path(r"C:\Collatz\audio_data\measured_K7_full_spectrum_chord.txt")

SR = 44100
DURATION = 30.0
F0 = 207.65  # G#3
TARGET_K = 7

# Decay calibration: α (per second)
ALPHA_PERRON = 0.05   # Perron mode: e^(-1.5) at t=30s ≈ 0.22 (gentle decay)
ALPHA_FAST = 3.0      # Fastest: e^(-3·t); 1/e at t≈0.33s, -40dB by t≈3.1s

# Soft envelope edges to suppress clicks
ATTACK_MS = 30
RELEASE_MS = 50


def voice_freq(arg_rad: float) -> float:
    f = F0 * (2.0 ** (arg_rad / (2.0 * np.pi)))
    while f < F0:
        f *= 2.0
    while f >= 2.0 * F0:
        f /= 2.0
    return f


def main() -> int:
    rows = []
    with open(INPUT_CSV) as f:
        for row in csv.DictReader(f):
            if int(row["k"]) == TARGET_K:
                rows.append({
                    "rank": int(row["rank"]),
                    "mag": float(row["magnitude"]),
                    "arg": float(row["argument_rad"]),
                    "is_complex": int(row["is_complex"]),
                })

    if not rows:
        print(f"No k={TARGET_K} rows found in {INPUT_CSV}")
        return 1

    rows.sort(key=lambda r: -r["mag"])

    mags = np.array([r["mag"] for r in rows])
    log_mags = np.log10(np.maximum(mags, 1e-30))
    log_max = float(log_mags.max())   # 0 for Perron
    log_min = float(log_mags.min())

    def voice_alpha(mag):
        lm = np.log10(max(mag, 1e-30))
        if abs(log_max - log_min) < 1e-12:
            return ALPHA_PERRON
        t = (log_max - lm) / (log_max - log_min)
        return ALPHA_PERRON + t * (ALPHA_FAST - ALPHA_PERRON)

    def voice_weight(mag):
        return 1.0 / (1.0 + abs(np.log10(max(mag, 1e-30))))

    n_samples = int(SR * DURATION)
    t = np.arange(n_samples) / SR
    mix = np.zeros(n_samples, dtype=np.float64)

    log_lines = []

    def log(msg):
        print(msg, flush=True)
        log_lines.append(msg)

    log(f"Rendering K_{TARGET_K} full-spectrum chord")
    log(f"  source: {INPUT_CSV}")
    log(f"  output: {OUTPUT_WAV}")
    log(f"  voices: {len(rows)}")
    log(f"  duration: {DURATION:.1f}s @ {SR} Hz")
    log(f"  fundamental: {F0:.2f} Hz (G#3)")
    log(f"  carrier: sine, mono")
    log("")
    log("Calibration:")
    log(f"  decay: α = {ALPHA_PERRON} (Perron, |λ|=1) → "
        f"{ALPHA_FAST} (smallest |λ|), linear in log10|λ|")
    log(f"  weight: w = 1 / (1 + |log10|λ||)")
    log(f"  log10|λ| span:  [{log_min:.4f}, {log_max:.4f}]")
    log("")
    log(f"  Perron (|λ|=1)            → α={ALPHA_PERRON}/s, e^(-α·30) = "
        f"{np.exp(-ALPHA_PERRON * DURATION):.4f}  (still audible)")
    log(f"  smallest (|λ|={mags.min():.3e}) → α={ALPHA_FAST}/s,  "
        f"e^(-α·1) = {np.exp(-ALPHA_FAST * 1.0):.4f},  "
        f"e^(-α·2) = {np.exp(-ALPHA_FAST * 2.0):.4f}  (1-2 s effective)")
    log("")

    header = (f"  {'rank':>4} {'|λ|':>13} {'arg':>9} {'f_voice':>9}  "
              f"{'α/s':>7} {'weight':>7}  {'cmplx':>5}")
    log(header)
    log("  " + "-" * (len(header) - 2))

    for r in rows:
        f_i = voice_freq(r["arg"])
        a_i = voice_alpha(r["mag"])
        w_i = voice_weight(r["mag"])
        voice = w_i * np.sin(2.0 * np.pi * f_i * t) * np.exp(-a_i * t)
        mix += voice
        log(f"  {r['rank']:>4} {r['mag']:>13.4e} {r['arg']:>+9.4f} "
            f"{f_i:>9.2f}  {a_i:>7.3f} {w_i:>7.4f}  "
            f"{'yes' if r['is_complex'] else 'no':>5}")

    # Soft attack and release
    attack = int(ATTACK_MS * SR / 1000)
    release = int(RELEASE_MS * SR / 1000)
    if attack > 0:
        mix[:attack] *= np.linspace(0.0, 1.0, attack) ** 2
    if release > 0:
        mix[-release:] *= np.linspace(1.0, 0.0, release) ** 2

    # Peak normalize to -1 dBFS
    peak_pre = float(np.max(np.abs(mix)))
    target = 10.0 ** (-1.0 / 20.0)  # -1 dBFS
    mix *= target / peak_pre
    peak_post = float(np.max(np.abs(mix)))

    log("")
    log(f"  peak before normalize: {peak_pre:.4f}")
    log(f"  peak after normalize:  {peak_post:.4f}  (-1 dBFS)")

    # Sanity: per-second RMS to confirm dynamic-range expectations
    log("")
    log("  per-second RMS (envelope verification):")
    log(f"    {'t (s)':>6}  {'RMS':>10}  {'dBFS':>8}")
    for sec in [0, 1, 2, 5, 10, 15, 20, 25, 29]:
        i0 = int(sec * SR)
        i1 = min(int((sec + 1) * SR), n_samples)
        if i0 >= n_samples:
            break
        rms = float(np.sqrt(np.mean(mix[i0:i1] ** 2)))
        db = 20 * np.log10(rms + 1e-30)
        log(f"    {sec:>6}  {rms:>10.4f}  {db:>+8.2f}")

    # Write 16-bit PCM
    OUTPUT_WAV.parent.mkdir(parents=True, exist_ok=True)
    mix_int16 = np.int16(np.clip(mix, -1.0, 1.0) * 32767.0)
    wavfile.write(OUTPUT_WAV, SR, mix_int16)
    log("")
    log(f"  wrote {OUTPUT_WAV}  ({mix_int16.nbytes/1e6:.2f} MB)")

    OUTPUT_LOG.write_text("\n".join(log_lines), encoding="utf-8")
    log(f"  wrote {OUTPUT_LOG}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
