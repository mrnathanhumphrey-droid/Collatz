"""
render_K7_impulse.py
====================
Construct impulse response h(t) of K_7 from its top-20 eigenvalues
(loaded from eigenvalue_spectra_consolidated.csv) and render as audio.

Two outputs:
  impulse_response_K7.wav            -- AM-modulated single carrier (G#3)
  impulse_response_K7_multivoice.wav -- per-eigenvalue carrier voices

Method per brief:
  h_i(t) = |lambda_i|^t * cos(arg(lambda_i) * t)   (logical step t)
  h(t)   = sum_i h_i(t)
  1 logical step = 100 ms.
"""
from __future__ import annotations

import csv
import os
import sys
import numpy as np
import scipy.io.wavfile as wavfile

sys.stdout.reconfigure(encoding="utf-8")

CSV_IN = r"C:\Collatz\audio_data\eigenvalue_spectra_consolidated.csv"
OUT_DIR = r"C:\Collatz\audio_data"

SAMPLE_RATE = 44100
LOGICAL_STEP_MS = 100
SAMPLES_PER_STEP = int(SAMPLE_RATE * LOGICAL_STEP_MS / 1000)  # 4410
N_STEPS = 80                          # 8 seconds total
CARRIER_FREQ = 207.65                 # G#3
TARGET_PEAK_DBFS = -1.0               # -1 dBFS normalization

target_amp_lin = 10 ** (TARGET_PEAK_DBFS / 20.0)


# -------- load top-20 K_7 eigenvalues --------

def load_k7_eigs():
    eigs = []
    with open(CSV_IN) as f:
        for row in csv.DictReader(f):
            if int(row["k"]) != 7:
                continue
            re = float(row["real_part"])
            im = float(row["imag_part"])
            eigs.append(complex(re, im))
    eigs.sort(key=lambda x: -abs(x))
    return eigs[:20]


# -------- impulse response h(t) at logical-step rate --------

def build_h_logical(eigs, n_steps):
    t = np.arange(n_steps, dtype=np.float64)
    h = np.zeros(n_steps)
    contributions = []
    for lam in eigs:
        rho = abs(lam)
        theta = float(np.angle(lam))
        # rho^t * cos(theta * t). For rho == 0 we keep the t=0 contribution
        # (cos(0)=1) and treat 0^positive as 0 cleanly.
        h_i = (rho ** t) * np.cos(theta * t)
        h += h_i
        contributions.append((lam, rho, theta, h_i))
    return t, h, contributions


# -------- interpolation logical -> audio --------

def interp_to_audio(t_logical, h_logical, n_steps):
    n_audio = n_steps * SAMPLES_PER_STEP
    t_audio = np.linspace(0, n_steps - 1, n_audio, endpoint=False)
    h_audio = np.interp(t_audio, t_logical, h_logical)
    return t_audio, h_audio, n_audio


# -------- main --------

def normalize_to_dbfs(x, target_lin=target_amp_lin):
    peak = float(np.max(np.abs(x)))
    if peak == 0:
        return x, 0.0
    scale = target_lin / peak
    return x * scale, peak


def to_int16(x):
    return np.clip(x * 32767.0, -32768, 32767).astype(np.int16)


def main():
    print("=" * 78)
    print("Impulse response of K_7 -> audio")
    print("=" * 78)
    print()

    eigs = load_k7_eigs()
    print(f"Loaded {len(eigs)} K_7 eigenvalues from consolidated CSV.")
    print()
    print(f"{'rank':>4}  {'|lam|':>14}  {'arg (rad)':>12}  "
          f"{'arg (deg)':>10}  Re      Im")
    excluded = []
    for i, lam in enumerate(eigs, 1):
        rho = abs(lam)
        theta = float(np.angle(lam))
        if not np.isfinite(rho) or not np.isfinite(theta):
            print(f"  [skip] rank {i}: non-finite eigenvalue")
            excluded.append((i, lam, "non-finite"))
            continue
        print(f"  {i:>3}  {rho:>14.6e}  {theta:>+12.6f}  {np.degrees(theta):>+10.3f}  "
              f"{lam.real:>+11.4e}  {lam.imag:>+11.4e}")
    print()

    eigs_use = [e for e in eigs if np.isfinite(abs(e))
                and np.isfinite(float(np.angle(e)))]

    # ---- h(t) at logical step ----
    t_log, h_log, contribs = build_h_logical(eigs_use, N_STEPS)
    print(f"h(t) at logical steps t=0..{N_STEPS-1}:")
    print(f"  h(0) = {h_log[0]:+.6f}  (sum of |lam_i|^0 * cos(0) = "
          f"{len(eigs_use)})")
    print(f"  h(1) = {h_log[1]:+.6f}")
    print(f"  h(2) = {h_log[2]:+.6f}")
    print(f"  h(5) = {h_log[5]:+.6f}")
    print(f"  h(10)= {h_log[10]:+.6f}  (Perron contribution dominates)")
    print(f"  h({N_STEPS-1}) = {h_log[-1]:+.6f}")
    print()

    # ---- interpolate to audio rate ----
    t_aud, h_aud, n_audio = interp_to_audio(t_log, h_log, N_STEPS)
    duration_s = n_audio / SAMPLE_RATE
    print(f"Audio rendering parameters:")
    print(f"  sample rate: {SAMPLE_RATE} Hz")
    print(f"  logical step: {LOGICAL_STEP_MS} ms = {SAMPLES_PER_STEP} samples")
    print(f"  total audio samples: {n_audio}, duration {duration_s:.2f}s")
    print()

    # ===== AM-MODULATED SINGLE CARRIER =====
    t_seconds = np.arange(n_audio) / SAMPLE_RATE
    carrier = np.cos(2 * np.pi * CARRIER_FREQ * t_seconds)
    am_signal = h_aud * carrier
    am_signal_norm, am_peak = normalize_to_dbfs(am_signal)
    am_int16 = to_int16(am_signal_norm)

    am_path = os.path.join(OUT_DIR, "impulse_response_K7.wav")
    wavfile.write(am_path, SAMPLE_RATE, am_int16)
    print(f"AM-modulated rendering:")
    print(f"  carrier: {CARRIER_FREQ} Hz (G#3)")
    print(f"  raw peak before normalization: {am_peak:.4f}")
    print(f"  normalized to -1 dBFS, clipped to int16")
    print(f"  saved {am_path}")
    print()

    # ===== MULTIVOICE: each eigenvalue gets own carrier =====
    # Map arg in [-pi, pi] to a musical frequency range. Use 100..1000 Hz.
    multi = np.zeros(n_audio)
    voice_log = []
    for lam in eigs_use:
        rho = abs(lam)
        theta = float(np.angle(lam))
        # frequency mapping: theta in [-pi, pi] -> [100, 1000] Hz (linear)
        freq = 100.0 + (theta + np.pi) / (2 * np.pi) * 900.0
        # decay envelope at logical rate, interpolated to audio
        env_log = rho ** t_log
        env_aud = np.interp(t_aud, t_log, env_log)
        voice = env_aud * np.cos(2 * np.pi * freq * t_seconds)
        multi += voice
        voice_log.append((lam, freq, float(np.max(np.abs(voice)))))

    multi_norm, multi_peak = normalize_to_dbfs(multi)
    multi_int16 = to_int16(multi_norm)
    multi_path = os.path.join(OUT_DIR, "impulse_response_K7_multivoice.wav")
    wavfile.write(multi_path, SAMPLE_RATE, multi_int16)

    print(f"Multivoice rendering:")
    print(f"  20 voices, frequency from arg in [100, 1000] Hz")
    print(f"  raw peak before normalization: {multi_peak:.4f}")
    print(f"  normalized to -1 dBFS, clipped to int16")
    print(f"  saved {multi_path}")
    print()
    print(f"  Per-voice frequencies and peak contributions:")
    for i, (lam, freq, vp) in enumerate(voice_log, 1):
        print(f"    rank {i:>2}: freq {freq:>7.1f} Hz, |lam| = {abs(lam):.4e}, "
              f"voice peak {vp:.4f}")
    print()

    # ---- dynamic range / report ----
    print("=" * 78)
    print("Report")
    print("=" * 78)
    print(f"AM file:")
    print(f"  total duration: {duration_s:.2f}s ({n_audio} samples)")
    print(f"  peak before normalization: {am_peak:.4f}")
    print(f"  dynamic range note: h(0) = {h_log[0]:.2f} dominates the impulse")
    print(f"    (sum of 20 |lam|^0 cos(0)); h(t>=1) drops by factor "
          f"{h_log[0]/max(abs(h_log[1]), 1e-30):.2f}.")
    print(f"    The Perron mode (lam=1) provides a constant +1 floor; the")
    print(f"    impulse spike at t=0 is followed by a sustained DC tone in the")
    print(f"    envelope, so the AM signal is a brief transient-into-steady carrier.")
    print()
    print(f"Multivoice file:")
    print(f"  total duration: {duration_s:.2f}s")
    print(f"  peak before normalization: {multi_peak:.4f}")
    print(f"  dynamic range note: 19 of 20 voices have |lam| <= 2e-3, so")
    print(f"    their envelopes decay to ~10^-12 within 4 logical steps")
    print(f"    (~400 ms). Only the Perron voice (|lam|=1, freq mapped from")
    print(f"    arg=0 -> 550 Hz) sustains. Audible result is a dense brief click")
    print(f"    at t=0 plus a sustained 550 Hz tone.")
    print()
    if excluded:
        print(f"  excluded eigenvalues (non-finite): {excluded}")
    else:
        print("  no eigenvalues excluded.")


if __name__ == "__main__":
    main()
