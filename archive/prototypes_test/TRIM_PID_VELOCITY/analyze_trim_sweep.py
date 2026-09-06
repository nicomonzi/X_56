#!/usr/bin/env python3
"""Analyse the MBDyn trim sweep without modifying the original results."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, sosfiltfilt, welch
from scipy.stats import linregress


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "trim_sweep_results"
OUT = RESULTS / "analysis"
BFF_SPEED = 42.3
LBF_TO_N = 4.4482216152605
LBFIN_TO_NM = 0.1129848290276167


def lowpass(values: np.ndarray, cutoff_hz: float, fs_hz: float) -> np.ndarray:
    sos = butter(4, cutoff_hz / (0.5 * fs_hz), btype="low", output="sos")
    return sosfiltfilt(sos, values)


def rms_envelope(values: np.ndarray, fs_hz: float, window_s: float = 1.0) -> np.ndarray:
    samples = max(3, int(round(window_s * fs_hz)))
    return np.sqrt(uniform_filter1d(values * values, size=samples, mode="nearest"))


def smooth_xy(x: np.ndarray, y: np.ndarray, xmax: float | None = None):
    mask = np.isfinite(x) & np.isfinite(y)
    if xmax is not None:
        mask &= x <= xmax
    xx = x[mask]
    yy = y[mask]
    grid = np.linspace(xx.min(), xx.max(), 500)
    return grid, PchipInterpolator(xx, yy)(grid)


def savefig(fig: plt.Figure, name: str) -> None:
    fig.savefig(OUT / name, dpi=220, bbox_inches="tight")
    plt.close(fig)


def load_cases():
    cases = {}
    for history in sorted(RESULTS.glob("V_*_mps/history.csv")):
        speed = float(history.parent.name.removeprefix("V_").removesuffix("_mps"))
        data = np.genfromtxt(history, delimiter=",", names=True)
        cases[speed] = data
    if not cases:
        raise RuntimeError(f"No history.csv files found below {RESULTS}")
    return cases


def load_six_component_reaction(case_directory: Path, time: np.ndarray):
    """Load local Fx,Fy,Fz,Mx,My,Mz of trim joint 23 from the MBDyn text output."""
    samples = []
    with (case_directory / "result.jnt").open(errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) < 7:
                continue
            try:
                if int(fields[0]) != 23:
                    continue
                samples.append([float(value.replace("D", "E")) for value in fields[1:7]])
            except ValueError:
                continue
    reaction = np.asarray(samples, dtype=float)
    count = min(len(time), len(reaction))
    reaction = reaction[-count:]
    reaction[:, :3] *= LBF_TO_N
    reaction[:, 3:] *= LBFIN_TO_NM
    return time[-count:], reaction


def analyse_case(speed: float, data: np.ndarray):
    t = data["time_s"]
    my = data["My_Nm"]
    dt = float(np.median(np.diff(t)))
    fs = 1.0 / dt

    my_slow = lowpass(my, 2.0, fs)
    baseline = lowpass(my_slow, 0.15, fs)
    oscillatory = my_slow - baseline
    envelope = rms_envelope(oscillatory, fs)

    final = t >= t[-1] - 1.0
    interior_final = (t >= t[-1] - 2.0) & (t <= t[-1] - 0.5)
    spectral = t >= 5.0

    freq, psd = welch(my[spectral] - np.mean(my[spectral]), fs=fs, nperseg=256)
    positive = freq > 0.05
    dominant_hz = float(freq[positive][np.argmax(psd[positive])])
    total_power = np.trapz(psd[positive], freq[positive])
    high = freq >= 10.0
    high_fraction = float(np.trapz(psd[high], freq[high]) / total_power)

    slow_freq, slow_psd = welch(
        oscillatory[spectral], fs=fs, nperseg=min(256, np.count_nonzero(spectral))
    )
    slow_band = (slow_freq >= 0.15) & (slow_freq <= 2.0)
    slow_dominant_hz = float(slow_freq[slow_band][np.argmax(slow_psd[slow_band])])

    fit = (t >= 5.0) & (t <= 13.0)
    floor = max(np.percentile(envelope[fit], 5) * 0.25, 1e-8)
    usable = fit & (envelope > floor)
    regression = linregress(t[usable], np.log(envelope[usable]))

    if speed >= BFF_SPEED:
        status = "beyond_BFF_reference"
    elif speed >= 39.0:
        status = "near_BFF_nonstationary"
    elif abs(np.mean(data["Fz_N"][final])) > 10.0:
        status = "incomplete_settling"
    else:
        status = "usable_trim"

    row = {
        "speed_mps": speed,
        "theta_deg": float(np.mean(data["theta_deg"][final])),
        "elevator_deg": float(np.mean(data["delta_elevator_deg"][final])),
        "Fz_mean_N": float(np.mean(data["Fz_N"][final])),
        "Fz_std_N": float(np.std(data["Fz_N"][final])),
        "My_mean_Nm": float(np.mean(my[final])),
        "My_raw_std_Nm": float(np.std(my[final])),
        "My_slow_std_Nm": float(np.std(my_slow[interior_final])),
        "My_ripple_rms_Nm": float(
            np.sqrt(np.mean((my[interior_final] - my_slow[interior_final]) ** 2))
        ),
        "dominant_raw_frequency_Hz": dominant_hz,
        "power_fraction_above_10Hz": high_fraction,
        "dominant_slow_frequency_Hz": slow_dominant_hz,
        "slow_envelope_rate_1_s": float(regression.slope),
        "slow_envelope_fit_R2": float(regression.rvalue**2),
        "status": status,
    }
    auxiliary = {
        "t": t,
        "my": my,
        "my_slow": my_slow,
        "oscillatory": oscillatory,
        "envelope": envelope,
        "freq": freq,
        "psd": psd,
        "fs": fs,
    }
    return row, auxiliary


def add_speed_context(ax):
    ax.grid(True, alpha=0.25)
    ax.set_xlim(19.5, 42.2)


def plot_trim_solution(rows):
    v = np.array([r["speed_mps"] for r in rows])
    theta = np.array([r["theta_deg"] for r in rows])
    elevator = np.array([r["elevator_deg"] for r in rows])

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 7.4), sharex=True)
    for ax, y, ylabel, color in [
        (axes[0], theta, r"Pitch trim $\theta$ [deg]", "#2166ac"),
        (axes[1], elevator, r"Elevator trim $\delta_e$ [deg]", "#1b7837"),
    ]:
        grid, interp = smooth_xy(v, y, xmax=42.0)
        ax.plot(grid, interp, color=color, lw=2.2)
        pre = v <= 42.0
        ax.scatter(v[pre], y[pre], s=28, color=color, edgecolor="white", zorder=3)
        ax.set_ylabel(ylabel)
        add_speed_context(ax)
    axes[1].set_xlabel("Airspeed [m/s]")
    axes[0].set_ylabel(r"$\theta$ [deg]")
    axes[1].set_ylabel(r"$\delta$ control surfaces [deg]")
    fig.suptitle("Trim variables")
    fig.tight_layout()
    savefig(fig, "01_trim_solution_interpolated.png")


def plot_quality(rows):
    v = np.array([r["speed_mps"] for r in rows])
    fz = np.array([r["Fz_mean_N"] for r in rows])
    slow_std = np.array([r["My_slow_std_Nm"] for r in rows])
    ripple = np.array([r["My_ripple_rms_Nm"] for r in rows])
    rates = np.array([r["slow_envelope_rate_1_s"] for r in rows])
    fit_quality = np.array([r["slow_envelope_fit_R2"] for r in rows])

    fig, axes = plt.subplots(3, 1, figsize=(9.2, 9.4), sharex=True)
    axes[0].axhline(0, color="0.35", lw=0.9)
    axes[0].plot(v, fz, "o", color="#2166ac", ms=4)
    grid, interp = smooth_xy(v, fz, xmax=42.0)
    axes[0].plot(grid, interp, color="#2166ac", lw=1.7)
    axes[0].set_ylabel(r"Mean $F_z$ [N]")

    pre = v <= 42.0
    for values, marker, color, label in [
        (slow_std, "o", "#1b7837", "Slow component (<2 Hz)"),
        (ripple, "s", "#762a83", "Raw high-frequency ripple"),
    ]:
        grid, interp = smooth_xy(v, np.maximum(values, 1e-3), xmax=42.0)
        axes[1].semilogy(grid, interp, color=color, lw=1.7, label=label)
        axes[1].semilogy(v[pre], np.maximum(values[pre], 1e-3), marker, ls="none",
                        color=color, ms=5)
    axes[1].set_ylabel(r"Moment variability [N m]")
    axes[1].legend(fontsize=8)

    good = fit_quality >= 0.25
    axes[2].axhline(0, color="0.35", lw=0.9)
    axes[2].scatter(v[good], rates[good], c=rates[good], cmap="coolwarm", vmin=-0.3,
                    vmax=0.3, edgecolor="black", linewidth=0.35, label=r"Envelope fit $R^2\geq0.25$")
    axes[2].scatter(v[~good], rates[~good], facecolor="none", edgecolor="0.55",
                    label="Unreliable exponential fit")
    axes[2].set_ylabel(r"Slow-envelope rate [s$^{-1}$]")
    axes[2].set_xlabel("Airspeed [m/s]")
    axes[2].legend(fontsize=8)

    for ax in axes:
        add_speed_context(ax)
    fig.suptitle("Trim convergence and constraint-moment diagnostics")
    fig.tight_layout()
    savefig(fig, "02_trim_quality_diagnostics.png")


def plot_signal_separation(aux):
    selected = [25.0, 31.0, 38.0, 41.0]
    fig, axes = plt.subplots(len(selected), 2, figsize=(12.0, 9.6), sharex=True)
    for row_idx, speed in enumerate(selected):
        a = aux[speed]
        # Omit the final half-second: zero-phase filtering has no future samples
        # there and its boundary transient is not physical.
        mask = (a["t"] >= 5.0) & (a["t"] <= a["t"][-1] - 0.5)
        for col in range(2):
            axes[row_idx, col].axhline(0, color="0.7", lw=0.7)
            axes[row_idx, col].grid(True, alpha=0.22)
        axes[row_idx, 0].plot(a["t"][mask], a["my"][mask], color="0.55", lw=0.55,
                              label="Raw joint reaction")
        axes[row_idx, 0].plot(a["t"][mask], a["my_slow"][mask], color="#2166ac", lw=1.7,
                              label="Zero-phase low-pass, 2 Hz")
        axes[row_idx, 1].plot(a["t"][mask], a["my_slow"][mask], color="#2166ac", lw=1.4)
        axes[row_idx, 0].set_ylabel(f"{speed:.0f} m/s\n$M_y$ [N m]")
        axes[row_idx, 1].set_ylabel("$M_y$ [N m]")
    axes[0, 0].set_title("Raw signal and slow component")
    axes[0, 1].set_title("Slow component only (independent scale)")
    axes[0, 0].legend(fontsize=8, loc="upper right")
    axes[-1, 0].set_xlabel("Time [s]")
    axes[-1, 1].set_xlabel("Time [s]")
    fig.suptitle("Constraint-moment separation: numerical ripple versus trim transient")
    fig.tight_layout()
    savefig(fig, "03_constraint_moment_signal_separation.png")


def plot_spectrum(aux, rows):
    selected = [25.0, 31.0, 38.0, 41.0]
    colors = plt.cm.viridis(np.linspace(0.1, 0.9, len(selected)))
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.6))
    for speed, color in zip(selected, colors):
        freq = aux[speed]["freq"]
        psd = aux[speed]["psd"]
        area = np.trapz(psd[freq > 0], freq[freq > 0])
        axes[0].semilogy(freq, psd / area, color=color, lw=1.35, label=f"{speed:.0f} m/s")
    axes[0].axvline(0.5, color="0.35", ls=":", label="PID measurement filter: 0.5 Hz")
    axes[0].axvspan(22, 25, color="#762a83", alpha=0.09, label="Near-Nyquist band")
    axes[0].set_xlim(0, 25)
    axes[0].set_ylim(bottom=1e-8)
    axes[0].set_xlabel("Frequency [Hz]")
    axes[0].set_ylabel("Normalised PSD [1/Hz]")
    axes[0].grid(True, which="both", alpha=0.22)
    axes[0].legend(fontsize=8)

    plotted_rows = [r for r in rows if r["speed_mps"] <= 42.0]
    v = np.array([r["speed_mps"] for r in plotted_rows])
    hf = 100 * np.array([r["power_fraction_above_10Hz"] for r in plotted_rows])
    dom = np.array([r["dominant_raw_frequency_Hz"] for r in plotted_rows])
    axes[1].plot(v, hf, "o-", color="#762a83", ms=4, label="Power above 10 Hz")
    axes[1].set_ylabel("High-frequency power [%]", color="#762a83")
    axes[1].tick_params(axis="y", labelcolor="#762a83")
    axes[1].set_xlabel("Airspeed [m/s]")
    axes[1].grid(True, alpha=0.22)
    second = axes[1].twinx()
    second.plot(v, dom, "s", color="#d95f02", ms=3.5, label="Dominant frequency")
    second.axhline(25, color="0.4", ls=":", lw=1)
    second.set_ylabel("Dominant raw frequency [Hz]", color="#d95f02")
    second.tick_params(axis="y", labelcolor="#d95f02")
    axes[1].set_xlim(19.5, 42.2)
    fig.suptitle("Raw constraint-moment spectrum (sampling frequency 50 Hz)")
    fig.tight_layout()
    savefig(fig, "04_constraint_moment_spectrum.png")


def plot_envelopes(aux, rows):
    selected = [30.0, 35.0, 38.0, 40.0, 42.0]
    row_map = {r["speed_mps"]: r for r in rows}
    colors = plt.cm.plasma(np.linspace(0.05, 0.9, len(selected)))
    fig, ax = plt.subplots(figsize=(9.5, 5.8))
    for speed, color in zip(selected, colors):
        a = aux[speed]
        t = a["t"]
        env = a["envelope"]
        reference = np.median(env[(t >= 5.0) & (t <= 6.0)])
        normalised = env / max(reference, 1e-12)
        shown = (t >= 5.0) & (t <= 13.0)
        rate = row_map[speed]["slow_envelope_rate_1_s"]
        r2 = row_map[speed]["slow_envelope_fit_R2"]
        ax.semilogy(t[shown], normalised[shown], color=color, lw=1.6,
                    label=f"{speed:.0f} m/s: rate {rate:+.3f} 1/s, $R^2$={r2:.2f}")
        fit = np.exp(rate * (t[shown] - 5.5))
        level = np.median(normalised[(t >= 5.0) & (t <= 6.0)])
        ax.semilogy(t[shown], level * fit, color=color, ls="--", lw=0.9, alpha=0.8)
    ax.axhline(1, color="0.5", lw=0.8)
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Normalised RMS envelope")
    ax.set_title("Low-frequency constraint-moment transient\n"
                 "solid: measured envelope; dashed: exponential trend")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    savefig(fig, "05_low_frequency_oscillation_envelope.png")


def plot_v30_complete(cases):
    speed = 30.0
    data = cases[speed]
    t, reaction = load_six_component_reaction(
        RESULTS / "V_030.0_mps", data["time_s"]
    )
    fs = 1.0 / float(np.median(np.diff(t)))
    slow = np.column_stack([lowpass(reaction[:, idx], 2.0, fs) for idx in range(6)])
    interior = t <= t[-1] - 0.5
    component_names = ["Fx_N", "Fy_N", "Fz_N", "Mx_Nm", "My_Nm", "Mz_Nm"]

    with (OUT / "v30_reaction_history_slow.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["time_s", *component_names])
        writer.writerows(zip(t[interior], *[slow[interior, idx] for idx in range(6)]))

    final_window = (t >= 13.0) & (t <= 14.5)
    with (OUT / "v30_reaction_summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "component",
                "initial_slow",
                "peak_abs_slow",
                "final_mean_slow_13_to_14p5s",
                "final_std_slow_13_to_14p5s",
                "final_std_raw_13_to_14p5s",
            ]
        )
        for idx, name in enumerate(component_names):
            writer.writerow(
                [
                    name,
                    slow[0, idx],
                    np.max(np.abs(slow[interior, idx])),
                    np.mean(slow[final_window, idx]),
                    np.std(slow[final_window, idx]),
                    np.std(reaction[final_window, idx]),
                ]
            )

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.3), sharex=True)
    axes[0].plot(t[interior], slow[interior, 2], color="#1b7837", lw=1.8)
    axes[1].plot(t[interior], slow[interior, 4], color="#762a83", lw=1.8)
    for ax in axes:
        ax.axhline(0, color="0.35", lw=0.8)
        ax.grid(True, alpha=0.25)
    axes[0].set_ylabel(r"$F_z$ [N]")
    axes[1].set_ylabel(r"$M_y$ [N m]")
    axes[1].set_xlabel("Time [s]")
    fig.suptitle("30 m/s - reaction forces")
    fig.tight_layout()
    savefig(fig, "06_v30_constraint_reactions.png")

    history_t = data["time_s"]
    history_interior = history_t <= history_t[-1] - 0.5
    fz_slow = lowpass(data["Fz_N"], 2.0, fs)
    my_slow = lowpass(data["My_Nm"], 2.0, fs)

    fig, axes = plt.subplots(2, 1, figsize=(10.0, 7.3), sharex=True)
    axes[0].plot(history_t, data["theta_deg"], color="#2166ac", lw=1.8)
    axes[0].set_ylabel(r"$\theta$ [deg]")
    axes[1].plot(history_t, data["delta_elevator_deg"], color="#1b7837", lw=1.8)
    axes[1].set_ylabel(r"$\delta$ control surfaces [deg]")
    axes[1].set_xlabel("Time [s]")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    fig.suptitle("30 m/s — trim variables")
    fig.tight_layout()
    savefig(fig, "07_v30_trim_time_history.png")

    baseline = lowpass(my_slow, 0.15, fs)
    envelope = rms_envelope(my_slow - baseline, fs)
    shown = (history_t >= 2.0) & (history_t <= 13.0)
    fit = (history_t >= 5.0) & (history_t <= 13.0)
    regression = linregress(history_t[fit], np.log(np.maximum(envelope[fit], 1e-9)))
    trend = np.exp(regression.intercept + regression.slope * history_t[shown])

    fig, axes = plt.subplots(2, 1, figsize=(9.5, 7.0), sharex=True)
    axes[0].plot(history_t[history_interior], my_slow[history_interior],
                 color="#762a83", lw=1.8)
    axes[0].axhline(0, color="0.25", lw=0.9)
    axes[0].set_ylabel(r"Slow $M_y$ [N m]")
    axes[0].set_title("Filtered reaction moment: convergence towards zero")
    axes[1].semilogy(history_t[shown], envelope[shown], color="#2166ac", lw=1.8,
                    label="Measured RMS envelope")
    axes[1].semilogy(history_t[shown], trend, color="#b2182b", ls="--", lw=1.4,
                    label=f"Exponential decay: {regression.slope:+.3f} 1/s")
    axes[1].set_ylabel("Oscillation envelope [N m]")
    axes[1].set_xlabel("Time [s]")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.suptitle("30 m/s — decay of the slow pitch-moment transient")
    fig.tight_layout()
    savefig(fig, "08_v30_my_decay.png")


def write_summary(rows):
    fields = list(rows[0])
    with (OUT / "trim_analysis_summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_report(rows):
    by_speed = {r["speed_mps"]: r for r in rows}
    usable = [r["speed_mps"] for r in rows if r["status"] == "usable_trim"]
    near = [r["speed_mps"] for r in rows if r["status"] == "near_BFF_nonstationary"]
    report = f"""# Analisi dello sweep di trim MBDyn

Questa analisi usa esclusivamente `TRIM_PID_VELOCITY/trim_sweep_results`.
La velocità BFF di **{BFF_SPEED:.1f} m/s** è trattata come dato esterno fornito
dall'utente, non come risultato ricavato da questi run.

## Descrizione della simulazione

Lo sweep è costituito da simulazioni aeroelastiche MBDyn di trim, con velocità
imposta e costante per ciascun caso. Ogni simulazione dura **15 s**, usa un
passo temporale di **0.02 s** (frequenza di campionamento 50 Hz) e comprende
gravità, aerodinamica, deformabilità modale e superfici mobili.

Il nodo di riferimento modale è collegato a terra mediante il `total pin joint`
23. Le tre traslazioni, rollio e imbardata sono vincolati; il pitch è imposto
dal comando prodotto dal PID di theta. Di conseguenza Fz e My sono reazioni
del vincolo necessarie a mantenere la configurazione prescritta, non forze
misurate durante un moto rigido completamente libero.

Il controllo viene abilitato a **t = 1 s** con una rampa regolare che termina a
**t = 2 s**. Prima dell'attivazione Fz sostiene essenzialmente il peso del
modello; successivamente i due controllori modificano theta e le superfici
mobili fino a ridurre le reazioni di trim.

### Funzione dei PID

Entrambi i controllori sono puramente integrali: `Kp = 0` e `Kd = 0`.

- Il PID `PID_THETA` usa la reazione verticale **Fz** come errore, con
  `Ki = 2e-5`. La sua uscita impone l'angolo di pitch theta, limitato
  nell'intervallo ±20 deg. Il suo obiettivo è portare Fz a zero, cioè fare in
  modo che la portanza equilibri il peso senza una reazione verticale residua.
- Il PID `PID_ELEVATOR` usa il momento vincolare **My**, con segno di retroazione
  `MY_SIGN = -1` e `Ki = 5e-5`. La sua uscita è il comando collettivo
  `ELEVATOR_DRIVE`, anch'esso limitato a ±20 deg, e mira a portare My a zero.

`ELEVATOR_DRIVE` è applicato simmetricamente a dieci giunti: BFL, BFR,
WF1L-WF4L e WF1R-WF4R. La variabile indicata nei grafici come
`delta control surfaces` rappresenta quindi il comando collettivo di tutte
queste superfici, non un singolo elevatore.

Poiché si tratta di controllori integrali con una misura filtrata, il percorso
verso il trim può presentare piccoli superamenti e inversioni. Tale transitorio
descrive l'algoritmo numerico di ricerca del trim e non deve essere interpretato
come una manovra reale dell'aeromobile.

### Filtri

Nel modello MBDyn, prima di entrare nei PID, sia Fz sia My attraversano un
**Butterworth passa-basso del secondo ordine**, con:

- frequenza di taglio: **0.5 Hz**;
- passo di campionamento: **0.02 s**;
- coefficienti ricorsivi:
  `a1 = 1.911197067426073`, `a2 = -0.914975834801434`;
- coefficienti del numeratore:
  `b0 = 0`, `b1 = 0.000944691843840151`,
  `b2 = 0.001889383687680302`, `b3 = 0.000944691843840151`.

Questo è un filtro causale appartenente alla simulazione: attenua il ripple
della reazione prima che raggiunga i controllori, introducendo anche ritardo di
fase vicino alla frequenza di taglio.

Nei soli grafici delle reazioni è stato inoltre applicato in post-processing un
**Butterworth passa-basso del quarto ordine a 2 Hz, zero-phase**. Questo secondo
filtro non modifica la simulazione e non è usato dai PID: serve esclusivamente
a separare il transitorio lento dal contenuto numerico a 22-25 Hz presente nel
momento vincolare grezzo. Per evitare gli artefatti di bordo del filtraggio
zero-phase, nelle figure viene omesso l'ultimo mezzo secondo. Per il grafico
dell'inviluppo viene anche rimossa la deriva sotto 0.15 Hz e calcolato un
inviluppo RMS su una finestra di 1 s.

## Risultati principali

- Sono stati ricostruiti **{len(rows)} casi**, da {rows[0]['speed_mps']:.0f} a
  {rows[-1]['speed_mps']:.0f} m/s, direttamente dagli `history.csv`.
- Con il criterio |media Fz| < 10 N e lontano dalla zona BFF, i trim più
  utilizzabili sono {min(usable):.0f}–{max(usable):.0f} m/s.
- Fra {min(near):.0f} e {max(near):.0f} m/s il transitorio lento diventa
  non stazionario; i punti non vanno interpretati come equilibri convergenti.
- I casi oltre 42 m/s non sono mostrati nei grafici; rimangono intatti nel CSV
  e nei risultati originali.
- A 31 m/s la deviazione standard del momento grezzo nell'ultimo secondo è
  {by_speed[31.0]['My_raw_std_Nm']:.2f} N m, ma quella della componente sotto
  2 Hz è {by_speed[31.0]['My_slow_std_Nm']:.2f} N m. La differenza è ripple
  ad alta frequenza, non errore quasi-statico di trim.
- Nei casi a bassa e media velocità fino al
  {100*by_speed[31.0]['power_fraction_above_10Hz']:.1f}% della potenza del
  momento grezzo può trovarsi sopra 10 Hz, con picchi tipici fra 22 e 25 Hz.

## Perché il momento vincolare oscilla

`My_Nm` è la reazione istantanea del vincolo che impone il pitch, quindi è un
moltiplicatore di Lagrange sensibile alle accelerazioni modali e al passo di
integrazione. Il modello campiona a 50 Hz (Nyquist 25 Hz), mentre il contenuto
dominante del ripple è spesso 22–25 Hz: soltanto circa due campioni per ciclo.
Questo comportamento è coerente con un modo computazionale quasi alternato
del vincolo, non con una reale oscillazione aerodinamica a quella frequenza.

Inoltre il PID non usa direttamente questo segnale grezzo: nel modello Fz e My
sono filtrati a 0.5 Hz prima del controllo. Per questo il grafico
`03_constraint_moment_signal_separation.png` mostra sia la reazione grezza sia
una stima zero-phase sotto 2 Hz. Quest'ultima è adatta a leggere il transitorio
meccanico lento, ma non viene presentata come una nuova simulazione o come
"dato reale" inventato.

## Superfici mobili effettivamente usate

Il trim non comanda una sola superficie di coda. Tutti i dieci giunti indicati
in `include/control_surfaces.mbd` ricevono lo stesso `ELEVATOR_DRIVE`:
**BFL, BFR, WF1L–WF4L e WF1R–WF4R**. La variabile `delta_elevator_deg` va quindi
letta come comando simmetrico collettivo delle superfici mobili sinistra/destra,
non come deflessione di un elevatore isolato. Non risultano comandi
differenziali di alettoni o un comando separato del timone in questi run.

## Interpretazione del decadimento

`05_low_frequency_oscillation_envelope.png` riporta l'inviluppo RMS della
componente lenta, dopo aver rimosso la deriva quasi-statica sotto 0.15 Hz.
Un tasso negativo indica decadimento del transitorio di trim; un tasso
positivo indica crescita/non convergenza. Le regressioni con basso R² sono
esplicitamente marcate come poco affidabili nel CSV e nel grafico diagnostico.

Questi run bloccano traslazioni e rotazioni e impongono il pitch tramite PID:
il decadimento osservato è quindi quello del **trim vincolato**, non lo
smorzamento di un modo body-freedom libero. Non è corretto stimare da questa
figura una nuova velocità di flutter o pretendere che il cambio di segno
coincida esattamente con 42.3 m/s.

## Verifica numerica consigliata

Per confermare definitivamente l'origine del ripple conviene ripetere un solo
caso rappresentativo (31 m/s) con `dt=0.01 s` e poi `dt=0.005 s`. Se il picco
22–25 Hz si sposta insieme alla frequenza di Nyquist o cala fortemente, la sua
origine numerica è confermata. La soluzione di trim va confrontata usando le
medie lente di Fz/My e non la deviazione standard della reazione grezza.

## File prodotti

1. `01_trim_solution_interpolated.png`: pitch ed elevatore, PCHIP senza
   estrapolazione verso/oltre la BFF.
2. `02_trim_quality_diagnostics.png`: residui, ripple e tasso dell'inviluppo.
3. `03_constraint_moment_signal_separation.png`: segnale grezzo e lento.
4. `04_constraint_moment_spectrum.png`: evidenza spettrale del ripple.
5. `05_low_frequency_oscillation_envelope.png`: decadimento/crescita misurati.
6. `trim_analysis_summary.csv`: metriche numeriche per tutti i 26 casi.
7. `06_v30_constraint_reactions.png`: componenti lente Fz e My del vincolo
   a 30 m/s.
8. `07_v30_trim_time_history.png`: theta e comando collettivo delle superfici
   mobili nel tempo.
9. `08_v30_my_decay.png`: convergenza filtrata di My e inviluppo esponenziale.
10. `v30_reaction_summary.csv` e `v30_reaction_history_slow.csv`: valori delle
    sei componenti e relative statistiche.
"""
    (OUT / "README.md").write_text(report)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.frameon": True,
            "figure.dpi": 120,
        }
    )
    cases = load_cases()
    rows = []
    auxiliary = {}
    for speed, data in cases.items():
        row, aux = analyse_case(speed, data)
        rows.append(row)
        auxiliary[speed] = aux
    rows.sort(key=lambda r: r["speed_mps"])
    write_summary(rows)
    plot_trim_solution(rows)
    plot_quality(rows)
    plot_signal_separation(auxiliary)
    plot_spectrum(auxiliary, rows)
    plot_envelopes(auxiliary, rows)
    plot_v30_complete(cases)
    write_report(rows)
    print(f"Analysed {len(rows)} cases. Results written to {OUT}")


if __name__ == "__main__":
    main()
