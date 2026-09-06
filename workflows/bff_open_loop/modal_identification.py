"""Frequency-agnostic multimodal Matrix Pencil identification utilities."""

from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

WINDOWS_S = (0.5, 0.75, 1.0, 1.5, 2.0, 2.5)
DETRENDS = ("raw", "linear", "poly2")
BASE_ORDERS = (4, 6, 8, 10, 12, 14)


def detrend_variant(time: np.ndarray, signal: np.ndarray, kind: str) -> np.ndarray:
    t = np.asarray(time, float) - float(time[0])
    y = np.asarray(signal, float)
    if kind == "raw":
        return y
    degree = 1 if kind == "linear" else 2
    return y - np.polyval(np.polyfit(t, y, degree), t)


def matrix_pencil(signal: np.ndarray, dt: float, order: int) -> tuple[np.ndarray, np.ndarray, float]:
    """Return continuous poles, complex amplitudes and normalized residual."""
    y = np.asarray(signal, complex)
    n = len(y)
    rows = max(order + 2, min(n // 2, 180))
    columns = n - rows + 1
    if n < 20 or order >= min(rows, columns - 1):
        return np.array([], complex), np.array([], complex), math.inf
    hankel = np.lib.stride_tricks.sliding_window_view(y, columns).T
    h0, h1 = hankel[:, :-1], hankel[:, 1:]
    u, singular, vh = np.linalg.svd(h0, full_matrices=False)
    if singular[0] <= 0.0:
        return np.array([], complex), np.array([], complex), math.inf
    numerical_rank = int(np.count_nonzero(singular > singular[0] * 1e-10))
    rank = min(order, numerical_rank)
    if rank < 2:
        return np.array([], complex), np.array([], complex), math.inf
    v = vh[:rank].conj().T
    operator = u[:, :rank].conj().T @ h1 @ v @ np.diag(1.0 / singular[:rank])
    discrete = np.linalg.eigvals(operator)
    keep = np.isfinite(discrete) & (np.abs(discrete) > 1e-12)
    discrete = discrete[keep]
    poles = np.log(discrete) / dt
    vandermonde = discrete[None, :] ** np.arange(n)[:, None]
    amplitudes = np.linalg.lstsq(vandermonde, y, rcond=None)[0]
    reconstructed = vandermonde @ amplitudes
    residual = float(np.linalg.norm(y - reconstructed) / max(np.linalg.norm(y), 1e-30))
    return poles, amplitudes, residual


def extract_pole_rows(
    signal_name: str,
    time: np.ndarray,
    signal: np.ndarray,
    window_s: float,
    local_valid: bool,
) -> list[dict]:
    rows: list[dict] = []
    dt = float(np.median(np.diff(time)))
    max_order = min(14, max(4, len(time) // 4))
    orders = [value for value in BASE_ORDERS if value <= max_order]
    for detrending in DETRENDS:
        value = detrend_variant(time, signal, detrending)
        scale = float(np.sqrt(np.mean(np.abs(value) ** 2)))
        if not np.isfinite(scale) or scale < 1e-12:
            continue
        normalized = value / scale
        for order in orders:
            poles, amplitudes, residual = matrix_pencil(normalized, dt, order)
            for pole, amplitude in zip(poles, amplitudes):
                frequency = float(abs(pole.imag) / (2.0 * math.pi))
                sigma = float(pole.real)
                # Keep one member of each complex pair and all real poles.
                if pole.imag < -1e-7 or frequency > 0.45 / dt or abs(sigma) > 30.0:
                    continue
                omega_n = math.hypot(sigma, 2.0 * math.pi * frequency)
                rows.append({
                    "signal": signal_name, "window_s": window_s,
                    "local_window_valid": local_valid, "detrend": detrending,
                    "order": order, "sigma_per_s": sigma, "frequency_hz": frequency,
                    "damping_ratio": None if omega_n == 0.0 else -sigma / omega_n,
                    "amplitude_real": float(amplitude.real * scale),
                    "amplitude_imag": float(amplitude.imag * scale),
                    "amplitude_abs": float(abs(amplitude) * scale),
                    # Pair amplitude relative to RMS-normalized input. This is
                    # robust to high-amplitude cancelling nuisance poles.
                    "participation": float(min(1.0, 2.0 * abs(amplitude))),
                    "fit_residual": residual,
                    "is_complex": bool(frequency >= 0.15),
                })
    return rows


def _near(a: dict, b: dict) -> bool:
    fmean = 0.5 * (a["frequency_hz"] + b["frequency_hz"])
    smean = 0.5 * (abs(a["sigma_per_s"]) + abs(b["sigma_per_s"]))
    return (
        abs(a["frequency_hz"] - b["frequency_hz"]) <= max(0.18, 0.07 * fmean)
        and abs(a["sigma_per_s"] - b["sigma_per_s"]) <= max(0.65, 0.20 * smean)
    )


def cluster_poles(rows: list[dict], signal_name: str) -> list[dict]:
    source = [row for row in rows if row["signal"] == signal_name and row["is_complex"] and row["local_window_valid"]]
    if not source:
        return []
    parent = list(range(len(source)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        a, b = find(i), find(j)
        if a != b:
            parent[b] = a

    for i in range(len(source)):
        for j in range(i):
            if _near(source[i], source[j]):
                union(i, j)
    grouped: dict[int, list[dict]] = defaultdict(list)
    for i, row in enumerate(source):
        grouped[find(i)].append(row)

    valid_windows = {row["window_s"] for row in source}
    valid_orders = {row["order"] for row in source}
    clusters: list[dict] = []
    for members in grouped.values():
        if len(members) < 3:
            continue
        frequencies = np.array([row["frequency_hz"] for row in members])
        sigmas = np.array([row["sigma_per_s"] for row in members])
        median_f, median_s = float(np.median(frequencies)), float(np.median(sigmas))
        core = [row for row in members if abs(row["frequency_hz"] - median_f) <= max(0.25, 0.10 * median_f) and abs(row["sigma_per_s"] - median_s) <= max(1.0, 0.35 * abs(median_s))]
        if len(core) < 3:
            continue
        frequencies = np.array([row["frequency_hz"] for row in core])
        sigmas = np.array([row["sigma_per_s"] for row in core])
        windows = sorted({row["window_s"] for row in core})
        orders = sorted({row["order"] for row in core})
        detrends = sorted({row["detrend"] for row in core})
        participation = float(np.median([row["participation"] for row in core]))
        f_std, sigma_std = float(np.std(frequencies)), float(np.std(sigmas))
        window_persistence = len(windows) / max(len(valid_windows), 1)
        order_persistence = len(orders) / max(len(valid_orders), 1)
        detrend_persistence = len(detrends) / len(DETRENDS)
        stability = math.exp(-f_std / max(0.12, 0.08 * float(np.median(frequencies)))) * math.exp(-sigma_std / max(0.7, 0.25 * abs(float(np.median(sigmas)))))
        participation_score = min(1.0, participation / 0.10)
        quality = 0.30 * window_persistence + 0.25 * order_persistence + 0.15 * detrend_persistence + 0.15 * stability + 0.15 * participation_score
        physical = bool(
            len(windows) >= min(2, len(valid_windows)) and len(orders) >= min(3, len(valid_orders))
            and len(detrends) >= 2 and participation >= 0.01
            and f_std <= max(0.22, 0.10 * float(np.median(frequencies)))
            and sigma_std <= max(1.2, 0.40 * abs(float(np.median(sigmas))))
        )
        clusters.append({
            "signal": signal_name, "cluster_id": -1,
            "frequency_hz": float(np.median(frequencies)), "frequency_std_hz": f_std,
            "sigma_per_s": float(np.median(sigmas)), "sigma_std_per_s": sigma_std,
            "damping_ratio": float(np.median([row["damping_ratio"] for row in core if row["damping_ratio"] is not None])),
            "participation_median": participation, "window_persistence": window_persistence,
            "order_persistence": order_persistence, "detrend_persistence": detrend_persistence,
            "windows_s": windows, "orders": orders, "detrends": detrends,
            "observations": len(core), "quality": quality, "physical": physical,
        })
    clusters.sort(key=lambda item: (item["physical"], item["quality"]), reverse=True)
    for index, cluster in enumerate(clusters, 1):
        cluster["cluster_id"] = index
    return clusters


def cluster_distance(a: dict, b: dict) -> tuple[float, float, float]:
    df = abs(a["frequency_hz"] - b["frequency_hz"])
    ds = abs(a["sigma_per_s"] - b["sigma_per_s"])
    cost = df / max(0.20, 0.08 * 0.5 * (a["frequency_hz"] + b["frequency_hz"])) + ds / max(0.9, 0.30 * 0.5 * (abs(a["sigma_per_s"]) + abs(b["sigma_per_s"])))
    return cost, df, ds


def match_bff_candidate(
    clusters: dict[str, list[dict]],
    reference_frequency_hz: float | None = None,
) -> dict:
    # Identify bending from a measured structural quantity, not from rigid-body
    # pitch rate.  The symmetric physical wing-tip displacement is the primary
    # observable; the SWB1 modal coordinate and its velocity must independently
    # contain the same persistent complex pair.
    required = ("symmetric_tip", "swb1", "swb1_velocity")
    if any(not clusters.get(name) for name in required):
        return {"accepted": False, "confidence": "none", "reason": "missing persistent clusters in one or more required signals", "matches": {}}
    hypotheses = []
    for tip_cluster in clusters["symmetric_tip"]:
        matches = {"symmetric_tip": tip_cluster}
        total_cost = 0.0
        compatible = True
        for name in required[1:]:
            options = [(cluster_distance(tip_cluster, candidate), candidate) for candidate in clusters[name]]
            (cost, df, ds), chosen = min(options, key=lambda item: item[0][0])
            if df > max(0.30, 0.12 * tip_cluster["frequency_hz"]) or ds > max(1.8, 0.50 * abs(tip_cluster["sigma_per_s"])):
                compatible = False
            matches[name] = chosen
            total_cost += cost
        qualities = [value["quality"] for value in matches.values()]
        participations = [value["participation_median"] for value in matches.values()]
        physical = all(matches[name]["physical"] for name in required)
        required_participation = min(matches[name]["participation_median"] for name in required)
        # Compare the mean normalized mismatch, not its raw sum: adding an
        # independent corroborating signal must not lower confidence merely
        # because another distance term exists.
        comparison_weight = float(len(required) - 1)
        score = float(
            np.mean(qualities)
            * math.exp(-total_cost / (3.0 * comparison_weight))
            * min(1.0, required_participation / 0.02)
        )
        raw_score = score
        if reference_frequency_hz is not None:
            # Adjacent-speed continuation prevents a higher structural mode
            # from stealing the tracked SWB1 branch when its local fit happens
            # to be cleaner in one case.
            tracking_scale = max(0.45, 0.18 * reference_frequency_hz)
            tracking_difference = abs(tip_cluster["frequency_hz"] - reference_frequency_hz)
            score *= math.exp(-(tracking_difference / tracking_scale) ** 2)
        else:
            tracking_difference = None
        hypotheses.append((score, compatible and physical, matches, total_cost, raw_score, tracking_difference))
    score, compatible_physical, matches, cost, raw_score, tracking_difference = max(hypotheses, key=lambda item: item[0])
    weights = np.array([max(value["quality"], 1e-6) for value in matches.values()])
    frequency = float(np.average([value["frequency_hz"] for value in matches.values()], weights=weights))
    sigma = float(np.average([value["sigma_per_s"] for value in matches.values()], weights=weights))
    omega_n = math.hypot(sigma, 2.0 * math.pi * frequency)
    rigid_matches = {}
    for name in ("q", "alpha", "pitch"):
        if not clusters.get(name):
            continue
        (_, df, ds), rigid = min(
            ((cluster_distance({"frequency_hz": frequency, "sigma_per_s": sigma}, value), value) for value in clusters[name]),
            key=lambda item: item[0][0],
        )
        if (
            rigid["physical"]
            and df <= max(0.35, 0.12 * frequency)
            and ds <= max(1.8, 0.50 * abs(sigma))
        ):
            rigid_matches[name] = rigid
    # Keep rigid participation as a diagnostic of body-freedom coupling, but
    # do not use q/alpha to select or veto the bending branch.  Those channels
    # are reserved for the independent short-period identification below; the
    # BFF observable is the persistent physical tip/SWB1/SWB1dot pole.
    body_freedom_coupling = bool("q" in rigid_matches or "alpha" in rigid_matches)
    accepted = bool(compatible_physical and raw_score >= 0.55)
    confidence = "high" if accepted and score >= 0.75 else "medium" if accepted else "low" if score >= 0.30 else "none"
    reason = (
        "persistent compatible tip/SWB1/SWB1dot pair"
        if accepted else
        "best structural hypothesis rejected: insufficient persistence or compatibility"
    )
    return {
        "accepted": accepted, "confidence": confidence, "quality_score": raw_score,
        "reason": reason, "frequency_hz": frequency, "sigma_per_s": sigma,
        "damping_ratio": -sigma / omega_n if omega_n else None, "matching_cost": cost,
        "tracking_reference_frequency_hz": reference_frequency_hz,
        "tracking_frequency_difference_hz": tracking_difference,
        "body_freedom_coupling_confirmed": body_freedom_coupling,
        "rigid_body_matches": {name: {key: value[key] for key in ("cluster_id", "frequency_hz", "sigma_per_s", "participation_median", "quality", "physical")} for name, value in rigid_matches.items()},
        "matches": {name: {key: value[key] for key in ("cluster_id", "frequency_hz", "sigma_per_s", "participation_median", "quality", "physical")} for name, value in matches.items()},
    }


def match_short_period_candidate(clusters: dict[str, list[dict]], bff_candidate: dict) -> dict:
    """Find the rigid longitudinal pair without using it to select the BFF.

    Pitch rate, effective angle of attack and pitch attitude are deliberately
    kept separate from the structural tip/SWB observables.  Near flutter the
    two hypotheses may converge; ``coupled_to_bff`` records that fact rather
    than forcing one label onto the common pole.
    """
    required = ("q", "alpha", "pitch")
    if any(not clusters.get(name) for name in required):
        return {
            "accepted": False, "confidence": "none",
            "reason": "missing persistent clusters in q, alpha or pitch",
            "matches": {},
        }
    hypotheses = []
    for q_cluster in clusters["q"]:
        matches = {"q": q_cluster}
        total_cost = 0.0
        compatible = True
        for name in required[1:]:
            options = [(cluster_distance(q_cluster, candidate), candidate) for candidate in clusters[name]]
            (cost, df, ds), chosen = min(options, key=lambda item: item[0][0])
            if df > max(0.35, 0.10 * q_cluster["frequency_hz"]) or ds > max(1.8, 0.50 * abs(q_cluster["sigma_per_s"])):
                compatible = False
            matches[name] = chosen
            total_cost += cost
        physical = all(matches[name]["physical"] for name in required)
        participation = min(matches[name]["participation_median"] for name in required)
        score = float(
            np.mean([matches[name]["quality"] for name in required])
            * math.exp(-total_cost / (3.0 * (len(required) - 1)))
            * min(1.0, participation / 0.02)
        )
        hypotheses.append((score, compatible and physical, matches, total_cost))
    score, compatible_physical, matches, cost = max(hypotheses, key=lambda item: item[0])
    weights = np.array([max(value["quality"], 1e-6) for value in matches.values()])
    frequency = float(np.average([value["frequency_hz"] for value in matches.values()], weights=weights))
    sigma = float(np.average([value["sigma_per_s"] for value in matches.values()], weights=weights))
    omega_n = math.hypot(sigma, 2.0 * math.pi * frequency)
    accepted = bool(compatible_physical and score >= 0.55)
    coupled = bool(
        accepted and bff_candidate.get("accepted", False)
        and abs(frequency - bff_candidate["frequency_hz"]) <= max(0.35, 0.12 * bff_candidate["frequency_hz"])
        and abs(sigma - bff_candidate["sigma_per_s"]) <= max(1.8, 0.50 * abs(bff_candidate["sigma_per_s"]))
    )
    structural_matches = {}
    for name in ("symmetric_tip", "swb1", "swb1_velocity"):
        if not clusters.get(name):
            continue
        (_, df, ds), structural = min(
            ((cluster_distance({"frequency_hz": frequency, "sigma_per_s": sigma}, value), value) for value in clusters[name]),
            key=lambda item: item[0][0],
        )
        if structural["physical"] and df <= max(0.35, 0.10 * frequency) and ds <= max(1.8, 0.50 * abs(sigma)):
            structural_matches[name] = structural
    structural_overlap = len(structural_matches) >= 2
    isolated = bool(accepted and not structural_overlap)
    confidence = "high" if accepted and score >= 0.75 else "medium" if accepted else "low" if score >= 0.30 else "none"
    if accepted and coupled:
        reason = "persistent q/alpha/pitch pair coincident with the structural BFF pole"
    elif accepted and structural_overlap:
        reason = "persistent q/alpha/pitch pair also present in structural channels; not an isolated short-period identification"
    elif accepted:
        reason = "persistent rigid-longitudinal pair in q, effective alpha and pitch"
    else:
        reason = "best q/alpha/pitch hypothesis rejected for insufficient persistence or compatibility"
    return {
        "accepted": accepted, "confidence": confidence, "quality_score": score,
        "reason": reason,
        "mode_label": "coupled_bff_short_period" if coupled else "rigid_flex_longitudinal_observable" if structural_overlap else "isolated_short_period_observable",
        "coupled_to_bff": coupled, "structural_overlap": structural_overlap,
        "short_period_isolated": isolated, "frequency_hz": frequency, "sigma_per_s": sigma,
        "damping_ratio": -sigma / omega_n if omega_n else None, "matching_cost": cost,
        "coincident_structural_matches": {name: {key: value[key] for key in ("cluster_id", "frequency_hz", "sigma_per_s", "participation_median", "quality", "physical")} for name, value in structural_matches.items()},
        "matches": {name: {key: value[key] for key in ("cluster_id", "frequency_hz", "sigma_per_s", "participation_median", "quality", "physical")} for name, value in matches.items()},
    }


def identify_multimodal(
    time: np.ndarray,
    signals: dict[str, np.ndarray],
    local_valid_by_window: dict[float, bool],
    tracking_reference_frequency_hz: float | None = None,
) -> tuple[list[dict], dict[str, list[dict]], dict, dict]:
    all_rows: list[dict] = []
    t0 = float(time[0])
    for window in WINDOWS_S:
        use = (time >= t0) & (time <= t0 + window + 1e-9)
        if np.count_nonzero(use) < 20:
            continue
        for name, signal in signals.items():
            all_rows.extend(extract_pole_rows(name, time[use], signal[use], window, local_valid_by_window.get(window, False)))
    clusters = {name: cluster_poles(all_rows, name) for name in signals}
    candidate = match_bff_candidate(clusters, tracking_reference_frequency_hz)
    short_period = match_short_period_candidate(clusters, candidate)
    return all_rows, clusters, candidate, short_period
