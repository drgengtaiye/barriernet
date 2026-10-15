from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
from typing import Optional, List, Tuple
from collections import deque
from dataclasses import dataclass

# Visium array coordinates (obs['array_row'], obs['array_col']) are "doubled column":
# only positions with row + col even hold a spot, and a spot's six neighbours sit at
# (0, +-2) and (+-1, +-1). Versions <= 0.2.0 used axial offsets here; on Visium data only
# the (1,-1)/(-1,1) steps ever matched, so each spot had at most two neighbours and
# rings grew along a single diagonal (see CHANGELOG.md, 0.2.1).
VISIUM_STEPS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
# True axial hexagonal coordinates (not what Visium exports), kept for other platforms.
AXIAL_STEPS = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
HEX_STEPS = VISIUM_STEPS  # backward-compatible name; now correct for Visium
_STEPS = {"visium": VISIUM_STEPS, "axial": AXIAL_STEPS}


def pick_first_existing(obs: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    for c in candidates:
        if c in obs.columns:
            return c
    return None


def ensure_tls_score(obs: pd.DataFrame) -> str:
    for c in ["TLS_niche_score", "TLS_core_score", "TLS_score_score", "TLS_score"]:
        if c in obs.columns:
            return c
    if all([(f"tg_{x}" in obs.columns) for x in ["Bcell", "Plasma", "DC"]]):
        obs["TLS_core_score"] = (pd.to_numeric(obs["tg_Bcell"], errors="coerce") +
                                 pd.to_numeric(obs["tg_Plasma"], errors="coerce") +
                                 pd.to_numeric(obs["tg_DC"], errors="coerce"))
        return "TLS_core_score"
    raise KeyError("Cannot find TLS score in obs.")


def build_hex_graph(obs: pd.DataFrame, coords: str = "visium") -> Tuple[list[list[int]], np.ndarray]:
    """Neighbour lists on the spot lattice, restricted to in-tissue spots.

    coords: "visium" (default; array_row/array_col as exported by Space Ranger) or
    "axial" (true axial hexagonal coordinates).
    """
    if coords not in _STEPS:
        raise ValueError("coords must be 'visium' or 'axial'")
    if ("array_row" not in obs.columns) or ("array_col" not in obs.columns):
        raise KeyError("Need obs['array_row'] and obs['array_col'] for hex-grid adjacency.")
    ar = pd.to_numeric(obs["array_row"], errors="coerce").fillna(-1).astype(int).to_numpy()
    ac = pd.to_numeric(obs["array_col"], errors="coerce").fillna(-1).astype(int).to_numpy()

    keep = np.ones(len(obs), dtype=bool)
    if "in_tissue" in obs.columns:
        keep = (pd.to_numeric(obs["in_tissue"], errors="coerce").fillna(0).astype(int).to_numpy() == 1)
    valid = keep & (ar >= 0) & (ac >= 0)

    if coords == "visium" and valid.any() and len(np.unique((ar[valid] + ac[valid]) % 2)) > 1:
        warnings.warn("array_row + array_col has mixed parity; these do not look like Visium "
                      "array coordinates. Check `coords`.")

    steps = _STEPS[coords]
    pos2i = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(ar, ac)) if valid[i]}
    nbrs = [[] for _ in range(len(obs))]
    for i, (r, c) in enumerate(zip(ar, ac)):
        if not valid[i]:
            continue
        for dr, dc in steps:
            j = pos2i.get((r + dr, c + dc))
            if j is not None:
                nbrs[i].append(j)

    deg = np.array([len(n) for n in nbrs])[valid]
    if deg.size and deg.mean() < 3:
        warnings.warn(f"Mean lattice degree is {deg.mean():.2f} (expected about 6 for a Visium "
                      f"section); the coordinate system is probably wrong (coords='{coords}').")
    return nbrs, keep


def lattice_degree(nbrs: list[list[int]], keep: np.ndarray) -> np.ndarray:
    return np.array([len(nbrs[i]) if keep[i] else 0 for i in range(len(nbrs))])


def bfs_dist(nbrs: list[list[int]], keep: np.ndarray, sources_mask: np.ndarray) -> np.ndarray:
    dist = np.full(len(nbrs), np.inf)
    q = deque()
    src = np.where(keep & sources_mask)[0]
    for s in src:
        dist[s] = 0
        q.append(s)
    while q:
        u = q.popleft()
        for v in nbrs[u]:
            if dist[v] > dist[u] + 1:
                dist[v] = dist[u] + 1
                q.append(v)
    return dist


def edge_dist(nbrs: list[list[int]], keep: np.ndarray) -> np.ndarray:
    """Hops to the nearest tissue-edge spot (a spot with fewer than six in-tissue neighbours)."""
    edge = keep & (lattice_degree(nbrs, keep) < 6)
    return bfs_dist(nbrs, keep, edge)


def local_maxima_peaks(xy: np.ndarray, score: np.ndarray, top_q: float = 0.95, k: int = 9) -> Tuple[np.ndarray, float, int]:
    """Spots at or above the `top_q` quantile that are not smaller than any of their k nearest
    spots in `xy` (the k include the spot itself, i.e. k - 1 neighbours)."""
    score = np.asarray(score, dtype=float)
    finite = score[np.isfinite(score)]
    thr = float(np.nanquantile(finite, top_q)) if finite.size else np.nan
    cand = np.isfinite(score) & (score >= thr)

    n = xy.shape[0]
    pk = np.zeros(n, dtype=bool)
    if n == 0 or not cand.any():
        return pk, thr, 0

    for i in np.where(cand)[0]:
        d = np.sum((xy - xy[i]) ** 2, axis=1)
        nn = np.argpartition(d, kth=min(k, n - 1))[:min(k, n)]
        nn = nn[cand[nn]]
        if nn.size and score[i] >= np.nanmax(score[nn]):
            pk[i] = True
    return pk, thr, int(pk.sum())


def lattice_maxima_peaks(nbrs: list[list[int]], keep: np.ndarray, score: np.ndarray,
                         top_q: float = 0.95) -> Tuple[np.ndarray, float, int]:
    """Spots at or above the `top_q` quantile that are not smaller than any lattice neighbour."""
    score = np.asarray(score, dtype=float)
    finite = score[np.isfinite(score) & keep]
    thr = float(np.nanquantile(finite, top_q)) if finite.size else np.nan
    cand = keep & np.isfinite(score) & (score >= thr)
    pk = np.zeros(len(score), dtype=bool)
    for i in np.where(cand)[0]:
        nb = [j for j in nbrs[i] if np.isfinite(score[j])]
        pk[i] = all(score[i] >= score[j] for j in nb)
    return pk, thr, int(pk.sum())


@dataclass
class SpatialKeys:
    tls: str
    ecm: Optional[str]
    tam: Optional[str]
    stf3: Optional[str]


def infer_spatial_keys(obs: pd.DataFrame) -> SpatialKeys:
    tls = ensure_tls_score(obs)
    ecm = pick_first_existing(obs, ["state_CAF_ECM", "prog_CAF_ECM", "score_CAF_ECM"])
    tam = pick_first_existing(obs, ["state_SPP1_TAM", "state_TAM_SUP", "prog_TAM_SUP", "score_TAM_SUP"])
    stf3 = pick_first_existing(obs, ["state_STF3_FAP_INHBA", "prog_STF3_FAP_INHBA", "prog_CAF_STF3_FAP_INHBA"])
    return SpatialKeys(tls=tls, ecm=ecm, tam=tam, stf3=stf3)


def assign_rings(obs: pd.DataFrame,
                 spatial_xy: np.ndarray,
                 tls_key: str,
                 peak_q: float = 0.95,
                 peak_k: int = 9,
                 far_min_dist: int = 3,
                 far_exclude_edge: bool = False,
                 edge_buffer: int = 0,
                 far_include_unreachable: bool = True,
                 peaks: Optional[np.ndarray] = None,
                 peak_method: str = "knn",
                 coords: str = "visium") -> pd.DataFrame:
    """TLS-centred zones: Hotspot (TLS peak), Ring1 (1 hop), Ring2 (2 hops), Far.

    far_min_dist            Far = at least this many hops from a TLS peak (default 3).
    far_exclude_edge        drop spots within `edge_buffer` hops of the tissue edge from Far
                            (edge = fewer than six in-tissue neighbours; buffer 0 = edge spots only).
    far_include_unreachable spots in tissue islands with no TLS peak (infinite distance) count
                            as Far (default True, as in the manuscript analyses).
    peaks                   optional precomputed boolean peak mask; otherwise peaks are called
                            on `tls_key` at the `peak_q` quantile, as local maxima among the
                            `peak_k` nearest spots in `spatial_xy` (peak_method="knn", default) or
                            among the six lattice neighbours (peak_method="lattice").
    Note: the TLS peaks stored with the v13 analyses were lattice maxima on the 0.2.0
    (broken) graph, i.e. against at most two diagonal neighbours; see CHANGELOG.md.
    """
    tls = pd.to_numeric(obs[tls_key], errors="coerce").to_numpy(dtype=float)
    nbrs, keep = build_hex_graph(obs, coords=coords)
    if peaks is not None:
        tls_pk = np.asarray(peaks).astype(bool)
        tls_thr, tls_n = np.nan, int(tls_pk.sum())
    elif peak_method == "knn":
        tls_pk, tls_thr, tls_n = local_maxima_peaks(spatial_xy, tls, top_q=peak_q, k=peak_k)
    elif peak_method == "lattice":
        tls_pk, tls_thr, tls_n = lattice_maxima_peaks(nbrs, keep, tls, top_q=peak_q)
    else:
        raise ValueError("peak_method must be 'knn' or 'lattice'")

    dist_edge = edge_dist(nbrs, keep)
    dist_tls = bfs_dist(nbrs, keep, tls_pk & keep)

    ring = np.full(len(obs), "NA", dtype=object)
    ring[keep] = "Other"
    ring[tls_pk & keep] = "Hotspot"
    ring[keep & (dist_tls == 1)] = "Ring1"
    ring[keep & (dist_tls == 2)] = "Ring2"

    far = keep & (dist_tls >= far_min_dist)
    if not far_include_unreachable:
        far &= np.isfinite(dist_tls)
    if far_exclude_edge:
        far &= dist_edge > edge_buffer
    if far.sum() < 30:
        warnings.warn(f"Only {int(far.sum())} Far spots; Ring - Far contrasts will be noisy.")
    ring[far] = "Far"

    return pd.DataFrame({
        "tls_peak": tls_pk.astype(int),
        "tls_peak_thr": tls_thr,
        "tls_peak_n": tls_n,
        "dist_tls_hops": dist_tls,
        "dist_edge_hops": dist_edge,
        "lattice_degree": lattice_degree(nbrs, keep),
        "ring": ring
    }, index=obs.index)


def ring_summary(obs: pd.DataFrame, ring_df: pd.DataFrame, keys: SpatialKeys) -> pd.DataFrame:
    df = pd.concat([obs[[keys.tls]].copy(), ring_df], axis=1)

    def add(key: Optional[str], label: str):
        if key is None:
            df[label] = np.nan
        else:
            df[label] = pd.to_numeric(obs[key], errors="coerce")

    add(keys.ecm, "ECM")
    add(keys.tam, "TAM")
    add(keys.stf3, "STF3")

    ring_order = ["Hotspot", "Ring1", "Ring2", "Far"]
    rows = []
    for r in ring_order:
        sub = df[df["ring"] == r]
        if sub.shape[0] == 0:
            continue
        rows.append({
            "ring": r,
            "n": int(sub.shape[0]),
            "TLS_mean": float(np.nanmean(sub[keys.tls].values)),
            "ECM_mean": float(np.nanmean(sub["ECM"].values)),
            "TAM_mean": float(np.nanmean(sub["TAM"].values)),
            "STF3_mean": float(np.nanmean(sub["STF3"].values)),
        })
    out = pd.DataFrame(rows)

    def get(r, col):
        s = out.loc[out["ring"] == r, col]
        return float(s.iloc[0]) if len(s) else np.nan

    for metric, col in [("ECM", "ECM_mean"), ("TAM", "TAM_mean"), ("STF3", "STF3_mean")]:
        out[f"{metric}_Ring1_minus_Far"] = get("Ring1", col) - get("Far", col)
        out[f"{metric}_Ring2_minus_Far"] = get("Ring2", col) - get("Far", col)

    return out


def barrier_index(section_table: pd.DataFrame,
                  ecm_col: str = "ECM_Ring1_minus_Far",
                  tam_col: str = "TAM_Ring1_minus_Far") -> pd.Series:
    """Section-level BarrierIndex: z(ECM Ring1-Far) + z(TAM Ring1-Far), each standardised
    across the sections in `section_table` (population SD, ddof = 0)."""
    z = lambda v: (v - v.mean()) / v.std(ddof=0)
    return z(section_table[ecm_col].astype(float)) + z(section_table[tam_col].astype(float))
