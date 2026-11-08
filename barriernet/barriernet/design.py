"""Planning a randomized comparison of two correlated scores (the manuscript's Fig. 2i), by normal theory."""
from __future__ import annotations
from typing import Sequence
import numpy as np
import pandas as pd
from scipy.stats import norm
from .calibrate import events_needed


def detectable_difference(events: float, rho: float, allocation: float = 0.5, power: float = 0.8, alpha: float = 0.05,
                          inflation: float = 1.0) -> float:
    """Smallest difference per SD between two scores' interactions detectable with ``events``."""
    z = norm.ppf(1 - alpha / 2) + norm.ppf(power)
    return z * np.sqrt(inflation * 2 * (1 - rho) / (events * allocation * (1 - allocation)))


def power(events: float, delta: float, rho: float, allocation: float = 0.5, alpha: float = 0.05, inflation: float = 1.0) -> float:
    se = np.sqrt(inflation * 2 * (1 - rho) / (events * allocation * (1 - allocation)))
    return float(norm.sf(norm.ppf(1 - alpha / 2) - delta / se) + norm.cdf(-norm.ppf(1 - alpha / 2) - delta / se))


def curve(rhos: Sequence[float] = tuple(np.round(np.arange(0.30, 0.96, 0.01), 2)), deltas: Sequence[float] = (0.10, 0.15, 0.20),
          allocation: float = 0.5, power_: float = 0.8, alpha: float = 0.05, inflation: float = 1.0) -> pd.DataFrame:
    return pd.DataFrame([dict(rho=r, delta=d, events=events_needed(d, r, allocation, power_, alpha, inflation)) for d in deltas for r in rhos])


__all__ = ["events_needed", "detectable_difference", "power", "curve"]
