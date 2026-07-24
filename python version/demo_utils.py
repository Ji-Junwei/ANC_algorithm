"""Shared utilities used by the MATLAB-style Python demonstration scripts."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.signal import firwin


def matlab_fir1_bandpass(
    order: int,
    low_hz: float,
    high_hz: float,
    fs: float,
) -> NDArray[np.float64]:
    """Approximate MATLAB ``fir1(order, [low high] / (Fs/2))``.

    MATLAB ``fir1`` and SciPy ``firwin`` both use an order-plus-one FIR length.
    A Hamming window and passband scaling are specified explicitly here.  Very
    small coefficient differences can still occur between library versions.
    """

    if not 0 < low_hz < high_hz < fs / 2:
        raise ValueError("Require 0 < low_hz < high_hz < fs/2.")
    if order < 1:
        raise ValueError("FIR order must be at least 1.")

    return firwin(
        numtaps=order + 1,
        cutoff=[low_hz, high_hz],
        pass_zero=False,
        window="hamming",
        scale=True,
        fs=fs,
    )


def awgn_measured(
    signal: ArrayLike,
    snr_db: float,
    rng: np.random.Generator,
) -> NDArray:
    """Add real AWGN using measured signal power, like MATLAB ``awgn``."""

    x = np.asarray(signal)
    signal_power = float(np.mean(np.abs(x) ** 2))
    if signal_power == 0.0:
        return x.copy()
    noise_power = signal_power / (10.0 ** (float(snr_db) / 10.0))
    noise = rng.standard_normal(x.shape) * np.sqrt(noise_power)
    return x + noise


def moving_average(signal: ArrayLike, window: int) -> NDArray[np.float64]:
    """Centered moving average for the demonstration plots."""

    x = np.asarray(signal, dtype=np.float64).reshape(-1)
    if x.size == 0:
        return x.copy()
    width = max(1, min(int(window), x.size))
    kernel = np.ones(width, dtype=np.float64)
    numerator = np.convolve(x, kernel, mode="same")
    denominator = np.convolve(np.ones_like(x), kernel, mode="same")
    return numerator / denominator


def time_vector(fs: float, duration_s: float) -> NDArray[np.float64]:
    """Match MATLAB ``0:1/Fs:T`` including the final endpoint."""

    sample_count = int(round(float(fs) * float(duration_s))) + 1
    return np.arange(sample_count, dtype=np.float64) / float(fs)
