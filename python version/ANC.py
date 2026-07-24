"""Single-channel FxLMS demonstration translated from ``ANC.m``."""

from __future__ import annotations

import argparse
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

from ANC_algorithm import ANC_algorithm
from demo_utils import matlab_fir1_bandpass, time_vector


def run_demo(
    *,
    Fs: int = 16_000,
    T: float = 30.0,
    wlen: int = 512,
    slen: int = 256,
    muw: float = 0.001,
    seed: int | None = None,
    plot: bool = True,
) -> dict[str, Any]:
    """Run the same synthetic single-channel setup as the MATLAB script."""

    t = time_vector(Fs, T)
    N = t.size
    rng = np.random.default_rng(seed)

    # Primary and secondary paths: MATLAB fir1(511, ...) and fir1(255, ...).
    Pri_path = matlab_fir1_bandpass(511, 20.0, 6000.0, Fs)
    Sec_path = matlab_fir1_bandpass(255, 20.0, 6000.0, Fs)
    if slen != Sec_path.size:
        raise ValueError(
            f"This demo creates a {Sec_path.size}-tap secondary path; "
            f"received slen={slen}."
        )

    # Broadband reference noise from 100 to 1000 Hz.
    noise = rng.standard_normal(N)
    fil = matlab_fir1_bandpass(127, 100.0, 1000.0, Fs)
    Ref = lfilter(fil, [1.0], noise)
    Dis = lfilter(Pri_path, [1.0], Ref)

    # SciPy has no direct equivalent of MATLAB dsp.FilteredXLMSFilter.  The
    # repository's own FxLMS class is therefore the canonical Python result.
    single_ANC = ANC_algorithm(wlen, slen, Sec_path, Dis, Ref)
    err, single_ANC = single_ANC.ANC_FxLMS(muw)

    if plot:
        sample_count = min(int(round(T * Fs)), N)
        plt.figure()
        plt.plot(t[:sample_count], Dis[:sample_count], label="Disturbance")
        plt.plot(t[:sample_count], err[:sample_count], label="Residual error signal")
        plt.title("Active Noise Control")
        plt.xlabel("Time (seconds)")
        plt.ylabel("Signal value")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    return {
        "t": t,
        "Ref": Ref,
        "Dis": Dis,
        "err": err,
        "controller": single_ANC,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration", type=float, default=30.0, help="Duration in s.")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_demo(T=args.duration, seed=args.seed, plot=not args.no_plot)
