"""Fully connected J x K x M FxLMS demo translated from ``MCFxLMS_tst.m``."""

from __future__ import annotations

import argparse
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

from MultiChannelFxLMS import MultiChannelFxLMS
from demo_utils import matlab_fir1_bandpass, time_vector


def run_demo(
    *,
    Fs: int = 16_000,
    T: float = 30.0,
    J: int = 6,
    K: int = 6,
    M: int = 6,
    wLen: int = 512,
    sLen: int = 256,
    muw: float = 1e-6,
    backend: str = "numpy",
    device: str = "auto",
    dtype: str = "float64",
    seed: int | None = None,
    plot: bool = True,
) -> dict[str, Any]:
    """Run the original fully connected synthetic ANC configuration."""

    t = time_vector(Fs, T)
    N = t.size
    rng = np.random.default_rng(seed)

    Pri_path = matlab_fir1_bandpass(511, 80.0, 5000.0, Fs)
    PrimaryPath = np.tile(Pri_path[None, None, :], (M, J, 1))
    Sec_path = matlab_fir1_bandpass(sLen - 1, 80.0, 5000.0, Fs)
    SecondaryPath = np.tile(Sec_path[None, None, :], (M, K, 1))

    noise = rng.standard_normal(N)
    fil = matlab_fir1_bandpass(63, 100.0, 1000.0, Fs)
    X = lfilter(fil, [1.0], noise)
    Ref = np.tile(X[None, :], (J, 1))

    Dis = np.zeros((M, N), dtype=np.float64)
    for m in range(M):
        for j in range(J):
            Dis[m, :] += lfilter(PrimaryPath[m, j, :], [1.0], Ref[j, :])

    if backend == "numpy":
        controller: Any = MultiChannelFxLMS(
            wLen, SecondaryPath, sLen, Ref, Dis, J, K, M
        )
        controller = controller.McFxLMS_controller(muw)
        e_CMANC = controller.Err
    elif backend == "torch":
        from torch_multichannel import TorchMultiChannelFxLMS

        controller = TorchMultiChannelFxLMS(
            wLen,
            SecondaryPath,
            sLen,
            Ref,
            Dis,
            J,
            K,
            M,
            device=device,
            dtype=dtype,
        )
        controller = controller.McFxLMS_controller(muw)
        e_CMANC = controller.Err.detach().cpu().numpy()
    else:
        raise ValueError("backend must be 'numpy' or 'torch'.")

    if plot:
        plot_results(t, T, Fs, Dis, e_CMANC, backend)

    return {
        "t": t,
        "Ref": Ref,
        "Dis": Dis,
        "error": e_CMANC,
        "PrimaryPath": PrimaryPath,
        "SecondaryPath": SecondaryPath,
        "controller": controller,
    }


def plot_results(
    t: np.ndarray,
    T: float,
    Fs: int,
    Dis: np.ndarray,
    error: np.ndarray,
    backend: str,
) -> None:
    M = Dis.shape[0]
    sample_count = min(int(round(T * Fs)), t.size)
    rows = int(np.ceil(M / 2))
    figure, axes = plt.subplots(rows, 2, squeeze=False, sharex=True)
    for m, axis in enumerate(axes.flat):
        if m >= M:
            axis.set_visible(False)
            continue
        axis.plot(t[:sample_count], Dis[m, :sample_count], label="Disturbance")
        axis.plot(t[:sample_count], error[m, :sample_count], label="MCFxLMS")
        axis.set_title(f"({chr(ord('a') + m)}) Error {m + 1}")
        axis.set_xlabel("Time (seconds)")
        axis.set_ylabel("Amplitude")
        axis.grid(True)
        if m == 0:
            axis.legend()
    figure.suptitle(f"Fully connected MCFxLMS ({backend})")
    figure.tight_layout()
    plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration", type=float, default=30.0, help="Duration in s.")
    parser.add_argument("--references", type=int, default=6, help="J.")
    parser.add_argument("--sources", type=int, default=6, help="K.")
    parser.add_argument("--errors", type=int, default=6, help="M.")
    parser.add_argument("--backend", choices=("numpy", "torch"), default="numpy")
    parser.add_argument("--device", default="auto", help="auto, cpu, cuda, cuda:0...")
    parser.add_argument("--dtype", choices=("float32", "float64"), default="float64")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--no-plot", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_demo(
        T=args.duration,
        J=args.references,
        K=args.sources,
        M=args.errors,
        backend=args.backend,
        device=args.device,
        dtype=args.dtype,
        seed=args.seed,
        plot=not args.no_plot,
    )
