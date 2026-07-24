"""1 x K x M MEFxLMS demonstration translated from ``MCANC_MEFxLMS.m``."""

from __future__ import annotations

import argparse
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lfilter

from McANC_SRMSE import McANC_SRMSE
from demo_utils import (
    awgn_measured,
    matlab_fir1_bandpass,
    moving_average,
    time_vector,
)


def _run_numpy(
    wLen: int,
    SecondaryPath: np.ndarray,
    sLen: int,
    K: int,
    M: int,
    Dis: np.ndarray,
    Ref: np.ndarray,
    muw: float,
) -> tuple[np.ndarray, McANC_SRMSE]:
    controller = McANC_SRMSE(wLen, SecondaryPath, sLen, K, M, Dis, Ref)
    specialized = {
        2: controller.McFxLMS_SRMSE_122,
        4: controller.McFxLMS_SRMSE_144,
        6: controller.McFxLMS_SRMSE_166,
    }
    method = specialized.get(K, controller.McFxLMS_SRMSE_ANC) if K == M else (
        controller.McFxLMS_SRMSE_ANC
    )
    return method(muw)


def _run_torch(
    wLen: int,
    SecondaryPath: np.ndarray,
    sLen: int,
    K: int,
    M: int,
    Dis: np.ndarray,
    Ref: np.ndarray,
    muw: float,
    device: str,
    dtype: str,
) -> tuple[np.ndarray, Any]:
    from torch_multichannel import TorchMcANC_SRMSE

    controller = TorchMcANC_SRMSE(
        wLen,
        SecondaryPath,
        sLen,
        K,
        M,
        Dis,
        Ref,
        device=device,
        dtype=dtype,
    )
    specialized = {
        2: controller.McFxLMS_SRMSE_122,
        4: controller.McFxLMS_SRMSE_144,
        6: controller.McFxLMS_SRMSE_166,
    }
    method = specialized.get(K, controller.McFxLMS_SRMSE_ANC) if K == M else (
        controller.McFxLMS_SRMSE_ANC
    )
    error, controller = method(muw)
    return error.detach().cpu().numpy(), controller


def run_demo(
    *,
    Fs: int = 16_000,
    T: float = 30.0,
    K: int = 6,
    M: int = 6,
    wLen: int = 512,
    sLen: int = 256,
    muw: float = 1e-6,
    snr_db: float = 40.0,
    backend: str = "numpy",
    device: str = "auto",
    dtype: str = "float64",
    seed: int | None = None,
    plot: bool = True,
) -> dict[str, Any]:
    """Run the MATLAB synthetic configuration with a selectable backend."""

    t = time_vector(Fs, T)
    N = t.size
    rng = np.random.default_rng(seed)

    Pri_path = matlab_fir1_bandpass(511, 80.0, 5000.0, Fs)
    Sec_path = matlab_fir1_bandpass(sLen - 1, 80.0, 5000.0, Fs)
    PrimaryPath = np.tile(Pri_path[None, :], (M, 1))
    SecondaryPath = np.tile(Sec_path[None, None, :], (M, K, 1))

    noise = rng.standard_normal(N)
    fil = matlab_fir1_bandpass(63, 100.0, 1000.0, Fs)
    Ref_clean = lfilter(fil, [1.0], noise)
    Dis = np.zeros((M, N), dtype=np.float64)
    for m in range(M):
        Dis[m, :] = lfilter(PrimaryPath[m, :], [1.0], Ref_clean)
    Ref = awgn_measured(Ref_clean, snr_db, rng)

    if backend == "numpy":
        e_CMANC, controller = _run_numpy(
            wLen, SecondaryPath, sLen, K, M, Dis, Ref, muw
        )
    elif backend == "torch":
        e_CMANC, controller = _run_torch(
            wLen,
            SecondaryPath,
            sLen,
            K,
            M,
            Dis,
            Ref,
            muw,
            device,
            dtype,
        )
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
    """Reproduce the MATLAB time-waveform and normalized-error plots."""

    M = Dis.shape[0]
    sample_count = min(int(round(T * Fs)), t.size)
    rows = int(np.ceil(M / 2))

    figure, axes = plt.subplots(rows, 2, squeeze=False, sharex=True)
    for m, axis in enumerate(axes.flat):
        if m >= M:
            axis.set_visible(False)
            continue
        axis.plot(t[:sample_count], Dis[m, :sample_count], label="Disturbance")
        axis.plot(t[:sample_count], error[m, :sample_count], label="MEFxLMS")
        axis.set_title(f"({chr(ord('a') + m)}) Error {m + 1}")
        axis.set_xlabel("Time (seconds)")
        axis.set_ylabel("Amplitude")
        axis.grid(True)
        if m == 0:
            axis.legend()
    figure.suptitle(f"MEFxLMS ({backend})")
    figure.tight_layout()

    figure, axes = plt.subplots(rows, 2, squeeze=False, sharex=True)
    eps = np.finfo(np.float64).tiny
    for m, axis in enumerate(axes.flat):
        if m >= M:
            axis.set_visible(False)
            continue
        dis_power = moving_average(Dis[m, :sample_count] ** 2, 2000)
        err_power = moving_average(error[m, :sample_count] ** 2, 2000)
        normalized_error = 10.0 * np.log10(
            np.maximum(err_power, eps) / np.maximum(dis_power, eps)
        )
        axis.plot(
            t[:sample_count],
            moving_average(normalized_error, 5000),
            label="MEFxLMS",
        )
        axis.set_title(f"({chr(ord('a') + m)}) Error {m + 1}")
        axis.set_xlabel("Time (seconds)")
        axis.set_ylabel("Normalized squared error (dB)")
        axis.set_ylim(top=5.0)
        axis.grid(True)
        if m == 0:
            axis.legend()
    figure.tight_layout()
    plt.show()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--duration", type=float, default=30.0, help="Duration in s.")
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
        K=args.sources,
        M=args.errors,
        backend=args.backend,
        device=args.device,
        dtype=args.dtype,
        seed=args.seed,
        plot=not args.no_plot,
    )
