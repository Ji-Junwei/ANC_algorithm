"""Single-reference multichannel FxLMS translated from ``McANC_SRMSE.m``.

Dimension order follows the MATLAB code:

- control filters ``Wc``: ``K x Lw``
- secondary paths ``SecP``: ``M x K x Ls``
- disturbance/error: ``M x N``
- reference: ``N``

The specialized 1x2x2, 1x4x4, and 1x6x6 MATLAB entry points are retained.
They call the arbitrary-channel implementation because the expanded MATLAB
equations are algebraically identical to the generic nested loops.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


class McANC_SRMSE:
    """Conventional 1-reference, K-source, M-error-sensor FxLMS."""

    def __init__(
        self,
        wLen: int,
        SecondaryPath: ArrayLike,
        sLen: int,
        Nums: int,
        Nume: int,
        Dis: ArrayLike,
        Ref: ArrayLike,
    ) -> None:
        if min(wLen, sLen, Nums, Nume) <= 0:
            raise ValueError("wLen, sLen, Nums, and Nume must be positive.")

        secondary_path = np.asarray(SecondaryPath)
        disturbance = np.asarray(Dis)
        reference = np.asarray(Ref).reshape(-1)
        dtype = np.result_type(
            secondary_path.dtype,
            disturbance.dtype,
            reference.dtype,
            np.float64,
        )

        self.wlen = int(wLen)
        self.slen = int(sLen)
        self.Nums = int(Nums)  # K: number of secondary sources
        self.Nume = int(Nume)  # M: number of error sensors
        self.Ref: NDArray = reference.astype(dtype, copy=True)
        self.Dis: NDArray = disturbance.astype(dtype, copy=True)
        self.SecP: NDArray = secondary_path.astype(dtype, copy=True)

        expected_secondary_shape = (self.Nume, self.Nums, self.slen)
        expected_disturbance_shape = (self.Nume, self.Ref.size)
        if self.SecP.shape != expected_secondary_shape:
            raise ValueError(
                f"SecondaryPath shape is {self.SecP.shape}; "
                f"expected {expected_secondary_shape} (M, K, Ls)."
            )
        if self.Dis.shape != expected_disturbance_shape:
            raise ValueError(
                f"Dis shape is {self.Dis.shape}; "
                f"expected {expected_disturbance_shape} (M, N)."
            )

        self.Wc: NDArray = np.zeros((self.Nums, self.wlen), dtype=dtype)
        self.yc: NDArray = np.zeros((self.Nums, self.Ref.size), dtype=dtype)

    def _require_square_channel_count(self, expected: int, method: str) -> None:
        if self.Nums != expected or self.Nume != expected:
            raise ValueError(
                f"{method} requires K=M={expected}; "
                f"received K={self.Nums}, M={self.Nume}."
            )

    def McFxLMS_SRMSE_122(self, muw: float) -> tuple[NDArray, "McANC_SRMSE"]:
        """MATLAB-compatible entry point for the 1 x 2 x 2 case."""

        self._require_square_channel_count(2, "McFxLMS_SRMSE_122")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_144(self, muw: float) -> tuple[NDArray, "McANC_SRMSE"]:
        """MATLAB-compatible entry point for the 1 x 4 x 4 case."""

        self._require_square_channel_count(4, "McFxLMS_SRMSE_144")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_166(self, muw: float) -> tuple[NDArray, "McANC_SRMSE"]:
        """MATLAB-compatible entry point for the 1 x 6 x 6 case."""

        self._require_square_channel_count(6, "McFxLMS_SRMSE_166")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_ANC(
        self, muw: float
    ) -> tuple[NDArray, "McANC_SRMSE"]:
        """Run arbitrary-channel SRMSE FxLMS using the MATLAB loop order."""

        N = self.Ref.size
        K = self.Nums
        M = self.Nume
        Lw = self.wlen
        Ls = self.slen

        e = np.zeros((M, N), dtype=self.Wc.dtype)
        xc = np.zeros(Lw, dtype=self.Wc.dtype)
        xs = np.zeros(Ls, dtype=self.Wc.dtype)
        ys = np.zeros((K, Ls), dtype=self.Wc.dtype)
        xf = np.zeros((K, M, Lw), dtype=self.Wc.dtype)

        for i in range(N):
            # Reference buffer for the K control filters.
            xc[1:] = xc[:-1]
            xc[0] = self.Ref[i]
            self.yc[:, i] = self.Wc @ xc

            # K control outputs pass through M x K secondary paths.
            ys[:, 1:] = ys[:, :-1]
            ys[:, 0] = self.yc[:, i]
            for mIdx in range(M):
                y_m = self.Wc.dtype.type(0)
                for kIdx in range(K):
                    y_m += np.dot(self.SecP[mIdx, kIdx, :], ys[kIdx, :])
                e[mIdx, i] = self.Dis[mIdx, i] - y_m

            # Update the K x M filtered-reference delay lines.
            xs[1:] = xs[:-1]
            xs[0] = self.Ref[i]
            for kIdx in range(K):
                for mIdx in range(M):
                    xfm = np.dot(xs, self.SecP[mIdx, kIdx, :])
                    xf[kIdx, mIdx, 1:] = xf[kIdx, mIdx, :-1]
                    xf[kIdx, mIdx, 0] = xfm

            # Sum the M error-sensor gradients for each control filter.
            for kIdx in range(K):
                grad_k = np.zeros(Lw, dtype=self.Wc.dtype)
                for mIdx in range(M):
                    grad_k += e[mIdx, i] * xf[kIdx, mIdx, :]
                self.Wc[kIdx, :] += float(muw) * grad_k

        return e, self
