"""Fully connected J x K x M FxLMS translated from MATLAB.

The implementation intentionally retains the MATLAB dimension order:

- ``Wc``: K x J x Lw
- ``SecP``: M x K x Ls
- ``Ref``: J x N
- ``Dis``, ``y``, ``Err``: M x N
- ``yc``: K x N
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


class MultiChannelFxLMS:
    """Conventional fully connected multichannel FxLMS controller."""

    def __init__(
        self,
        wLen: int,
        Secondarypath: ArrayLike,
        sLen: int,
        ref: ArrayLike,
        dis: ArrayLike,
        numref: int,
        numsec: int,
        numerr: int,
    ) -> None:
        if min(wLen, sLen, numref, numsec, numerr) <= 0:
            raise ValueError(
                "wLen, sLen, numref, numsec, and numerr must be positive."
            )

        secondary_path = np.asarray(Secondarypath)
        reference = np.asarray(ref)
        disturbance = np.asarray(dis)
        dtype = np.result_type(
            secondary_path.dtype,
            reference.dtype,
            disturbance.dtype,
            np.float64,
        )

        self.J = int(numref)
        self.K = int(numsec)
        self.M = int(numerr)
        self.wlen = int(wLen)
        self.slen = int(sLen)
        self.SecP: NDArray = secondary_path.astype(dtype, copy=True)
        self.Ref: NDArray = reference.astype(dtype, copy=True)
        self.Dis: NDArray = disturbance.astype(dtype, copy=True)

        if self.Ref.ndim != 2 or self.Ref.shape[0] != self.J:
            raise ValueError(
                f"ref shape is {self.Ref.shape}; expected (J, N) with J={self.J}."
            )
        N = self.Ref.shape[1]
        expected_secondary_shape = (self.M, self.K, self.slen)
        expected_disturbance_shape = (self.M, N)
        if self.SecP.shape != expected_secondary_shape:
            raise ValueError(
                f"Secondarypath shape is {self.SecP.shape}; "
                f"expected {expected_secondary_shape} (M, K, Ls)."
            )
        if self.Dis.shape != expected_disturbance_shape:
            raise ValueError(
                f"dis shape is {self.Dis.shape}; "
                f"expected {expected_disturbance_shape} (M, N)."
            )

        self.Wc: NDArray = np.zeros(
            (self.K, self.J, self.wlen), dtype=dtype
        )
        self.yc: NDArray = np.zeros((self.K, N), dtype=dtype)
        self.y: NDArray = np.zeros((self.M, N), dtype=dtype)
        self.Err: NDArray = np.zeros((self.M, N), dtype=dtype)

    def McFxLMS_controller(self, muw: float) -> "MultiChannelFxLMS":
        """Run the fully connected FxLMS iteration."""

        N = self.Ref.shape[1]

        xc = np.zeros((self.J, self.wlen), dtype=self.Wc.dtype)
        ys = np.zeros((self.K, self.slen), dtype=self.Wc.dtype)
        xf = np.zeros((self.J, self.slen), dtype=self.Wc.dtype)
        Xsf = np.zeros(
            (self.J, self.K, self.M, self.wlen), dtype=self.Wc.dtype
        )

        for n in range(N):
            # MATLAB: xc = [Ref(:,n), xc(:,1:end-1)]
            xc[:, 1:] = xc[:, :-1]
            xc[:, 0] = self.Ref[:, n]

            # K control signals, each formed from all J references.
            for kk in range(self.K):
                self.yc[kk, n] = np.sum(self.Wc[kk, :, :] * xc)

            # M anti-noise signals through the M x K secondary paths.
            ys[:, 1:] = ys[:, :-1]
            ys[:, 0] = self.yc[:, n]
            for mm in range(self.M):
                self.y[mm, n] = np.sum(self.SecP[mm, :, :] * ys)

            self.Err[:, n] = self.Dis[:, n] - self.y[:, n]

            # Filter every reference through every M x K secondary path.
            xf[:, 1:] = xf[:, :-1]
            xf[:, 0] = self.Ref[:, n]
            res = np.zeros((self.J, self.K, self.M), dtype=self.Wc.dtype)
            for jj in range(self.J):
                for kk in range(self.K):
                    for mm in range(self.M):
                        res[jj, kk, mm] = np.dot(
                            xf[jj, :], self.SecP[mm, kk, :]
                        )

            Xsf[:, :, :, 1:] = Xsf[:, :, :, :-1]
            Xsf[:, :, :, 0] = res

            # MATLAB Xsf1 order: K x J x Lw x M.
            delta = np.zeros_like(self.Wc)
            for mm in range(self.M):
                delta += self.Err[mm, n] * np.transpose(
                    Xsf[:, :, mm, :], (1, 0, 2)
                )
            self.Wc += float(muw) * delta

        return self
