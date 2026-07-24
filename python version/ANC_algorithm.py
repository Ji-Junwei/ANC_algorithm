"""Single-channel FxLMS controller translated from ``ANC_algorithm.m``.

The class and public attribute names intentionally follow the MATLAB source so
that the two implementations can be compared line by line.  NumPy uses
zero-based indexing, but the signal flow and sign convention are unchanged:

    e(n) = d(n) - y(n)
    Wc(n + 1) = Wc(n) + mu * xf(n) * e(n)
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


class ANC_algorithm:
    """Single-channel time-domain FxLMS controller.

    Parameters use the same meaning as the MATLAB constructor:

    - ``wLen``: control-filter length.
    - ``sLen``: secondary-path FIR length.
    - ``SecondaryPath``: estimated/real secondary path used by the original
      example, with shape ``(sLen,)``.
    - ``Dis``: disturbance signal.
    - ``Ref``: reference signal.
    """

    def __init__(
        self,
        wLen: int,
        sLen: int,
        SecondaryPath: ArrayLike,
        Dis: ArrayLike,
        Ref: ArrayLike,
    ) -> None:
        if wLen <= 0 or sLen <= 0:
            raise ValueError("wLen and sLen must be positive integers.")

        secondary_path = np.asarray(SecondaryPath).reshape(-1)
        disturbance = np.asarray(Dis).reshape(-1)
        reference = np.asarray(Ref).reshape(-1)
        dtype = np.result_type(
            secondary_path.dtype,
            disturbance.dtype,
            reference.dtype,
            np.float64,
        )

        self.wlen = int(wLen)
        self.slen = int(sLen)
        self.SecP: NDArray = secondary_path.astype(dtype, copy=True)
        self.Dis: NDArray = disturbance.astype(dtype, copy=True)
        self.Ref: NDArray = reference.astype(dtype, copy=True)

        if self.SecP.size != self.slen:
            raise ValueError(
                f"SecondaryPath has {self.SecP.size} taps; expected sLen={self.slen}."
            )
        if self.Dis.size != self.Ref.size:
            raise ValueError("Dis and Ref must contain the same number of samples.")

        self.N = self.Dis.size
        self.Wc: NDArray = np.zeros(self.wlen, dtype=dtype)
        self.yc: NDArray = np.zeros(self.N, dtype=dtype)

    def ANC_FxLMS(self, muw: float) -> tuple[NDArray, "ANC_algorithm"]:
        """Run the FxLMS iteration and return ``(error, self)``.

        Returning ``self`` is redundant in Python, but preserves the MATLAB
        call pattern ``[e, obj] = ANC_FxLMS(obj, muw)``.
        """

        e = np.zeros(self.N, dtype=self.Wc.dtype)

        # MATLAB names and dimensions are retained below.
        xc = np.zeros(self.wlen, dtype=self.Wc.dtype)  # control-filter x buffer
        xs = np.zeros(self.slen, dtype=self.Wc.dtype)  # secondary-path x buffer
        xf = np.zeros(self.wlen, dtype=self.Wc.dtype)  # filtered-reference buffer
        ys = np.zeros(self.slen, dtype=self.Wc.dtype)  # secondary-path y buffer

        for i in range(self.N):
            # xc = [Ref(i); xc(1:end-1)]
            xc[1:] = xc[:-1]
            xc[0] = self.Ref[i]
            self.yc[i] = np.dot(self.Wc, xc)

            # Pass the control signal through the secondary path.
            ys[1:] = ys[:-1]
            ys[0] = self.yc[i]
            y = np.dot(ys, self.SecP)

            # Keep the MATLAB sign convention exactly: e = d - y.
            e[i] = self.Dis[i] - y

            # Filter the reference through the secondary-path estimate.
            xs[1:] = xs[:-1]
            xs[0] = self.Ref[i]
            fx = np.dot(xs, self.SecP)
            xf[1:] = xf[:-1]
            xf[0] = fx

            # FxLMS control-filter update.
            self.Wc += float(muw) * xf * e[i]

        return e, self
