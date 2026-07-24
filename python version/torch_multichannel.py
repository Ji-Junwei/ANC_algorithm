"""Optional PyTorch/CUDA extensions for the two multichannel controllers.

The NumPy classes remain the canonical MATLAB-style translations.  These
PyTorch classes preserve the same state variables, dimensions, sample-by-
sample adaptation, and sign convention while vectorizing channel operations
with ``torch.einsum``.

FxLMS is recursive over time because every weight update affects the next
sample.  Therefore, the time loop cannot simply be processed as one batch.
CUDA is most likely to help when J/K/M and the FIR lengths are large; for small
systems, per-sample GPU kernel-launch overhead can make NumPy faster.
"""

from __future__ import annotations

from typing import Any

try:
    import torch
except ImportError as exc:  # pragma: no cover - depends on optional package
    raise ImportError(
        "torch_multichannel.py requires PyTorch. Install a PyTorch build that "
        "matches your CPU/CUDA environment before importing this module."
    ) from exc


def resolve_device(device: str | torch.device = "auto") -> torch.device:
    """Resolve ``auto`` to CUDA when available, otherwise CPU."""

    if isinstance(device, torch.device):
        resolved = device
    elif device == "auto":
        resolved = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    else:
        resolved = torch.device(device)

    if resolved.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError(
            f"CUDA device {resolved} was requested, but torch.cuda.is_available() "
            "is False."
        )
    return resolved


def resolve_dtype(dtype: str | torch.dtype = "float64") -> torch.dtype:
    """Map a readable dtype option to a PyTorch dtype.

    ``float64`` is the default because MATLAB and the NumPy translation use
    double precision.  ``float32`` is usually faster on consumer GPUs.
    """

    if isinstance(dtype, torch.dtype):
        if dtype not in (torch.float32, torch.float64):
            raise ValueError("Only torch.float32 and torch.float64 are supported.")
        return dtype

    options = {
        "float32": torch.float32,
        "float64": torch.float64,
        "single": torch.float32,
        "double": torch.float64,
    }
    try:
        return options[dtype.lower()]
    except (AttributeError, KeyError) as exc:
        raise ValueError("dtype must be float32/single or float64/double.") from exc


class TorchMcANC_SRMSE:
    """PyTorch/CUDA version of the 1 x K x M SRMSE FxLMS controller."""

    def __init__(
        self,
        wLen: int,
        SecondaryPath: Any,
        sLen: int,
        Nums: int,
        Nume: int,
        Dis: Any,
        Ref: Any,
        *,
        device: str | torch.device = "auto",
        dtype: str | torch.dtype = "float64",
    ) -> None:
        if min(wLen, sLen, Nums, Nume) <= 0:
            raise ValueError("wLen, sLen, Nums, and Nume must be positive.")

        self.device = resolve_device(device)
        self.dtype = resolve_dtype(dtype)
        self.wlen = int(wLen)
        self.slen = int(sLen)
        self.Nums = int(Nums)
        self.Nume = int(Nume)

        self.SecP = torch.as_tensor(
            SecondaryPath, dtype=self.dtype, device=self.device
        ).clone()
        self.Dis = torch.as_tensor(
            Dis, dtype=self.dtype, device=self.device
        ).clone()
        self.Ref = torch.as_tensor(
            Ref, dtype=self.dtype, device=self.device
        ).reshape(-1).clone()

        expected_secondary_shape = (self.Nume, self.Nums, self.slen)
        expected_disturbance_shape = (self.Nume, self.Ref.numel())
        if tuple(self.SecP.shape) != expected_secondary_shape:
            raise ValueError(
                f"SecondaryPath shape is {tuple(self.SecP.shape)}; "
                f"expected {expected_secondary_shape} (M, K, Ls)."
            )
        if tuple(self.Dis.shape) != expected_disturbance_shape:
            raise ValueError(
                f"Dis shape is {tuple(self.Dis.shape)}; "
                f"expected {expected_disturbance_shape} (M, N)."
            )

        self.Wc = torch.zeros(
            (self.Nums, self.wlen), dtype=self.dtype, device=self.device
        )
        self.yc = torch.zeros(
            (self.Nums, self.Ref.numel()),
            dtype=self.dtype,
            device=self.device,
        )

    def _require_square_channel_count(self, expected: int, method: str) -> None:
        if self.Nums != expected or self.Nume != expected:
            raise ValueError(
                f"{method} requires K=M={expected}; "
                f"received K={self.Nums}, M={self.Nume}."
            )

    def McFxLMS_SRMSE_122(
        self, muw: float
    ) -> tuple[torch.Tensor, "TorchMcANC_SRMSE"]:
        self._require_square_channel_count(2, "McFxLMS_SRMSE_122")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_144(
        self, muw: float
    ) -> tuple[torch.Tensor, "TorchMcANC_SRMSE"]:
        self._require_square_channel_count(4, "McFxLMS_SRMSE_144")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_166(
        self, muw: float
    ) -> tuple[torch.Tensor, "TorchMcANC_SRMSE"]:
        self._require_square_channel_count(6, "McFxLMS_SRMSE_166")
        return self.McFxLMS_SRMSE_ANC(muw)

    def McFxLMS_SRMSE_ANC(
        self, muw: float
    ) -> tuple[torch.Tensor, "TorchMcANC_SRMSE"]:
        """Run SRMSE FxLMS on the selected CPU or CUDA device."""

        N = self.Ref.numel()
        K = self.Nums
        M = self.Nume
        Lw = self.wlen
        Ls = self.slen
        mu = torch.as_tensor(muw, dtype=self.dtype, device=self.device)

        e = torch.zeros((M, N), dtype=self.dtype, device=self.device)
        xc = torch.zeros(Lw, dtype=self.dtype, device=self.device)
        xs = torch.zeros(Ls, dtype=self.dtype, device=self.device)
        ys = torch.zeros((K, Ls), dtype=self.dtype, device=self.device)
        xf = torch.zeros((K, M, Lw), dtype=self.dtype, device=self.device)

        # The update is explicitly manual; autograd is neither needed nor used.
        with torch.no_grad():
            for i in range(N):
                xc[1:].copy_(xc[:-1].clone())
                xc[0] = self.Ref[i]
                self.yc[:, i].copy_(self.Wc @ xc)

                ys[:, 1:].copy_(ys[:, :-1].clone())
                ys[:, 0].copy_(self.yc[:, i])
                anti_noise = torch.einsum("mkl,kl->m", self.SecP, ys)
                e[:, i].copy_(self.Dis[:, i] - anti_noise)

                xs[1:].copy_(xs[:-1].clone())
                xs[0] = self.Ref[i]
                current_xf = torch.einsum("l,mkl->km", xs, self.SecP)
                xf[:, :, 1:].copy_(xf[:, :, :-1].clone())
                xf[:, :, 0].copy_(current_xf)

                delta = torch.einsum("m,kml->kl", e[:, i], xf)
                self.Wc.add_(mu * delta)

        return e, self

    def cpu_numpy(self) -> dict[str, Any]:
        """Return detached NumPy copies of the main controller state."""

        return {
            "Wc": self.Wc.detach().cpu().numpy(),
            "yc": self.yc.detach().cpu().numpy(),
        }


class TorchMultiChannelFxLMS:
    """PyTorch/CUDA version of fully connected J x K x M FxLMS."""

    def __init__(
        self,
        wLen: int,
        Secondarypath: Any,
        sLen: int,
        ref: Any,
        dis: Any,
        numref: int,
        numsec: int,
        numerr: int,
        *,
        device: str | torch.device = "auto",
        dtype: str | torch.dtype = "float64",
    ) -> None:
        if min(wLen, sLen, numref, numsec, numerr) <= 0:
            raise ValueError(
                "wLen, sLen, numref, numsec, and numerr must be positive."
            )

        self.device = resolve_device(device)
        self.dtype = resolve_dtype(dtype)
        self.J = int(numref)
        self.K = int(numsec)
        self.M = int(numerr)
        self.wlen = int(wLen)
        self.slen = int(sLen)

        self.SecP = torch.as_tensor(
            Secondarypath, dtype=self.dtype, device=self.device
        ).clone()
        self.Ref = torch.as_tensor(
            ref, dtype=self.dtype, device=self.device
        ).clone()
        self.Dis = torch.as_tensor(
            dis, dtype=self.dtype, device=self.device
        ).clone()

        if self.Ref.ndim != 2 or self.Ref.shape[0] != self.J:
            raise ValueError(
                f"ref shape is {tuple(self.Ref.shape)}; "
                f"expected (J, N) with J={self.J}."
            )
        N = self.Ref.shape[1]
        expected_secondary_shape = (self.M, self.K, self.slen)
        expected_disturbance_shape = (self.M, N)
        if tuple(self.SecP.shape) != expected_secondary_shape:
            raise ValueError(
                f"Secondarypath shape is {tuple(self.SecP.shape)}; "
                f"expected {expected_secondary_shape} (M, K, Ls)."
            )
        if tuple(self.Dis.shape) != expected_disturbance_shape:
            raise ValueError(
                f"dis shape is {tuple(self.Dis.shape)}; "
                f"expected {expected_disturbance_shape} (M, N)."
            )

        self.Wc = torch.zeros(
            (self.K, self.J, self.wlen),
            dtype=self.dtype,
            device=self.device,
        )
        self.yc = torch.zeros(
            (self.K, N), dtype=self.dtype, device=self.device
        )
        self.y = torch.zeros(
            (self.M, N), dtype=self.dtype, device=self.device
        )
        self.Err = torch.zeros(
            (self.M, N), dtype=self.dtype, device=self.device
        )

    def McFxLMS_controller(self, muw: float) -> "TorchMultiChannelFxLMS":
        """Run fully connected FxLMS on the selected device."""

        N = self.Ref.shape[1]
        mu = torch.as_tensor(muw, dtype=self.dtype, device=self.device)
        xc = torch.zeros(
            (self.J, self.wlen), dtype=self.dtype, device=self.device
        )
        ys = torch.zeros(
            (self.K, self.slen), dtype=self.dtype, device=self.device
        )
        xf = torch.zeros(
            (self.J, self.slen), dtype=self.dtype, device=self.device
        )
        Xsf = torch.zeros(
            (self.J, self.K, self.M, self.wlen),
            dtype=self.dtype,
            device=self.device,
        )

        with torch.no_grad():
            for n in range(N):
                xc[:, 1:].copy_(xc[:, :-1].clone())
                xc[:, 0].copy_(self.Ref[:, n])
                self.yc[:, n].copy_(
                    torch.einsum("kjl,jl->k", self.Wc, xc)
                )

                ys[:, 1:].copy_(ys[:, :-1].clone())
                ys[:, 0].copy_(self.yc[:, n])
                self.y[:, n].copy_(
                    torch.einsum("mkl,kl->m", self.SecP, ys)
                )
                self.Err[:, n].copy_(self.Dis[:, n] - self.y[:, n])

                xf[:, 1:].copy_(xf[:, :-1].clone())
                xf[:, 0].copy_(self.Ref[:, n])
                current_xsf = torch.einsum("jl,mkl->jkm", xf, self.SecP)
                Xsf[:, :, :, 1:].copy_(Xsf[:, :, :, :-1].clone())
                Xsf[:, :, :, 0].copy_(current_xsf)

                delta = torch.einsum("m,jkml->kjl", self.Err[:, n], Xsf)
                self.Wc.add_(mu * delta)

        return self

    def cpu_numpy(self) -> dict[str, Any]:
        """Return detached NumPy copies of the main controller state."""

        return {
            "Wc": self.Wc.detach().cpu().numpy(),
            "yc": self.yc.detach().cpu().numpy(),
            "y": self.y.detach().cpu().numpy(),
            "Err": self.Err.detach().cpu().numpy(),
        }
