# Python version

This folder is a MATLAB-oriented Python translation of the algorithms in
`../MATLAB version/`. The NumPy implementations are the canonical translations;
the PyTorch/CUDA implementation is an optional multichannel extension.

## MATLAB-to-Python mapping

| MATLAB file | Python file | Purpose |
|---|---|---|
| `ANC_algorithm.m` | `ANC_algorithm.py` | Single-channel time-domain FxLMS class |
| `ANC.m` | `ANC.py` | Single-channel synthetic example |
| `McANC_SRMSE.m` | `McANC_SRMSE.py` | Single-reference, K-source, M-error FxLMS |
| `MCANC_MEFxLMS.m` | `MCANC_MEFxLMS.py` | 1 x K x M example |
| `MultiChannelFxLMS.m` | `MultiChannelFxLMS.py` | Fully connected J x K x M FxLMS |
| `MCFxLMS_tst.m` | `MCFxLMS_tst.py` | Fully connected multichannel example |
| — | `torch_multichannel.py` | Optional PyTorch CPU/CUDA extension |

The Python class names, public properties, method names, dimension order, and
sample-by-sample update order intentionally remain close to MATLAB. Python uses
zero-based indexing, and buffer shifts are performed in place to avoid creating
a new array at every sample.

## Preserved conventions

- Single channel: `Wc (Lw,)`, `SecP (Ls,)`.
- SRMSE: `Wc (K, Lw)`, `SecP (M, K, Ls)`, error `(M, N)`.
- Fully connected: `Wc (K, J, Lw)`, `SecP (M, K, Ls)`,
  reference `(J, N)`, error `(M, N)`.
- The original sign convention is preserved:

  ```text
  e(n) = d(n) - y(n)
  W(n + 1) = W(n) + mu * gradient(n)
  ```

- The specialized `McFxLMS_SRMSE_122`, `144`, and `166` entry points are
  retained. They call the generic channel-loop implementation because the
  expanded MATLAB equations are algebraically identical.

## Installation

Core NumPy/SciPy version:

```bash
python -m pip install -r requirements.txt
```

Optional PyTorch version:

```bash
python -m pip install -r requirements-torch.txt
```

For CUDA, install the PyTorch build that matches the local NVIDIA driver and
CUDA runtime. The code itself selects CUDA with `--device auto` when
`torch.cuda.is_available()` is true.

## Run

Run commands from this folder:

```bash
python ANC.py
python MCANC_MEFxLMS.py --backend numpy
python MCFxLMS_tst.py --backend numpy
```

Optional PyTorch/CUDA:

```bash
python MCANC_MEFxLMS.py --backend torch --device auto --dtype float64
python MCFxLMS_tst.py --backend torch --device cuda --dtype float32
```

The defaults preserve the MATLAB 30-second examples. For a quick check, reduce
the duration:

```bash
python MCANC_MEFxLMS.py --duration 0.1 --no-plot
```

## GPU scope and limitation

FxLMS updates are recursive in time: the weight update at sample `n` changes
the output at sample `n+1`. The PyTorch implementation therefore retains the
sample loop and accelerates the channel/tap calculations inside each sample.
CUDA is most useful for larger J/K/M and longer filters. For small systems,
kernel-launch and device overhead can make the NumPy CPU version faster.

Use `float64` when comparing against MATLAB because MATLAB uses double precision
by default. Use `float32` only when speed/memory is more important and the
resulting numerical difference has been checked.

## Tests

```bash
python -m unittest discover -s tests -v
```

PyTorch parity tests are skipped automatically if PyTorch is not installed.

## Known library-level difference

MATLAB `fir1` is represented by SciPy `firwin` with the same order-plus-one
length, Hamming window, and passband scaling. Small FIR coefficient differences
may remain between MATLAB and SciPy versions. The adaptive-filter equations are
not changed.
