# ANC_algorithm
Conventional ANC algorithms including single-channel FxLMS and multichannel
FxLMS, with parallel MATLAB and Python implementations.

## Repository structure

- `MATLAB version/`: original MATLAB implementation.
- `python version/`: MATLAB-oriented NumPy/SciPy translation plus optional
  PyTorch/CUDA multichannel extensions.

## Algorithms

- `ANC_algorithm`: single-channel FxLMS.
- `McANC_SRMSE`: multichannel ANC (1 x K x M) with one reference, K secondary
  sources, and M error sensors; MEFxLMS/SRMSE.
- `MultiChannelFxLMS`: fully connected multichannel FxLMS with arbitrary
  J references, K secondary sources, and M error sensors.

See [`python version/README.md`](python%20version/README.md) for Python
dependencies, file mapping, run commands, tensor dimensions, tests, and CUDA
notes.
