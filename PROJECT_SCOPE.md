# Project Scope

## Purpose

`ANC_algorithm` is the public reference repository for conventional active noise control algorithms.

It provides compact, readable MATLAB implementations together with closely matched Python translations for studying, validating, and extending single-channel and conventional multichannel FxLMS algorithms.

## Implementation authority

The repository contains two parallel implementation tracks:

```text
MATLAB version/
python version/
```

Their roles are:

- `MATLAB version/` contains the original algorithm implementations and is the primary mathematical and signal-flow reference.
- The NumPy/SciPy code in `python version/` is the canonical Python translation and should remain comparable to the MATLAB implementation.
- `python version/torch_multichannel.py` is an optional PyTorch CPU/CUDA extension for the multichannel algorithms. It is not a replacement for the NumPy reference implementation.

When MATLAB and Python behaviour differ, first determine whether the difference comes from indexing, array orientation, numeric precision, or a library-level implementation detail before changing the adaptive-filter equations.

## Algorithms in scope

### Single-channel FxLMS

Primary files:

```text
MATLAB version/ANC_algorithm.m
MATLAB version/ANC.m
python version/ANC_algorithm.py
python version/ANC.py
```

Core dimensions:

```text
Wc:   (Lw,)
SecP: (Ls,)
```

### Single-reference multichannel SRMSE / MEFxLMS

Primary files:

```text
MATLAB version/McANC_SRMSE.m
MATLAB version/MCANC_MEFxLMS.m
python version/McANC_SRMSE.py
python version/MCANC_MEFxLMS.py
```

System structure:

```text
1 reference × K secondary sources × M error sensors
```

Core dimensions:

```text
Wc:   (K, Lw)
SecP: (M, K, Ls)
error:(M, N)
```

The specialized `122`, `144`, and `166` entry points are retained for compatibility. Their Python implementations may call the generic channel-loop implementation when the equations are algebraically equivalent.

### Fully connected multichannel FxLMS

Primary files:

```text
MATLAB version/MultiChannelFxLMS.m
MATLAB version/MCFxLMS_tst.m
python version/MultiChannelFxLMS.py
python version/MCFxLMS_tst.py
```

System structure:

```text
J references × K secondary sources × M error sensors
```

Core dimensions:

```text
Wc:       (K, J, Lw)
SecP:     (M, K, Ls)
reference:(J, N)
error:    (M, N)
```

## Core algorithm invariants

Unless a task explicitly changes the mathematical definition, preserve:

```text
e(n) = d(n) - y(n)
W(n + 1) = W(n) + mu * gradient(n)
```

Also preserve:

- sample-by-sample recursive adaptation;
- the existing reference, secondary-source, error-sensor, and tap dimension order;
- the secondary-path filtering used to form the filtered reference;
- the public MATLAB-oriented class, property, and method naming in the Python translations where practical;
- double precision when comparing Python results with MATLAB;
- the distinction between the physical secondary path and its estimate when an implementation explicitly introduces both.

FxLMS is recursive in time because each weight update affects later samples. GPU implementations may vectorize channel and tap calculations, but they must not silently replace the time-recursive algorithm with an unrelated batch optimization.

## In scope

- Conventional single-channel FxLMS.
- Conventional centralized multichannel FxLMS.
- MATLAB examples and Python equivalents.
- NumPy/SciPy translations that preserve the original signal flow.
- Optional PyTorch CPU/CUDA acceleration for multichannel channel/tap calculations.
- Dimension, dtype, input, and device validation.
- MATLAB-to-Python numerical parity tests.
- Clear examples, run instructions, and algorithm documentation.
- Additional conventional ANC reference algorithms when their scope and mathematical definition are documented.

## Out of scope

- IMC feedback ANC and other private or application-specific feedback-controller development.
- Output-constrained MOV-FxLMS research maintained as a separate project.
- Distributed ANC communication, network, federation, or coprocessor algorithms maintained in separate repositories.
- Silent incorporation of paper-specific algorithms into a general reference implementation.
- Fixed-point, Q-format, firmware, hardware-driver, or production DSP deployment unless explicitly introduced as a separate extension.
- Private acoustic measurements, credentials, machine-specific paths, or generated experiment outputs.
- Claims that an optional GPU implementation is faster without a documented benchmark for the relevant dimensions and hardware.

## MATLAB and Python synchronization

For changes to an existing translated algorithm:

1. Identify the original MATLAB file and the corresponding Python file.
2. State whether the change is a bug fix, parity correction, validation improvement, optimization, or intentional algorithm extension.
3. Preserve line-by-line comparability where it does not conflict with correctness or Python safety.
4. Update both implementations when the mathematical algorithm itself changes, unless the task is explicitly language-specific.
5. Keep Python-only validation, typing, or device handling separate from the adaptive-filter equations.
6. Document intentional differences between MATLAB, NumPy, and PyTorch.

## Numerical and library conventions

- MATLAB uses double precision by default; use NumPy `float64` and PyTorch `float64` for parity checks.
- PyTorch `float32` may be used for performance experiments only after numerical differences are checked.
- MATLAB `fir1` and SciPy `firwin` can produce small coefficient differences despite matched order, window, and scaling choices.
- Such library-level differences must not be misreported as changes to the FxLMS equations.
- Array shape validation should fail clearly rather than silently reshape incompatible multichannel data.

## Validation expectations

Changes should use the strongest applicable checks, including:

```bash
cd "python version"
python -m unittest discover -s tests -v
```

As applicable, also verify:

- MATLAB and NumPy output parity on deterministic short signals;
- NumPy and PyTorch parity in `float64`;
- expected tensor and array dimensions;
- sign convention and weight-update order;
- finite controller weights, outputs, and error signals;
- CPU/CUDA device behaviour;
- example scripts with short duration and plotting disabled;
- documentation and file mapping after files are added, renamed, or removed.

Report every modified, added, renamed, or deleted file and clearly separate static validation from results requiring MATLAB, CUDA, or local runtime data.
