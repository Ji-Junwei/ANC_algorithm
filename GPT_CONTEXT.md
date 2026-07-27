# GPT Context

## Repository identity

This is the public `ANC-Research/ANC_algorithm` repository.

It contains conventional single-channel and centralized multichannel active noise control reference algorithms with parallel MATLAB and Python implementations. Do not treat it as an IMC-feedback, output-constrained, distributed-communication, or application-specific research repository.

Because this repository is public, do not add private repository names, unpublished methods, private data descriptions, internal results, credentials, or future confidential plans.

## Required reading before making changes

Before modifying code, read the latest version of:

1. `PROJECT_SCOPE.md`;
2. the root `README.md`;
3. `python version/README.md` for Python-related work;
4. the target MATLAB and/or Python implementation;
5. the corresponding example script;
6. relevant tests and requirements files.

Always use the latest requested branch. Do not rely only on previous conversations, cached code, or similarly named algorithms in another repository.

## Implementation authority

- `MATLAB version/` contains the original reference implementations.
- The NumPy/SciPy files in `python version/` are the canonical Python translations.
- `python version/torch_multichannel.py` is an optional acceleration extension, not the authoritative replacement for NumPy.
- When changing an adaptive-filter equation, inspect the corresponding MATLAB and Python implementations together.
- Keep unrelated algorithms and examples read-only.

## Task isolation

At the start of each task, identify:

- the exact algorithm: single-channel FxLMS, 1×K×M SRMSE/MEFxLMS, or J×K×M fully connected FxLMS;
- the target language: MATLAB, NumPy, PyTorch, or multiple implementations;
- the authoritative source and corresponding translation;
- the expected input/output dimensions;
- whether the task changes equations, validation, examples, performance, tests, or documentation.

Do not assume that a change requested for one multichannel formulation applies to the other.

## Core invariants

Preserve these unless the user explicitly changes them:

```text
e(n) = d(n) - y(n)
W(n + 1) = W(n) + mu * gradient(n)
```

Also preserve:

- sample-by-sample recursive FxLMS adaptation;
- current dimension order for references, secondary sources, error sensors, and FIR taps;
- filtered-reference construction through the secondary-path estimate;
- MATLAB-oriented Python names where they support direct comparison;
- NumPy `float64` and PyTorch `float64` for MATLAB parity;
- specialized SRMSE compatibility entry points (`122`, `144`, `166`).

Do not convert the algorithm into a batch optimizer merely to improve GPU utilization.

## MATLAB-to-Python changes

For translated algorithms:

1. compare the corresponding MATLAB and Python files before editing;
2. distinguish zero-based indexing and array-layout differences from mathematical differences;
3. preserve public class, property, and method mappings where practical;
4. update tests when dimensions, return values, validation, or numeric behaviour change;
5. document intentional deviations from MATLAB;
6. do not change both implementations solely for stylistic consistency.

A Python-only safety improvement, such as shape validation or device selection, does not automatically require a MATLAB change. A mathematical algorithm change normally requires both implementations to be reviewed and updated.

## PyTorch and CUDA changes

When modifying the optional PyTorch implementation:

- keep NumPy as the parity reference;
- preserve manual, sample-recursive weight updates;
- do not use autograd unless a new algorithm explicitly requires it;
- verify CPU behaviour before CUDA-specific claims;
- test `float64` parity and clearly label `float32` differences;
- avoid claiming speedups without dimensions, filter lengths, hardware, dtype, and timing methodology;
- retain clear errors when CUDA is requested but unavailable.

## Change discipline

- Make the smallest coherent change that satisfies the request.
- Do not add paper-specific, private, or unrelated ANC algorithms without explicit scope approval.
- Do not copy code or parameters from other repositories without documenting the source and version.
- Do not commit generated figures, caches, local virtual environments, private data, or machine-specific absolute paths.
- Report every modified, added, renamed, or deleted file.
- Explain equation changes and map them to exact implementation methods.
- Update the root README or `python version/README.md` when file mappings, algorithms, run commands, dependencies, dimensions, or supported backends change.
- Add a concise entry to `CHANGELOG.md` for significant user-visible changes.

## Verification

For Python changes, normally run from `python version/`:

```bash
python -m unittest discover -s tests -v
```

Also use short deterministic checks where appropriate:

```bash
python ANC.py
python MCANC_MEFxLMS.py --duration 0.1 --no-plot --backend numpy
python MCFxLMS_tst.py --duration 0.1 --no-plot --backend numpy
```

For PyTorch work, compare NumPy and PyTorch in `float64` on the same short inputs. CUDA-dependent results must be reported separately when CUDA is unavailable.

For MATLAB changes, identify the scripts that must be run locally and clearly separate code inspection from MATLAB execution results.

## Response requirements

For every completed modification, provide:

- files changed;
- the purpose of each change;
- the affected algorithm and dimensions;
- the equation or implementation mapping;
- MATLAB/Python/PyTorch compatibility impact;
- verification performed and results;
- limitations or local runtime checks still required.
