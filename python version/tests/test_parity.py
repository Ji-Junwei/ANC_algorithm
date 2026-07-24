"""Small deterministic checks for the MATLAB-style Python translations."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

PYTHON_VERSION_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PYTHON_VERSION_DIR))

from ANC_algorithm import ANC_algorithm
from McANC_SRMSE import McANC_SRMSE
from MultiChannelFxLMS import MultiChannelFxLMS


def literal_single_channel(
    wlen: int,
    secp: np.ndarray,
    disturbance: np.ndarray,
    reference: np.ndarray,
    mu: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Independent literal loop matching the equations in ANC_algorithm.m."""

    wc = np.zeros(wlen, dtype=np.float64)
    yc = np.zeros(reference.size, dtype=np.float64)
    error = np.zeros(reference.size, dtype=np.float64)
    xc = np.zeros(wlen, dtype=np.float64)
    xs = np.zeros(secp.size, dtype=np.float64)
    xf = np.zeros(wlen, dtype=np.float64)
    ys = np.zeros(secp.size, dtype=np.float64)

    for n in range(reference.size):
        xc = np.concatenate(([reference[n]], xc[:-1]))
        yc[n] = wc @ xc
        ys = np.concatenate(([yc[n]], ys[:-1]))
        error[n] = disturbance[n] - ys @ secp
        xs = np.concatenate(([reference[n]], xs[:-1]))
        xf = np.concatenate(([xs @ secp], xf[:-1]))
        wc = wc + mu * xf * error[n]

    return error, wc, yc


class TestSingleChannel(unittest.TestCase):
    def test_matches_literal_matlab_loop(self) -> None:
        reference = np.array([1.0, -0.25, 0.5, 0.125, -0.75])
        disturbance = np.array([0.4, -0.1, 0.2, -0.3, 0.15])
        secondary = np.array([0.7, -0.2])
        expected_e, expected_w, expected_yc = literal_single_channel(
            3, secondary, disturbance, reference, 0.05
        )

        controller = ANC_algorithm(3, 2, secondary, disturbance, reference)
        error, returned = controller.ANC_FxLMS(0.05)

        self.assertIs(returned, controller)
        np.testing.assert_allclose(error, expected_e, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(controller.Wc, expected_w, rtol=0.0, atol=1e-14)
        np.testing.assert_allclose(
            controller.yc, expected_yc, rtol=0.0, atol=1e-14
        )

    def test_identity_path_converges(self) -> None:
        rng = np.random.default_rng(1)
        reference = rng.standard_normal(500)
        controller = ANC_algorithm(
            1, 1, np.array([1.0]), reference.copy(), reference
        )
        error, _ = controller.ANC_FxLMS(0.05)

        early_power = np.mean(error[:50] ** 2)
        late_power = np.mean(error[-50:] ** 2)
        self.assertLess(late_power, early_power * 1e-4)


class TestSingleReferenceMultichannel(unittest.TestCase):
    def test_first_sample_axis_and_gradient(self) -> None:
        reference = np.array([2.0])
        disturbance = np.array([[1.0], [3.0]])
        secondary = np.array(
            [
                [[0.5], [-0.25]],
                [[0.1], [0.75]],
            ]
        )
        mu = 0.2
        controller = McANC_SRMSE(
            3, secondary, 1, 2, 2, disturbance, reference
        )
        error, controller = controller.McFxLMS_SRMSE_122(mu)

        np.testing.assert_allclose(error[:, 0], disturbance[:, 0])
        expected_tap_zero = np.empty(2)
        for k in range(2):
            expected_tap_zero[k] = (
                mu
                * reference[0]
                * np.sum(disturbance[:, 0] * secondary[:, k, 0])
            )
        np.testing.assert_allclose(controller.Wc[:, 0], expected_tap_zero)
        np.testing.assert_allclose(controller.Wc[:, 1:], 0.0)

    def test_specialized_122_matches_generic(self) -> None:
        rng = np.random.default_rng(10)
        reference = rng.standard_normal(12)
        disturbance = rng.standard_normal((2, 12))
        secondary = rng.standard_normal((2, 2, 3)) * 0.1

        specialized = McANC_SRMSE(
            4, secondary, 3, 2, 2, disturbance, reference
        )
        generic = McANC_SRMSE(4, secondary, 3, 2, 2, disturbance, reference)
        e_specialized, specialized = specialized.McFxLMS_SRMSE_122(1e-3)
        e_generic, generic = generic.McFxLMS_SRMSE_ANC(1e-3)

        np.testing.assert_allclose(e_specialized, e_generic)
        np.testing.assert_allclose(specialized.Wc, generic.Wc)
        np.testing.assert_allclose(specialized.yc, generic.yc)

    def test_one_by_one_identity_path_converges(self) -> None:
        rng = np.random.default_rng(2)
        reference = rng.standard_normal(500)
        disturbance = reference[None, :].copy()
        controller = McANC_SRMSE(
            1,
            np.ones((1, 1, 1)),
            1,
            1,
            1,
            disturbance,
            reference,
        )
        error, _ = controller.McFxLMS_SRMSE_ANC(0.05)

        self.assertLess(
            np.mean(error[0, -50:] ** 2),
            np.mean(error[0, :50] ** 2) * 1e-4,
        )


class TestFullyConnectedMultichannel(unittest.TestCase):
    def test_first_sample_axis_and_gradient(self) -> None:
        reference = np.array([[1.0], [-2.0]])
        disturbance = np.array([[0.5], [-1.5]])
        secondary = np.array(
            [
                [[0.2, 0.1], [-0.3, 0.2]],
                [[0.4, -0.2], [0.6, 0.05]],
            ]
        )
        mu = 0.1
        controller = MultiChannelFxLMS(
            3, secondary, 2, reference, disturbance, 2, 2, 2
        )
        controller = controller.McFxLMS_controller(mu)

        np.testing.assert_allclose(controller.Err[:, 0], disturbance[:, 0])
        expected = np.zeros((2, 2))
        for k in range(2):
            for j in range(2):
                expected[k, j] = (
                    mu
                    * reference[j, 0]
                    * np.sum(disturbance[:, 0] * secondary[:, k, 0])
                )
        np.testing.assert_allclose(controller.Wc[:, :, 0], expected)
        np.testing.assert_allclose(controller.Wc[:, :, 1:], 0.0)

    def test_one_by_one_identity_path_converges(self) -> None:
        rng = np.random.default_rng(3)
        reference = rng.standard_normal((1, 500))
        disturbance = reference.copy()
        controller = MultiChannelFxLMS(
            1,
            np.ones((1, 1, 1)),
            1,
            reference,
            disturbance,
            1,
            1,
            1,
        ).McFxLMS_controller(0.05)

        self.assertLess(
            np.mean(controller.Err[0, -50:] ** 2),
            np.mean(controller.Err[0, :50] ** 2) * 1e-4,
        )


@unittest.skipUnless(
    importlib.util.find_spec("torch") is not None,
    "Optional PyTorch is not installed.",
)
class TestTorchParity(unittest.TestCase):
    def test_srms_torch_cpu_matches_numpy_float64(self) -> None:
        from torch_multichannel import TorchMcANC_SRMSE

        rng = np.random.default_rng(20)
        reference = rng.standard_normal(10)
        disturbance = rng.standard_normal((2, 10))
        secondary = rng.standard_normal((2, 2, 3)) * 0.05

        numpy_controller = McANC_SRMSE(
            4, secondary, 3, 2, 2, disturbance, reference
        )
        numpy_error, numpy_controller = (
            numpy_controller.McFxLMS_SRMSE_ANC(5e-4)
        )

        torch_controller = TorchMcANC_SRMSE(
            4,
            secondary,
            3,
            2,
            2,
            disturbance,
            reference,
            device="cpu",
            dtype="float64",
        )
        torch_error, torch_controller = (
            torch_controller.McFxLMS_SRMSE_ANC(5e-4)
        )

        np.testing.assert_allclose(
            torch_error.cpu().numpy(), numpy_error, rtol=1e-13, atol=1e-13
        )
        np.testing.assert_allclose(
            torch_controller.Wc.cpu().numpy(),
            numpy_controller.Wc,
            rtol=1e-13,
            atol=1e-13,
        )

    def test_fully_connected_torch_cpu_matches_numpy_float64(self) -> None:
        from torch_multichannel import TorchMultiChannelFxLMS

        rng = np.random.default_rng(30)
        reference = rng.standard_normal((2, 8))
        disturbance = rng.standard_normal((2, 8))
        secondary = rng.standard_normal((2, 2, 3)) * 0.05

        numpy_controller = MultiChannelFxLMS(
            4, secondary, 3, reference, disturbance, 2, 2, 2
        ).McFxLMS_controller(5e-4)
        torch_controller = TorchMultiChannelFxLMS(
            4,
            secondary,
            3,
            reference,
            disturbance,
            2,
            2,
            2,
            device="cpu",
            dtype="float64",
        ).McFxLMS_controller(5e-4)

        np.testing.assert_allclose(
            torch_controller.Err.cpu().numpy(),
            numpy_controller.Err,
            rtol=1e-13,
            atol=1e-13,
        )
        np.testing.assert_allclose(
            torch_controller.Wc.cpu().numpy(),
            numpy_controller.Wc,
            rtol=1e-13,
            atol=1e-13,
        )


if __name__ == "__main__":
    unittest.main()
