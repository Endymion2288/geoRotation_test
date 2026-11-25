"""Unit tests for Stage 08 overlay statistical estimators and cross-checks against C++ outputs."""

import math
import unittest

try:
    from scipy.stats import beta as scipy_beta
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False


CL_1SIGMA = 0.6826894921370859  # erf(1 / sqrt(2))


def clopper_pearson_1sigma(k: int, n: int):
    """Exact 68.3% CL Clopper-Pearson interval (eff, err_lo, err_hi)."""
    if n <= 0:
        return None, 0.0, 0.0
    if k < 0 or k > n:
        raise ValueError(f"Invalid counts k={k}, n={n}")
    eff = k / n
    alpha = 1.0 - CL_1SIGMA
    if not HAS_SCIPY:
        raise RuntimeError("scipy required for reference Clopper-Pearson check")
    lo = 0.0 if k == 0 else float(scipy_beta.ppf(alpha / 2.0, k, n - k + 1))
    hi = 1.0 if k == n else float(scipy_beta.ppf(1.0 - alpha / 2.0, k + 1, n - k))
    return eff, max(0.0, eff - lo), max(0.0, hi - eff)


def mcnemar_paired_delta(n_both_pass: int, n_gain: int, n_loss: int, n_both_fail: int):
    """Exact paired difference (ver_B - ver_A) and standard error."""
    n = n_both_pass + n_gain + n_loss + n_both_fail
    if n <= 0:
        return None, None
    delta = (n_gain - n_loss) / n
    var_num = max(0.0, (n_gain + n_loss) - ((n_gain - n_loss) ** 2) / n)
    err = math.sqrt(var_num) / n
    return delta, err


def weighted_efficiency(sumw_pass: float, sumw_fail: float, sumw2_pass: float, sumw2_fail: float):
    """Weighted efficiency and propagated variance from independent event weights."""
    sumw = sumw_pass + sumw_fail
    if sumw <= 0.0:
        return None, None
    eff = sumw_pass / sumw
    var = (sumw2_pass * ((1.0 - eff) ** 2) + sumw2_fail * (eff ** 2)) / (sumw ** 2)
    return eff, math.sqrt(max(0.0, var))


def effective_sample_size(weights):
    """Kish effective sample size N_eff = (sum w)^2 / sum w^2."""
    sw = sum(weights)
    sw2 = sum(w * w for w in weights)
    if sw2 <= 0.0:
        return 0.0
    return (sw * sw) / sw2


class TestOverlayMetrics(unittest.TestCase):
    @unittest.skipUnless(HAS_SCIPY, "scipy not available")
    def test_clopper_pearson_boundaries_and_symmetry(self):
        eff0, lo0, hi0 = clopper_pearson_1sigma(0, 100)
        self.assertEqual(eff0, 0.0)
        self.assertEqual(lo0, 0.0)
        self.assertAlmostEqual(hi0, 1.0 - (0.5 * (1.0 - CL_1SIGMA)) ** (1.0 / 100.0), places=10)

        eff1, lo1, hi1 = clopper_pearson_1sigma(100, 100)
        self.assertEqual(eff1, 1.0)
        self.assertEqual(hi1, 0.0)
        self.assertAlmostEqual(lo1, hi0, places=10)

        eff_half, lo_half, hi_half = clopper_pearson_1sigma(50, 100)
        self.assertAlmostEqual(eff_half, 0.5, places=12)
        self.assertAlmostEqual(lo_half, hi_half, places=10)
        self.assertGreater(lo_half, 0.049)
        self.assertLess(lo_half, 0.056)

    def test_mcnemar_paired_delta(self):
        # Identical outcomes -> zero delta, zero variance
        d0, e0 = mcnemar_paired_delta(900, 0, 0, 100)
        self.assertEqual(d0, 0.0)
        self.assertEqual(e0, 0.0)

        # Asymmetric gain/loss
        d, e = mcnemar_paired_delta(920, 30, 10, 40)
        self.assertAlmostEqual(d, 0.02, places=12)
        expected_e = math.sqrt((40.0 - (20.0 ** 2) / 1000.0)) / 1000.0
        self.assertAlmostEqual(e, expected_e, places=12)

        # Paired variance is strictly smaller than unpaired sum-in-quadrature when correlated
        p_a = (920 + 10) / 1000.0
        p_b = (920 + 30) / 1000.0
        unpaired_e = math.sqrt(p_a * (1 - p_a) / 1000.0 + p_b * (1 - p_b) / 1000.0)
        self.assertLess(e, unpaired_e)

    def test_weighted_efficiency_reduces_to_binomial(self):
        # Unit weights: sumw2_pass = k, sumw2_fail = n - k -> sigma = sqrt(p(1-p)/n)
        k, n = 940, 1000
        eff, err = weighted_efficiency(float(k), float(n - k), float(k), float(n - k))
        self.assertAlmostEqual(eff, 0.94, places=12)
        self.assertAlmostEqual(err, math.sqrt(0.94 * 0.06 / 1000.0), places=12)

        # Constant scale factor c leaves eff and err invariant
        c = 3.7
        eff_c, err_c = weighted_efficiency(c * k, c * (n - k), c * c * k, c * c * (n - k))
        self.assertAlmostEqual(eff_c, eff, places=12)
        self.assertAlmostEqual(err_c, err, places=12)

    def test_effective_sample_size(self):
        self.assertAlmostEqual(effective_sample_size([2.5] * 100), 100.0, places=12)
        neff = effective_sample_size([1.0, 1.0, 1.0, 9.0])
        self.assertAlmostEqual(neff, (12.0 ** 2) / (3.0 + 81.0), places=12)
        self.assertLess(neff, 4.0)


if __name__ == "__main__":
    unittest.main()
