"""Unit tests for Stage 08 overlay coordinate transforms, manifest pairing, and schedd ranking."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import tempfile
import unittest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

import submit_overlay_2track as sov  # noqa: E402


class TestOverlayFrameAndManifest(unittest.TestCase):
    def test_acts_to_faser_orthogonality_det_and_involution(self) -> None:
        # Basis vectors in FASER frame
        ex = sov.faser_to_acts(1.0, 0.0, 0.0)
        ey = sov.faser_to_acts(0.0, 1.0, 0.0)
        ez = sov.faser_to_acts(0.0, 0.0, 1.0)
        self.assertEqual(ex, (0.0, 0.0, -1.0))
        self.assertEqual(ey, (0.0, 1.0, 0.0))
        self.assertEqual(ez, (1.0, 0.0, 0.0))

        # Determinant of R_{F->A} must be +1 (proper rotation, right-handed)
        det = (
            ex[0] * (ey[1] * ez[2] - ey[2] * ez[1])
            - ex[1] * (ey[0] * ez[2] - ey[2] * ez[0])
            + ex[2] * (ey[0] * ez[1] - ey[1] * ez[0])
        )
        self.assertAlmostEqual(det, 1.0, places=12)

        # Round-trip & norm preservation on realistic FASER track kinematics
        sample_vecs = [
            (-34.32, -38.21, 9.94),
            (-31.43, -23.02, 2453.53),
            (112.9, 417.9, 95046.7),
            (4038.2, -6643.5, 2472873.8),
        ]
        for v_f in sample_vecs:
            v_a = sov.faser_to_acts(*v_f)
            v_back = sov.acts_to_faser(*v_a)
            for a, b in zip(v_f, v_back):
                self.assertAlmostEqual(a, b, places=12)
            self.assertAlmostEqual(math.dist((0, 0, 0), v_f), math.dist((0, 0, 0), v_a), places=9)
            # Transverse radius in FASER = hypot(x_F, y_F) = hypot(y_A, z_A)
            self.assertAlmostEqual(math.hypot(v_f[0], v_f[1]), math.hypot(v_a[1], v_a[2]), places=12)

        # Verify against real smoke-test values from 008749 Event 1
        self.assertEqual(sov.acts_to_faser(9.94, -38.21, 34.31), (-34.31, -38.21, 9.94))
        self.assertEqual(sov.acts_to_faser(92983.0, 407.4, -109.7), (109.7, 407.4, 92983.0))

    def test_manifest_322_common_mask_and_323_total(self) -> None:
        with tempfile.TemporaryDirectory() as td:
            root = Path(td)
            overlay_root = root / "overlay"
            out_root = root / "out"
            for spec in sov.OVERLAY_DATASETS:
                ds = spec["name"]
                for sub in ("OverlaidRDOs", "OverlaidNTPs", "OverlaidNTPs_bug"):
                    (overlay_root / ds / sub).mkdir(parents=True)
                for i in range(9):
                    for j in range(9):
                        if ds == "10004X" and (i, j) == (2, 7):
                            continue  # missing RDO + baseline
                        rdo = overlay_root / ds / "OverlaidRDOs" / f"myOverlay_PosBin{i}_NegBin{j}-RDO.pool.root"
                        rdo.write_bytes(b"rdo")
                        if ds == "010698" and (i, j) == (2, 3):
                            continue  # RDO exists, baseline missing
                        for sub in ("OverlaidNTPs", "OverlaidNTPs_bug"):
                            (overlay_root / ds / sub / f"myOverlay_PosBin{i}_NegBin{j}-RDO.pool-ntuple.root").write_bytes(b"ntp")

            tasks, meta = sov.build_overlay_manifest(
                overlay_root,
                out_root,
                cond_tag="OFLCOND-FASER-05",
                geom="TI12MC04",
                geometry_rotation=True,
                material_map=root / "mat.json",
                calypso_source=root / "src",
                calypso_run=root / "run",
                calypso_commit="218813c2",
                overlay_reco_script=root / "reco.py",
                acts_to_faser_bin=root / "conv",
                subset="all",
            )
            self.assertEqual(len(tasks), 323)
            self.assertEqual(meta["total_shards"], 323)
            self.assertEqual(meta["common_baseline_shards"], 322)
            self.assertEqual(meta["extra_non_baseline_shards"], 1)
            extra = [t for t in tasks if not t["in_common_baseline_mask"]]
            self.assertEqual(len(extra), 1)
            self.assertEqual(extra[0]["task_id"], "010698_PosBin2_NegBin3")
            self.assertIsNone(extra[0]["baseline_jul8_ntp"])
            self.assertIsNone(extra[0]["baseline_may30_ntp"])
            self.assertEqual(len({t["task_id"] for t in tasks}), 323)

            # Tier-2 36-shard control subset check
            tasks36, meta36 = sov.build_overlay_manifest(
                overlay_root,
                out_root,
                cond_tag="OFLCOND-FASER-05",
                geom="TI12MC04",
                geometry_rotation=False,
                material_map=root / "mat.json",
                calypso_source=root / "src",
                calypso_run=root / "run",
                calypso_commit="218813c2",
                overlay_reco_script=root / "reco.py",
                acts_to_faser_bin=root / "conv",
                subset="tier2_36",
            )
            self.assertEqual(len(tasks36), 36)
            self.assertEqual(meta36["common_baseline_shards"], 36)

    def test_schedd_ranking_filters_non_ordinary_and_excluded(self) -> None:
        sample = """
batfink02.cern.ch batfink02.cern.ch 10 0 10044
bigbird08.cern.ch bigbird08.cern.ch 672 104925 10053
bigbird20.cern.ch bigbird20.cern.ch 287 23 10055
bigbird21.cern.ch bigbird21.cern.ch 608 235 10055
bigbird101.cern.ch bigbird101.cern.ch 5 0 20373
scbird01.cern.ch scbird01.cern.ch 0 0 2000
bigbird15.cern.ch bigbird15.cern.ch 71 2186 10055
"""
        ranked = sov.rank_bigbird_schedds(sample)
        self.assertEqual([r["name"] for r in ranked], ["bigbird20.cern.ch", "bigbird21.cern.ch", "bigbird15.cern.ch", "bigbird08.cern.ch"])

        ranked_excl = sov.rank_bigbird_schedds(sample, exclude_schedds=["bigbird20.cern.ch"])
        self.assertEqual(ranked_excl[0]["name"], "bigbird21.cern.ch")


if __name__ == "__main__":
    unittest.main()
