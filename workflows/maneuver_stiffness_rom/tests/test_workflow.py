from __future__ import annotations

import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "nastran"))
sys.path.insert(0, str(ROOT.parent / "maneuver_bff"))

from prestress_rom import build, fit_slope, force_text, load_matrix_cases, read_fem_matrix, MatrixCase
from render_mbdyn_case import render
from project_fixed_basis import FemData, project
from prepare_cases import ALTER_NAME, convert_deck, has_param, parameter_names
from maneuver_case import ManeuverPoint, render_point


class PrestressRomTests(unittest.TestCase):
    def test_sol103_statsub_deck_conversion(self) -> None:
        source = ("NASTRAN Q4TAPER=0.7\nSOL 103\nCEND\nTITLE=OLD\n"
                  "SUBCASE 1\n SPC=620\n LOAD=613\nSUBCASE 2\n STATSUB=1\n SPC=620\n METHOD=1\n"
                  " VECTOR(SORT1,REAL,PRINT,PLOT)=901\n"
                  "BEGIN BULK\nPARAM    POST         -1\nENDDATA\n")
        self.assertTrue(has_param(source, "POST"))
        deck = convert_deck(source, 1.0, 1.0)
        self.assertTrue(deck.startswith("NASTRAN Q4TAPER=0.7\n"))
        self.assertIn(f"INCLUDE '{ALTER_NAME}'", deck)
        self.assertIn("STATSUB=1", deck)
        names = parameter_names(deck[deck.index("BEGIN BULK"):])
        self.assertEqual(names.count("POST"), 1)
        self.assertEqual(names.count("FOLLOWK"), 1)
        self.assertIn("PARAM,FOLLOWK,NO", deck)
        self.assertIn("VECTOR(SORT1,REAL,PLOT)=ALL", deck)
        self.assertEqual(deck.count("SPC=620"), 1)
        self.assertIn("prestressed modes must remain free-free", deck)
        self.assertNotIn("QRMETH", deck)
        self.assertNotIn("SOL 106", deck)

    def test_fixed_basis_projection(self) -> None:
        rng = np.random.default_rng(4)
        psi_flat, _ = np.linalg.qr(rng.normal(size=(12, 3)))
        transform = np.asarray([[1.0, 0.2], [-0.1, 0.8], [0.3, -0.2]])
        phi_flat = psi_flat @ transform
        nodes = np.asarray([10, 20])
        baseline_modes = phi_flat.reshape(2, 6, 2).transpose(2, 0, 1)
        snapshot_modes = psi_flat.reshape(2, 6, 3).transpose(2, 0, 1)
        baseline = FemData(
            Path("baseline"), nodes, baseline_modes, np.eye(2), np.eye(2), np.ones((2, 6)),
        )
        stiffness = np.diag([4.0, 9.0, 16.0])
        snapshot = FemData(
            Path("snapshot"), nodes, snapshot_modes, np.eye(3), stiffness, np.ones((2, 6)),
        )
        projected, report = project(baseline, snapshot, [1, 2], set())
        np.testing.assert_allclose(projected, transform.T @ stiffness @ transform, atol=1e-12)
        self.assertLess(report["maximum_weighted_relative_residual"], 1e-12)

    def test_anchored_matrix_fit(self) -> None:
        baseline = np.diag(np.arange(1.0, 7.0))
        slope = np.arange(36.0).reshape(6, 6)
        slope = 0.5 * (slope + slope.T)
        cases = [
            MatrixCase(1.0, baseline, {}, Path("n1"), 0.0),
            MatrixCase(1.3, baseline + 0.3 * slope, {}, Path("n13"), 0.0),
            MatrixCase(1.6, baseline + 0.6 * slope, {}, Path("n16"), 0.0),
        ]
        fitted, _, residuals = fit_slope(cases, 1.0)
        np.testing.assert_allclose(fitted, slope)
        self.assertLess(max(item["relative_residual"] for item in residuals), 1e-14)

    def test_force_sign_and_equilibrium(self) -> None:
        config = json.loads((ROOT / "config.json").read_text())
        text = force_text(config, np.eye(6))
        self.assertIn("-model::drive(PRESTRESS_DELTA_N_DRIVE,Time)", text)
        self.assertIn("(Var-Q_INIT_7)", text)
        self.assertIn("list, 6, 7, 8, 9, 10, 11, 12", text)
        self.assertNotIn("\",,", text)

    def test_rejects_wrong_basis(self) -> None:
        base_config = json.loads((ROOT / "config.json").read_text())
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            fem = temp_path / "dummy.fem"
            fem.write_text("dummy")
            digest = hashlib.sha256(b"dummy").hexdigest()
            config = {
                **base_config,
                "model": {**base_config["model"], "baseline_fem_sha256": digest},
                "matrix_cases": [{"load_factor": 1.0, "directory": "case"}],
            }
            case = temp_path / "case"
            case.mkdir()
            np.savetxt(case / "kh_fixed_basis.csv", np.eye(6), delimiter=",")
            (case / "metadata.json").write_text(json.dumps({
                "load_factor": 1.0,
                "basis_id": "wrong",
                "modes": config["model"]["selected_modes"],
                "fem_sha256": digest,
                "units": "unit",
                "projection": "Phi0.T @ K @ Phi0",
            }))
            with self.assertRaisesRegex(ValueError, "basis_id"):
                load_matrix_cases(config, temp_path)

    def test_build_and_render_end_to_end(self) -> None:
        config = json.loads((ROOT / "config.json").read_text())
        fem = (ROOT / config["model"]["baseline_fem"]).resolve()
        full_k = read_fem_matrix(fem, 10)
        index = np.asarray([mode - 1 for mode in config["model"]["selected_modes"]])
        k1 = full_k[np.ix_(index, index)]
        slope = np.diag(np.linspace(-10.0, 15.0, 6))
        with tempfile.TemporaryDirectory() as temp:
            temp_path = Path(temp)
            for tag, load_factor, matrix in (
                ("n1", 1.0, k1),
                ("n16", 1.6, k1 + 0.6 * slope),
            ):
                case = temp_path / tag
                case.mkdir()
                np.savetxt(case / "kh_fixed_basis.csv", matrix, delimiter=",")
                (case / "metadata.json").write_text(json.dumps({
                    "load_factor": load_factor,
                    "basis_id": config["model"]["basis_id"],
                    "modes": config["model"]["selected_modes"],
                    "fem_sha256": config["model"]["baseline_fem_sha256"],
                    "units": "lbf/in modal generalized stiffness",
                    "projection": "Phi0.T @ K_tangent @ Phi0",
                }))
            config["model"]["baseline_fem"] = str(fem)
            config["matrix_cases"] = [
                {"load_factor": 1.0, "directory": str(temp_path / "n1")},
                {"load_factor": 1.6, "directory": str(temp_path / "n16")},
            ]
            config_path = temp_path / "config.json"
            config_path.write_text(json.dumps(config))
            generated = temp_path / "generated"
            report = build(config_path, generated)
            np.testing.assert_allclose(report["k_h_n"], slope, atol=1e-11)
            self.assertFalse(report["linearity_validated"])
            base_case = temp_path / "base.mbd"
            base_case.write_text(
                '# MANEUVER_METADATA {"family": "pullup"}\n'
                "begin: initial value;\n"
                "    forces: 2;\n"
                "    force: NASTRAN_DLM_ROM_FORCE, modal, MODAL_JOINT;\n"
                "    # Idealized active modal damper\n"
            )
            rendered = temp_path / "rendered.mbd"
            render(base_case, rendered, config_path, "commanded", generated)
            text = rendered.read_text()
            self.assertIn("forces: 3;", text)
            self.assertIn("force: PRESTRESS_ROM_FORCE, modal", text)
            self.assertNotIn("\",,", text)

    def test_dive_pullup_uses_exogenous_pitch_rate_load_schedule(self) -> None:
        point = ManeuverPoint(
            "dive_pullup", 66.75, 1.0, excited=False,
            pitch_angle_deg=10.22625, pitch_rate_deg_s=5.05,
            nominal_load_factor=1.6,
        )
        with tempfile.TemporaryDirectory() as temp:
            source = Path(temp) / "source.mbd"
            output = Path(temp) / "output.mbd"
            source.write_text(render_point(point))
            render(
                source, output, ROOT / "config.json", "maneuver_pitch_rate",
                ROOT / "generated",
            )
            text = output.read_text()
            self.assertIn("PRESTRESS_ROM_INJECTED", text)
            self.assertIn('"nz_source": "maneuver_pitch_rate"', text)
            self.assertIn("forces: 3;", text)
            self.assertIn("force: PRESTRESS_ROM_FORCE, modal, MODAL_JOINT", text)
            self.assertIn(
                "DIVE_PULLUP_ENABLE*VINF*model::drive(MANEUVER_Q_COMMAND_DRIVE,Time)/GRAVITY",
                text,
            )
            self.assertIn("((Time<SAS_OFF_START)||(Time>=SAS_ON_START))", text)

            wrong = Path(temp) / "wrong.mbd"
            wrong.write_text(render_point(ManeuverPoint("pullup", 66.75, 1.6, excited=False)))
            with self.assertRaisesRegex(ValueError, "family=dive_pullup"):
                render(
                    wrong, output, ROOT / "config.json", "maneuver_pitch_rate",
                    ROOT / "generated",
                )

if __name__ == "__main__":
    unittest.main()
