"""preCICE 2.x adapter between the MBDyn nodal socket and DUST."""

from __future__ import annotations

import csv
from pathlib import Path
import numpy as np
import precice


WRITE_FIELDS = ("Position", "Velocity", "Rotation", "AngularVelocity")
READ_FIELDS = ("Force", "Moment")
CG_IN = np.array([163.1873833858090, 0.1105295710880, 101.2397973588480])


class CouplingAdapter:
    """Drive the serial-implicit MBDyn/DUST coupling until final time."""

    def __init__(self, mbdyn, config_file: str, diagnostics_file: str | None = None,
                 state_function=None) -> None:
        self.mbdyn = mbdyn
        self.state_function = state_function or (lambda _time: "SMOKE")
        self._diagnostics_handle = None
        self._diagnostics = None
        if diagnostics_file is not None:
            path = Path(diagnostics_file)
            path.parent.mkdir(parents=True, exist_ok=True)
            self._diagnostics_handle = path.open("w", newline="")
            self._diagnostics = csv.writer(self._diagnostics_handle)
            self._diagnostics.writerow([
                "time_s", "state", "coupling_iterations",
                "left_tip_dz_in", "right_tip_dz_in", "symmetric_tip_dz_in",
                "antisymmetric_tip_dz_in", "left_tip_vz_inps", "right_tip_vz_inps",
                "total_fx_lbf", "total_fy_lbf", "total_fz_lbf",
                "global_mx_cg_lbfin", "global_my_cg_lbfin", "global_mz_cg_lbfin",
            ])
        self.is_v3 = hasattr(precice, "Participant")
        if self.is_v3:
            self.interface = precice.Participant("MBDyn", config_file, 0, 1)
            self.mesh_id = "MBDynNodes"
            self.vertex_ids = self.interface.set_mesh_vertices(
                self.mesh_id, self.mbdyn.reference_nodes
            )
            self.data_ids = {name: name for name in (*WRITE_FIELDS, *READ_FIELDS)}
            if self.interface.requires_initial_data():
                self._write_kinematics()
            self.interface.initialize()
            self.precice_dt = self.interface.get_max_time_step_size()
        else:
            self.interface = precice.Interface("MBDyn", config_file, 0, 1)
            self.mesh_id = self.interface.get_mesh_id("MBDynNodes")
            self.vertex_ids = self.interface.set_mesh_vertices(
                self.mesh_id, self.mbdyn.reference_nodes
            )
            self.data_ids = {
                name: self.interface.get_data_id(name, self.mesh_id)
                for name in (*WRITE_FIELDS, *READ_FIELDS)
            }
            self.precice_dt = self.interface.initialize()
            write_initial = precice.action_write_initial_data()
            if self.interface.is_action_required(write_initial):
                self._write_kinematics()
                self.interface.mark_action_fulfilled(write_initial)
            self.interface.initialize_data()

    def _write_kinematics(self) -> None:
        for name in WRITE_FIELDS:
            if self.is_v3:
                self.interface.write_data(
                    self.mesh_id, name, self.vertex_ids, self.mbdyn.data[name]
                )
            else:
                self.interface.write_block_vector_data(
                    self.data_ids[name], self.vertex_ids, self.mbdyn.data[name]
                )

    def _read_loads(self) -> None:
        for name in READ_FIELDS:
            if self.is_v3:
                self.mbdyn.data[name] = self.interface.read_data(
                    self.mesh_id, name, self.vertex_ids, self.precice_dt
                )
            else:
                self.mbdyn.data[name] = self.interface.read_block_vector_data(
                    self.data_ids[name], self.vertex_ids
                )

    def _write_diagnostics(self, time_s: float, iterations: int) -> None:
        if self._diagnostics is None:
            return
        # Structural coupling order: left tip=19, right tip=38 (zero based).
        displacement = self.mbdyn.data["Position"] - self.mbdyn.reference_nodes
        left_dz = float(displacement[19, 2])
        right_dz = float(displacement[38, 2])
        left_vz = float(self.mbdyn.data["Velocity"][19, 2])
        right_vz = float(self.mbdyn.data["Velocity"][38, 2])
        force = np.sum(self.mbdyn.data["Force"], axis=0)
        arms = self.mbdyn.data["Position"] - CG_IN
        moment = np.sum(
            self.mbdyn.data["Moment"] + np.cross(arms, self.mbdyn.data["Force"]),
            axis=0,
        )
        self._diagnostics.writerow([
            f"{time_s:.9f}", self.state_function(time_s), iterations,
            f"{left_dz:.12e}", f"{right_dz:.12e}",
            f"{0.5 * (left_dz + right_dz):.12e}",
            f"{0.5 * (left_dz - right_dz):.12e}",
            f"{left_vz:.12e}", f"{right_vz:.12e}",
            *(f"{value:.12e}" for value in force),
            *(f"{value:.12e}" for value in moment),
        ])
        self._diagnostics_handle.flush()

    def run(self, structural_dt: float) -> None:
        write_checkpoint = None if self.is_v3 else precice.action_write_iteration_checkpoint()
        read_checkpoint = None if self.is_v3 else precice.action_read_iteration_checkpoint()
        saved = None
        simulation_time = 0.0
        iteration = 0

        try:
            while self.interface.is_coupling_ongoing():
                iteration += 1
                write_required = (self.interface.requires_writing_checkpoint()
                                  if self.is_v3 else
                                  self.interface.is_action_required(write_checkpoint))
                if write_required:
                    saved = {
                        name: self.mbdyn.data[name].copy()
                        for name in WRITE_FIELDS
                    }
                    if not self.is_v3:
                        self.interface.mark_action_fulfilled(write_checkpoint)

                # Apply the last DUST loads to the structural trial iteration.
                if self.mbdyn.send_loads(converged=False):
                    break
                if self.mbdyn.receive_kinematics():
                    break

                self._write_kinematics()
                dt = min(structural_dt, self.precice_dt)
                if self.is_v3:
                    self.interface.advance(dt)
                    self.precice_dt = self.interface.get_max_time_step_size()
                else:
                    self.precice_dt = self.interface.advance(dt)
                self._read_loads()

                read_required = (self.interface.requires_reading_checkpoint()
                                 if self.is_v3 else
                                 self.interface.is_action_required(read_checkpoint))
                if read_required:
                    if saved is None:
                        raise RuntimeError("preCICE requested a missing checkpoint")
                    for name, values in saved.items():
                        self.mbdyn.data[name] = values
                    if not self.is_v3:
                        self.interface.mark_action_fulfilled(read_checkpoint)
                else:
                    # Tell MBDyn to accept the converged trial and start a new step.
                    if self.mbdyn.send_loads(converged=True):
                        break
                    if self.mbdyn.receive_kinematics():
                        break
                    simulation_time += dt
                    self._write_diagnostics(simulation_time, iteration)
                    print(
                        f"Coupled time {simulation_time:8.3f} s "
                        f"(last implicit iteration {iteration})",
                        flush=True,
                    )
                    iteration = 0
        finally:
            self.interface.finalize()
            self.mbdyn.close()
            if self._diagnostics_handle is not None:
                self._diagnostics_handle.close()
